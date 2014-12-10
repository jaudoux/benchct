#! /usr/bin/perl

#author Nicolas Philippe
#goal   post treatment of chimeras with stricter CRAC parameters
#       


use strict;
use warnings;
use POSIX;

my $debug = 0;

use constant BREAK_MIN_THRESHOLD => 0.75;
use constant IS_STRANDED => 'FALSE';
use constant MIN_SUPPORT_IN_CHIMERA => 2;
use constant MAX_SINGLE_WINDOWS => 5;
use constant MAX_LOCATION_OUTSIDE => 20;
use constant MAX_SUPPORT_IN_CHIMERA => 500;
use constant MAX_VARIATION_IN_CHIMERA => 0.95;


if (@ARGV < 2){
    print STDERR "Usage: postProcessChimera.pl <file.chimera> <k-mer length> [--min-break chimera_break_min_length] [--stranded is_stranded] [--min-support min_support_in_chimera] [--max-support max_support_in_chimera] [--max-single max_single_windows] [--max-variation max_variation_in_chimera] [--max-location max_location_outside] [--output-reads file_out flag]

                  flags description for the filtration:
                      1: If the chimera has a tipical profil, i.e locations on two different chr.
                      2: If the two locations are both on same chr and strand but the distance is too large for a splice.  
                      4: If the two locations are both on same chr and strand but not in the good order. 
                      8: If the chimera is inverted like 4 but with a repeated factor in the read. 
                     16: If the two locations are on same chr but on different strand.  
                     32: If there are no single location before and after the chimera break, 
                         which creates ambuiguity (extension until max_single_windows)  
                     64: If the chimera break is not clean (it is segmented by falses locations).
                    128: If the chimera break is not long enough (< break_min_length).
                    256: If a k-mer is not enough supported (< min_support_in_chimera).
                    512: If a k-mer is too supported (> max_support_in_chimera).
                   1024: there is too much variation inside chimera break. 
                   2048: there are too much locations before or after the chimera break. 
                   4096: the reads that contain chimera junction are all the same
             
                  For example, a chimera flag of 0 means that the chimera is the once at this position 
                  and it has a fairly long clean break with its two locations on two differents chr.

                  chimera_break_min_length: the break length from which the break is not long enough compared to the theoretical k-1 
                                   (by default k-mer length - 1 - ".BREAK_MIN_THRESHOLD.").
                  is_stranded: TRUE if the reads are stranded otherwise FALSE (by default ".IS_STRANDED.").
                  min_support_in_chimera: the minimal support for a k-mer in the chimera break (by default ".MIN_SUPPORT_IN_CHIMERA.").
                  max_support_in_chimera: the maximal support considered for a k-mer in the chimera break 
                                          (by default ".MAX_SUPPORT_IN_CHIMERA.").
                  max_single_windows: the windows length to search a single location before 
                                      and after the chimera break (by default ".MAX_SINGLE_WINDOWS.").
                  max_variation_in_chimera: the max percent variation between max_support_in_chimera et min_support_in_chimera
                                            (by default".MAX_VARIATION_IN_CHIMERA."). 
                  max_location_outside: the max location level after and before the chimera break (by default".MAX_LOCATION_OUTSIDE."). 
                  
                  print all chimera features on STDOUT with the following format: 
                      read_id\tflag_filter\tnb_occ_chimera\tchr1\tloc1\tstrand1\tchr2\tloc2\tstrand2
                  print statistics on STDERR
     
                  if --output-reads is activated, the corresponding reads of CRAC associated to the chimera 
                  with the flags description chosen are saved in the file. 
\n";
    exit 1;
}




my $break_min_length = int($ARGV[1] * BREAK_MIN_THRESHOLD);
my $isStranded = IS_STRANDED;
my $min_support_in_chimera = MIN_SUPPORT_IN_CHIMERA;
my $max_support_in_chimera = MAX_SUPPORT_IN_CHIMERA;
my $max_single_windows = MAX_SINGLE_WINDOWS;
my $max_variation_in_chimera = MAX_VARIATION_IN_CHIMERA;
my $max_location_outside = MAX_LOCATION_OUTSIDE;
my ($file_reads, $flag);


for (my $i=2 ; $i <= @ARGV ; $i+=2){
    if (defined $ARGV[$i]){
	if ($ARGV[$i] eq '--min-break'){
	    $break_min_length = $ARGV[$i+1];    
	}elsif ($ARGV[$i] eq '--stranded'){
	    $isStranded = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--min-support'){
	    $min_support_in_chimera = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--max-support'){
	    $max_support_in_chimera = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--max-single'){
	    $max_single_windows = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--max-variation'){
	    $max_variation_in_chimera = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--max-location'){
	    $max_location_outside = $ARGV[$i+1];
	}elsif ($ARGV[$i] eq '--output-reads'){
	    $file_reads = $ARGV[$i+1];
	    if (!defined $ARGV[$i+2]){
		print STDERR "you must specified flag after --output-reads $file_reads\n";
		exit(2);
	    }else{
		$flag = $ARGV[$i+2];
	    }
	    $i++;
	}else{
	    print STDERR "parameter $ARGV[$i] does not exist!\n";
	    exit(1);
	}
    }
}

if (defined $file_reads){
    open(OUT, ">$file_reads") or die("unabled to open $file_reads: $!");
}


sub reverse_complement($){
    my $dna = shift;
    
    # reverse the DNA sequence
    my $revcomp = reverse($dna);
    
    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub reverse_tab($) {
    my $string = shift;
    my @tab = split(/,/,$string);
    my $newString;
    for (my $i=$#tab ; $i > 0 ; $i--){
	$newString .= $tab[$i];
	$newString .= ",";
    }
    $newString .= $tab[0];
    return $newString;
}


sub reverseRead($$$$$$$$$$$){
    my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq,$read_length) = @_;
    my $tmp = $chr1;
    $chr1 = $chr2;
    $chr2 = $tmp;
    $tmp = $pos1;
    $pos1 = $pos2;
    $pos2 = $tmp;
    $tmp = -$strand2;
    $strand2 = -$strand1;
    $strand1 = $tmp;
    $tmp = $read_length - $posMut;
    $posMut = $tmp;
    $tmp = reverse_tab($profilSupport);
    $profilSupport = $tmp;
    $tmp = reverse_tab($profilLoc);
    $profilLoc = $tmp;
    $tmp = reverse_complement($seq);
    $seq = $tmp;
    return ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq)
}

sub getFeaturesFromBreak($$$$$$$){
    my ($posMut,$profilSupport,$profilLoc,$klength,$max_single_windows,$break_min_length,$read_length) = @_;
    my @locs = split(/,/,$profilLoc);
    my @support = split(/,/,$profilSupport);
    my $nb_break = 0;
    my $found = 0;
    my $minLength = 0;
    my $minSupport = $max_support_in_chimera;
    my $maxSupport = 0;
    my $averageSupport = 0;
    my $maxLocation = 1;
    my $tmp = 0;
    my $break_l = 1;
    my $foundSingleStart = 0;
    my $foundSingleEnd = 0;
    my $posStart = ($posMut - $klength);
    if ($posStart < 0){
	$posStart = 0 
    }

    # we adjust the chimera break start
    while ($locs[$posStart] > 0){
	$posStart++;
    }
    # we save the start of the break 
    my $posStartBreak = $posStart;
    # we extend according to the windows length
    $posStart -= $max_single_windows;
    if ($posStart < 0){
	$posStart = 0 
    }
    # posMut is the end of the break
    my $posEnd = $posMut + $max_single_windows;
    if ($posEnd > ($read_length - $klength)){
	$posEnd = $read_length - $klength; 
    }
      
    for (my $i=$posStart ; $i<= $posEnd ; $i++){
	# start of a new break
	if (($locs[$i] == 0) && !$found && $i >= $posStartBreak && $i <= $posMut){
	    if ($support[$i] < $minSupport){
		$minSupport = $support[$i];
	    }
	    if ($support[$i] > $maxSupport){
		$maxSupport = $support[$i];
	    }
	    $averageSupport = $support[$i];
	    $tmp = 1;
	    $found = 1;
	    $nb_break += 1;
	# Inside the chimera break but there is a location and if tmp > 0, we get out the precedent 
        # break and compute its length to know if it is the shortest
	}elsif (($locs[$i] != 0) && $i >= $posStartBreak && $i <= $posMut){
	    if ($tmp > $minLength && $tmp > 0){
		$minLength = $tmp;
	    }
	    $tmp = 0;
	    $found = 0;
	# inside the break, we save the minimal support from a k-mer
	}elsif($found && $i >= $posStartBreak && $i <= $posMut){
	    if ($support[$i] < $minSupport){
		$minSupport = $support[$i];
	    }
	    if ($support[$i] > $maxSupport){
		$maxSupport = $support[$i];
	    }
	    $averageSupport += $support[$i];
	    $tmp += 1;
	    $break_l += 1;
	}
	
        #check the max location
	if ($locs[$i] > $maxLocation){
	    $maxLocation = $locs[$i];
	}
	
	# search if there is a single loc before the chimera break    
	if (($locs[$i] == 1) && $i <= ($posStart + $max_single_windows) && !$foundSingleStart){
	    $foundSingleStart = 1;
	# search if there is a single loc after the chimera break        
	}elsif (($locs[$i] == 1) && $i > $posMut && !$foundSingleEnd){
	    $foundSingleEnd = 1;
	}
	
    }

    # compute average support inside chimera break
    $averageSupport = floor($averageSupport/$break_l);
    
    # check if the last break length is shorter than the former
    if ($tmp > $minLength && $tmp > 0){
	$minLength = $tmp;
    }

    my $foundSingle;
    # iff single on the left and on the right of chimera break
    if ($foundSingleStart && $foundSingleEnd){
	$foundSingle = 1;
    }else{
	$foundSingle = 0;
    }

    return ($nb_break,$minLength,$minSupport,$maxSupport,$averageSupport,$foundSingle,$maxLocation);
}

########################################################
open(IN,$ARGV[0]) or die("impossible to open $ARGV[0] : $!");
my %hash;
my $read_length;

# First Step: a classical filtration
while (<IN>){
    chomp;
    my ($read_id,$type,$chr1,$strand1,$pos1,$chr2,$strand2,$pos2,$posMut,$singleLoc,$seq,$profilSupport,$profilLoc) = $_ =~ /^(\d+)\s+(\S+)\s+(\S+)\|(\S+)?,(\S+)?\s+(\S+)\|(\S+)?,(\S+)\s+pos_junction=(\S+)\s+(.*?)\s+([CATG]+)\s+(\S+)\s+(\S+)$/ or next;
    $read_length = length($seq);
      
    if ($debug){
	print STDERR "$_\n";
    }
    
    # is the reads are not stranded, the antisens reads and the sens reads 
    # are considered as the same.
    if ($isStranded eq 'FALSE'){
	if (!defined($hash{"$chr1,$pos1,$chr2,$pos2"}) && !defined($hash{"$chr2,$pos2,$chr1,$pos1"})){
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} = 0;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'MSUPPORT'} = 1;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ANCHOR'} = 1;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} = $posMut;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ID'} = $read_id;
	}elsif (!defined($hash{"$chr1,$pos1,$chr2,$pos2"}) && defined($hash{"$chr2,$pos2,$chr1,$pos1"})){
	    ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq) = reverseRead($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq,$read_length);
	}
    }else{
	# when RNA-Seq is oriented, the second paired-read give the orientation
	if ($read_id % 2 == 0){
	    ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq) = reverseRead($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$posMut,$profilSupport,$profilLoc,$seq,$read_length);
	}
	if (!defined($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"})){
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} = 0;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'MSUPPORT'} = 1;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ANCHOR'} = 1;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} = $posMut;
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ID'} = $read_id;
	}
    }

    # if several reads cover the chimera, we increase the counter for each read 
    unless (($read_id % 2 == 0 && $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ID'} == ($read_id+1))
	    || ($read_id % 2 == 1 && $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ID'} == ($read_id-1))){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} += 1;
    }

    # if the chimera junction is not the same on the reads, that is to say reads are not anchored
    if ($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} ne $posMut){ 
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ANCHOR'} = 0;
    }
    
    #save the old value
    my $oldValue = 0;   
    if (defined($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'})){
	$oldValue = $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'};
    }
    
    # if the chimera is already define, we take it if its break is more in the middle.
    if ($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} == 1 || ($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} > $posMut && $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} < ($read_length - $ARGV[1] - $posMut))){ 
	
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'POSMUT'} = $posMut;
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ID'} = $read_id;
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 0;
	if (defined $file_reads){
	    my $chimera = "$read_id $type $chr1\|$strand1,$pos1 $chr2\|$strand2,$pos2 pos_junction=$posMut $singleLoc $seq $profilSupport $profilLoc"; 
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'READS'} = $chimera;
	}
    }

    
    if ($chr1 eq $chr2){
	if ($strand1 == $strand2){
	    if (($strand1 == 1 && $pos1 > $pos2) || ($strand1 == -1 && $pos2 > $pos1)){
		if (abs($pos2 - $pos1) >= $read_length){
		    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 4;
		}else{
		    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 8;
		}
	    }else{
		$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 2;
	    }
	}else{
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 16;
	}
    }else{
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = 1;
    }
    
    my ($nb_break,$minLength,$supportMin, $supportMax, $averageSupport, $foundSingle, $maxLocation) 
	= getFeaturesFromBreak($posMut,$profilSupport,$profilLoc,$ARGV[1],$max_single_windows, $break_min_length, $read_length);
    my $current_variation = 1 - ($supportMin/$supportMax);
    
    #update the chimera support 
    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'MSUPPORT'} = $averageSupport;
    
    if ($debug){	
	print STDERR "\nnb_break: $nb_break, min_break: $minLength, min_support: $supportMin, max_support: $supportMax, average_support: $averageSupport, hasSingle: $foundSingle, variation_in_support: $current_variation, max_location: $maxLocation\n-------\n";
    }
    
    if (!$foundSingle){
	# print STDERR "\nnb_break: $nb_break, min_break: $minLength, min_support: $supportMin, max_support: $supportMax, average_support: $averageSupport, hasSingle: $foundSingle, variation_in_support: $current_variation, max_location: $maxLocation\n-------\n";
	# print STDERR "$_\n";
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 32;
    }
    if ($nb_break > 1){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 64;
    }
    if ($minLength < $break_min_length){
	# print STDERR "\nnb_break: $nb_break, min_break: $minLength, min_support: $supportMin, max_support: $supportMax, average_support: $averageSupport, hasSingle: $foundSingle, variation_in_support: $current_variation, max_location: $maxLocation\n-------\n";
	# print STDERR "$_\n";
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 128;
    }
    if ($supportMax > $max_support_in_chimera){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 512;
    }
    if ($max_variation_in_chimera < $current_variation){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 1024;
    }
    if ($maxLocation > $max_location_outside){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 2048;
    }
    
    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'BREAK'} = $minLength;    
    
    #check the min between oldValue and newValue
    if ($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} > $oldValue && ($oldValue != 0) ){
	$hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} = $oldValue;
    }
}
close(IN);

#Second step: filtration of multiple chimeras on the same location or the ANCHOR cases 
my %hashFilterMulti;
foreach my $key (keys %hash){
    my @elt = keys %{$hash{$key}};
    my ($chr1,$pos1,$chr2,$pos2) = split(/,/,$key);
    foreach my $elt (@elt){
	my ($strand1,$strand2) = split(/,/,$elt);
        
	# check nb_reads cover the chimera junction
	if ($hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} < $min_support_in_chimera){
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} += 256;
	}
	
        # check the anchor 
	if ( $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'ANCHOR'} == 1 && $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'OCCNB'} > 1){
	    $hash{"$chr1,$pos1,$chr2,$pos2"}{"$strand1,$strand2"}{'VALUE'} +=  4096;
	}
    }
}

print "read_id\tflag_filter\tspanning_reads\tchr1\tloc1\tstrand1\tchr2\tloc2\tstrand2\n";
my ($nb_chimeras,$nb_reads) = (0,0);
foreach my $key (keys %hash){
    my @elt = keys %{$hash{$key}};
    my @allkeys = split(/,/,$key);
    foreach my $elt (@elt){
	if (defined $file_reads){
	    my $isGoodFlag = ($hash{$key}{$elt}{'VALUE'} & $flag);
	    if ($isGoodFlag && $hash{$key}{$elt}{'VALUE'} <= $flag){
		$nb_chimeras++;
		$nb_reads += $hash{$key}{$elt}{'OCCNB'};	
		print OUT $hash{$key}{$elt}{'READS'}."\t"."spanning_junction=".$hash{$key}{$elt}{'MSUPPORT'}."\t"."spanning_reads=".$hash{$key}{$elt}{'OCCNB'}."\t"."chimera_type=".$hash{$key}{$elt}{'VALUE'}."\t"."break_length=".$hash{$key}{$elt}{'BREAK'}."\n";
	    }
	}
	my @allelts = split(/,/,$elt);
	print $hash{$key}{$elt}{'ID'}."\t".$hash{$key}{$elt}{'VALUE'}."\t".$hash{$key}{$elt}{'OCCNB'}."\t".$allkeys[0]."\t".$allkeys[1]."\t".$allelts[0]."\t".$allkeys[2]."\t".$allkeys[3]."\t".$allelts[1]."\n";
    }
}
close(OUT);
if ($debug){
    if ($isStranded eq 'TRUE'){
	print STDERR "RNA-seq is stranded:\n";
    }else{
	print STDERR "RNA-seq is not stranded:\n";
    } 
    print STDERR "\tchimeras: $nb_chimeras\tnb_reads: $nb_reads\n";
}
