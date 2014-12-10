#!/usr/bin/perl -w
use strict;

#use constant FIXED_THRESHOLD => 100;
use constant VARIANT_THRESHOLD => 100;
use constant END_EXON => 0;
use constant START_EXON => 1;

my $stranded = 0;

if (@ARGV  < 3) {
    print STDERR "Usage: $0 [--stranded] <gff file> <k-mer size> <junction upper bound> <[crac,topHat] result>+\n";
    print STDERR "result on STDOUT\n"; 
    exit 1;
}



# Binary search, array passed by reference

# search array of integers a for given integer x
# return index where found or -1 if not found
sub bsearch {
    my ($x, $a, $idx) = @_;            # search for x in array a (at index idx in subarray)
    my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
    my $i;                       # index of probe
    while ($l <= $u) {
	$i = int(($l + $u)/2);
	if ($a->[$i]->[$idx] < $x) {
	    $l = $i+1;
	}
	elsif ($a->[$i]->[$idx] > $x) {
	    $u = $i-1;
	} 
	else {
	    return $i; # found
	}
    }
    # Not found, get the closest
    if ($i > 0 && abs($x-$a->[$i]->[$idx]) > abs($x - $a->[$i-1]->[$idx])) {
	return $i-1;
    } elsif ($i < @$a-1 && abs($x-$a->[$i]->[$idx]) > abs($x - $a->[$i+1]->[$idx])) {
	return $i+1;
    } else {
	return $i;
    }
}

my %hashloc;

sub noMultiplicity($$$$) {
    my ($chr1, $exon1, $chr2, $exon2) = @_;
#     my @keys = sort {$a <=> $b} keys %{$hashloc{"$chr1,$chr2"}};
#     if ($#keys < 0) {
# 	return 1;
#     }
#     my $i = bsearch($exon1, \@keys);
#     for (my $delta = -1; $delta <= 1; $delta++) {
# 	if (($delta+$i >= 0 && $delta+$i <= $#keys)
# 	    && ((($keys[$delta+$i] <= $exon1)
# 		 && $exon1 <= $hashloc{"$chr1,$chr2"}->{$keys[$delta+$i]})
# 		||
# 		(($keys[$delta+$i] <= $exon2)
# 		 && $exon2 <= $hashloc{"$chr1,$chr2"}->{$keys[$delta+$i]}))
# 	    ) {
# 	    return 0;
# 	}
#     }
#     return 1;
   return ! defined($hashloc{"$chr1,$exon1,$chr2,$exon2"});
}

sub recordJunction($$$$) {
    my ($chr1, $exon1, $chr2, $exon2) = @_;
    $hashloc{"$chr1,$exon1,$chr2,$exon2"}=1;
#     $hashloc{"$chr1,$chr2"}->{$exon1} = $exon2;
}

# In fact not gene but transcript
sub findSameGene($$) {
    my ($start, $end) = @_;
    my %elementsStart;
    my %elementsEnd;
    @elementsStart{@{$start->{GENE}}} = ();
    if (ref($end) eq 'HASH') {
	# Parameter is an exon
	@elementsEnd{@{$end->{GENE}}} = ();
    } elsif (ref($end) eq 'ARRAY') {
	# Parameter is an array of genes
	@elementsEnd{@$end} = ();
    } else {
	# What is this parameter??
	print STDERR 'Wrong parameter given as $end',"\n";
	exit 2;
    }
    my %elements;
    my @inter;
    foreach my $elem (keys %elementsStart, keys %elementsEnd) {
	$elements{$elem}++;
	if ($elements{$elem} == 2) {
	    push @inter, $elem;
	}
    }
    return @inter;    
}

sub crossWithRefSeq{
    my ($expressedGenes, $knownExons, $startExons, $infoExons, $tagNumber, $loc, $loc_end, $chr, $chr_end, $strand, $strand_end, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,$bound_max, $recursive_call) = @_;
    my %knownExons = %$knownExons; 
    my %startExons = %$startExons; 

    my @exons;
    $exons[END_EXON] = \%$knownExons;
    $exons[START_EXON] = \%startExons;
    

    my $tmp;
#     if ($loc > $loc_end) {
# 	$tmp = $loc;
# 	$loc = $loc_end;
# 	$loc_end = $tmp;
	
# 	$tmp = $chr;
# 	$chr = $chr_end;
# 	$chr_end = $tmp;
#     }

    if (! defined $recursive_call) {
	$recursive_call = 0;
    }
    
    if (! defined(@{$knownExons{$chr}}) && defined(@{$knownExons{'chr'.$chr}})) {
	$chr = 'chr'.$chr;
	$chr_end = 'chr'.$chr_end;
    }
#     my %expressedGenes = %$expressedGenes;
    # new junction imply new transcript 
    if (! defined(@{$knownExons{$chr}}) || ! defined(@{$knownExons{$chr_end}})) {
	print "newTranscript [$chr:$loc;$chr_end:$loc_end] $tagNumber\n";
	if(noMultiplicity($chr,$loc,$chr_end,$loc_end)){
	    $newTranscript++;
	    recordJunction($chr,$loc,$chr_end,$loc_end);
	}
    } else {
	my ($indexStart, $indexEnd);
	my $isChimera = 0;
	my $nbNiceExonsMet = 0;
	my $border_first = END_EXON;
	my $border_last;

	# We don't where the exons are supposed to be so we need to search in several different ways
	# The most likely is that both exons are on the same strand so one of the position should be
	# at the end of an exon and the other at the start of an exon.
	# However we also allow for two positions corresponding to the start (resp. to the end) of an exon.
	# But in that case that is not a ``classic'' splice.
	$indexStart = bsearch($loc, \@{$exons[END_EXON]->{$chr}}, END_EXON);
	if (abs($loc - $knownExons{$chr}[$indexStart]->[END_EXON]) > $bound_max) {
	    $indexStart = bsearch($loc, \@{$exons[START_EXON]->{$chr}}, START_EXON);
	    $border_first = START_EXON;
	}
# 	if (abs($loc - $exons[$border_first]->{$chr}[$indexStart]->[$border_first]) <= $bound_max) {
	    # We have found a nice candidate so far.
	    # Try to see if we can find a nice candidate with loc_end
	    $border_last = ($border_first + 1) % 2;
	    $indexEnd = bsearch($loc_end, \@{$exons[$border_last]->{$chr_end}}, $border_last);
	    
	    if (abs($loc_end - $exons[$border_last]->{$chr_end}[$indexEnd]->[$border_last]) 
		> $bound_max) {
		$border_last = ($border_last + 1) % 2;
		$indexEnd = bsearch($loc_end, \@{$exons[$border_last]->{$chr_end}}, $border_last);
	    }
# 	}

	my $exonStart = $exons[$border_first]->{$chr}[$indexStart]->[$border_first];
	my $exonEnd = $exons[$border_last]->{$chr_end}[$indexEnd]->[$border_last];
	my @current_gene;


	if ($chr eq $chr_end  && $border_last != $border_first
	    && ($strand == $strand_end  
		&& (! $stranded || (($strand == 1 && $border_first == END_EXON)
				    || ($strand == -1 && $border_last == END_EXON))))
	    && 	(@current_gene = findSameGene(\%{$infoExons->{$chr.'@'.$exonStart}},
					     \%{$infoExons->{$chr_end.'@'.$exonEnd}}))) {
	    #  We have a quite normal splice
	    # search how many exons separate both
	    my $lastExon;
	    my $firstBorder = END_EXON;
	    my $lastBorder = START_EXON;
	    my $startingExon;
	    if ($border_last == START_EXON) {
		# We take the index corresponding to the start exon
		$lastExon = $indexEnd;
		$startingExon = $indexStart;
	    } else {
		$lastExon = $indexStart;
		$startingExon = $indexEnd;
	    }
# 	    print "loc = $loc, loc_end = $loc_end\n";

	    my $currentExon = $lastExon;
	    do {
		$currentExon--;
		# print "lastExon = $lastExon -> pos end = ".$exons[$lastBorder]->{$chr}[$lastExon]->[END_EXON]." pos start = ".$exons[$lastBorder]->{$chr}[$lastExon]->[START_EXON]."\n";
		# print "currentExon = $currentExon -> pos end = ".$exons[$lastBorder]->{$chr}[$currentExon]->[END_EXON]." pos start = ".$exons[$lastBorder]->{$chr}[$currentExon]->[START_EXON]."\n";
		if ($exons[$lastBorder]->{$chr}[$currentExon]->[$lastBorder] 
		    < $exons[$lastBorder]->{$chr}[$lastExon]->[$lastBorder]
		    && $exons[$lastBorder]->{$chr}[$currentExon]->[END_EXON] 
		    < $exons[$lastBorder]->{$chr}[$lastExon]->[START_EXON]
		    && $exons[$firstBorder]->{$chr}[$startingExon]->[END_EXON] 
		    < $exons[$lastBorder]->{$chr}[$currentExon]->[START_EXON]) {
                    
                    my @intersect = findSameGene(\%{$infoExons->{$chr.'@'.$exons[$lastBorder]->{$chr}[$currentExon]->[START_EXON]}},
                                                 \@current_gene);
                    if (@intersect) {
                        # Not all the transcripts have this exon, so we don't
                        # necessarily have to take it into account
                        # But we must update our set of allowable transcripts
                        if (@intersect < @current_gene) {
                            my %delh;
                            @delh{@intersect} = ();
                            @current_gene = grep ! exists $delh{$_}, @current_gene;
                        } else {
                            # if it is strictly comprised between the two exons of our junction
                            $nbNiceExonsMet++;
                            $lastExon = $currentExon;
                        }
                    }
		}
	    } while ($currentExon > 0 
		     && $exons[$lastBorder]->{$chr}[$currentExon]->[START_EXON] 
		     > $exons[$firstBorder]->{$chr}[$startingExon]->[END_EXON]);
# 	    print "startingExon = $startingExon -> pos = ". $exons[$firstBorder]->{$chr}[$startingExon]->[END_EXON]."\n";
	} else {
	    $isChimera = 1;
	}

	#mark the junction
	#$infoExons->{$chr.'@'.$exonStart}{'STATUT'}= 1;
	#$infoExons->{$chr_end.'@'.$exonEnd}{'STATUT'}= 1;

#	print "exonStart = $exonStart, exonEnd = $exonEnd\n";

	# classifying according to the specificity	
	if (abs($exonStart-$loc) <= VARIANT_THRESHOLD 
	    && abs($exonEnd-$loc_end) <= VARIANT_THRESHOLD) {
	    if (abs($exonStart-$loc) <= $bound_max
		&& abs($exonEnd-$loc_end) <= $bound_max) {
		if (! $isChimera 
		    && @current_gene) {
		    $expressedGenes->{$current_gene[0]}++;
		    if ($nbNiceExonsMet == 0) {
			print "overlap $tagNumber $chr:$exonStart-$exonEnd\n";
			if(noMultiplicity($chr,$exonStart,$chr_end,$exonEnd)){
			    $overlap++;
			}
		    } elsif ($nbNiceExonsMet >= 1) {
			print "splice $nbNiceExonsMet $tagNumber  $chr:$exonStart-$exonEnd\n";
			if(noMultiplicity($chr,$exonStart,$chr_end,$exonEnd)){
			    $splice++;
			}
		    } else  {
			# Just impossible... legacy stuff
			print "newExon [$chr:$loc;$loc_end] $tagNumber\n";
			if(noMultiplicity($chr,$exonStart,$chr_end,$exonEnd)){
			    $newExon++;
			}
		    } 
		    recordJunction($chr,$exonStart,$chr_end,$exonEnd);
		} else {
                    if ($chr eq $chr_end && $loc > $loc_end && $recursive_call == 0 && ! $stranded) {
		#	print "recursive call\n";
                        ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = 
			    crossWithRefSeq($expressedGenes, $knownExons, $startExons, $infoExons, $tagNumber, $loc_end, $loc, $chr_end, $chr, $strand_end, $strand, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,
					    $bound_max,$recursive_call+1);
                    } else {
                        print "chimera  [$chr:$loc;$chr_end:$loc_end] $tagNumber\n";
                        if(noMultiplicity($chr,$exonStart,$chr_end,$exonEnd)){
                            $intra++;
			    recordJunction($chr,$exonStart,$chr_end,$exonEnd);
                        }
                    }
		}
	    }elsif (abs($exonStart-$loc) <= $bound_max 
		    || abs($exonEnd-$loc_end) <= $bound_max 
		    ) {
		print "variant [$chr:$loc;$loc_end] $chr:$exonStart-$exonEnd $tagNumber\n";
		if(noMultiplicity($chr,$loc,$chr_end,$loc_end)){
		    $variant++;
		}
		recordJunction($chr,$loc,$chr_end,$loc_end);
	    }else{
		print "false positive [$chr:$loc;$loc_end] $chr:$exonStart-$exonEnd $tagNumber\n";
		if(noMultiplicity($chr,$loc,$chr_end,$loc_end)){
		    $falsePositive++;
		}
		recordJunction($chr,$loc,$chr_end,$loc_end);
	    } 
	}else{#print "recursive_call = $recursive_call, stranded = $stranded\n";
	    if ($recursive_call <= 1 && ! $stranded) {
		 ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = 
		     crossWithRefSeq($expressedGenes, $knownExons, $startExons, $infoExons, $tagNumber, $loc_end, $loc, $chr_end, $chr,$strand_end, $strand, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,
				     $bound_max, $recursive_call+1);
	    } else  {
		if (abs($exonStart-$loc) <= VARIANT_THRESHOLD 
		    || abs($exonEnd-$loc_end) <= VARIANT_THRESHOLD) {
		    print "newTranscript [$chr:$loc;$chr_end:$loc_end] $tagNumber\n";
		} else {
		    print "newTotalTranscript [$chr:$loc;$chr_end:$loc_end] $tagNumber\n";
		}
		if(noMultiplicity($chr,$loc,$chr_end,$loc_end)){
		    $newTranscript++;
		}
		recordJunction($chr,$loc,$chr_end,$loc_end);
	    }
	}
    }
    return ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive);
}

sub statSplicing{
    my ($expressedGenes,$overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = @_;
    my %expressedGenes = %$expressedGenes;
    print scalar(keys %expressedGenes), " genes\n";
    foreach my $gene (sort {$expressedGenes{$b} <=> $expressedGenes{$a}} keys %expressedGenes) {
	print $gene,": ",$expressedGenes{$gene},"\n";
    }
    
    my $total = $overlap+$splice+$newExon+$intra+$variant+$newTranscript+$falsePositive;
    printf "overlap: %d (%.2f%%)\n",$overlap,$overlap*100/$total;
    printf "splice: %d (%.2f%%)\n",$splice,$splice*100/$total;
    printf "newExon: %d (%.2f%%)\n",$newExon,$newExon*100/$total;
    printf "chimera: %d (%.2f%%)\n",$intra,$intra*100/$total;
    printf "variant: %d (%.2f%%)\n",$variant,$variant*100/$total;
    printf "newTranscript: %d (%.2f%%)\n",$newTranscript,$newTranscript*100/$total;
    printf "falsePositive: %d (%.2f%%)\n",$falsePositive,$falsePositive*100/$total;
}
#-------------------------------------------------------------------------
# end function
#-------------------------------------------------------------------------

if ($ARGV[0] eq "--stranded") {
    $stranded = 1;
    shift @ARGV;
}
my $db = shift @ARGV;
my $threshold = shift @ARGV;
my $bound_max = shift @ARGV;
my @crac = @ARGV;
my %knownExons;
my %startExons;
my $infoExons = {};
my @line;
my %geneNames;

# Read the DB file and fill the hashTables knownExons and infoExons
open(DB, $db) or die("Unable to open file $db: $!\n");
while (<DB>) {
    @line = split /,/;
    my ($gene, $startExon, $endExon, $strand, $chr) = 
	($line[1], $line[3], $line[4], $line[7], $line[17]);

    if ($startExon =~ /^[0-9]+$/ && $endExon =~ /^[0-9]+$/) {
# 	if ($stranded && $strand == -1) {
# 	    # The information is reversed on reverse strand
# 	    my $tmp = $endExon;
# 	    $endExon = $startExon;
# 	    $startExon = $tmp;
# 	}
	push @{$knownExons{$chr}}, [$endExon, $startExon];
	push @{$startExons{$chr}}, [$endExon, $startExon];
        # In fact not gene but transcript
	$geneNames{$chr.'@'.$startExon}{$gene}=1;
# 	push @{$infoExons->{$chr.'@'.$startExon}{'GENE'}}, $gene;
	$infoExons->{$chr.'@'.$startExon}{'STRAND'}= $strand;
	$infoExons->{$chr.'@'.$startExon}{'TYPE'}= 0; # START
	#$infoExons->{$chr.'@'.$startExon}{'STATUT'}= 0;
	$geneNames{$chr.'@'.$endExon}{$gene}=1;
# 	push @{$infoExons->{$chr.'@'.$endExon}{'GENE'}}, $gene;
	$infoExons->{$chr.'@'.$endExon}{'STRAND'}= $strand;
	$infoExons->{$chr.'@'.$endExon}{'TYPE'}= 1; # END
	#$infoExons->{$chr.'@'.$endExon}{'STATUT'}= 0;

    }
}
close(DB);

foreach my $gene (keys %geneNames) {
    @{$infoExons->{$gene}{'GENE'}} = keys %{$geneNames{$gene}};
}

foreach my $chr (keys %knownExons) {
    # Sorting limits of exons by ascending order
    @{$knownExons{$chr}} = sort {$a->[END_EXON] <=> $b->[END_EXON]} @{$knownExons{$chr}};
    @{$startExons{$chr}} = sort {$a->[START_EXON] <=> $b->[START_EXON]} @{$startExons{$chr}};
}



# my @test = @{$knownExons{'Y'}};
# foreach my $lim (@test) {
#     print $lim," ";
# }
# print "\n";

# print $infoExons->{'Y@26959330'}{'GENE'}."\t".$infoExons->{'Y@26959330'}{'TYPE'}."\t".$infoExons->{'Y@26959330'}{'STRAND'}."\n";

my %expressedGenes;
my $cluster = {};
my ($splice, $overlap, $intra, $newTranscript, $newExon, $variant, $falsePositive) = (0,0,0,0,0,0,0);
my ($tagNumber, $chr, $chr_end, $strand, $loc, $gap_length, $loc_end, $junction_name, $score, $nb_ins, $nb_del, $strand_end, $posShift, $gapLength);

foreach my $crac (@crac) {
# Read the CRAC file and search the good interval in %infoExons 
    open(CRAC, $crac) or die("Unable to open file $crac: $!\n");
    while (<CRAC>){
#	if (($tagNumber, $chr, $strand, $loc, $gap_length) = /^(\d+)\s+\S+\s+(\S+)?\|(\S+)?,(\d+).*?gap_length=(\d+)/){
	if (($tagNumber, $chr, $strand, $loc, $gap_length) = /^(\d+)\s+[a-z]+\s+(\S+)?\|(\S+)?,(\d+).*?gap_length=(\d+)/){	    
            # CRAC, overlap,splice case
	    if ($stranded && $strand == -1) {
		$loc_end = $loc;
		$loc = $loc_end + $gap_length;
	    } else {
		$loc_end = $loc+$gap_length;
	    }
# 	    print "CRAC start=$loc, end = $loc_end\n";
	    ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = crossWithRefSeq(\%expressedGenes,\%knownExons,\%startExons,$infoExons,$tagNumber,$loc,$loc_end,$chr,$chr,$strand, $strand, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,$bound_max);
	} 
#	elsif(($tagNumber, $chr, $strand, $loc, $nb_ins, $nb_del) = /^(\d+)\s+\S+\s+(\S+)?\|(\S+)?,(\d+).*nb_ins=(\d+).*nb_del=(\d+)/){
	elsif(($tagNumber, $chr, $strand, $loc, $nb_ins, $nb_del) = /^(\d+)\s+[a-z]+\s+(\S+)?\|(\S+)?,(\d+).*nb_ins=(\d+).*nb_del=(\d+)/){
            # CRAC, {genome, biotag}indel, snp,
	    if ($stranded && $strand == -1) {
		$loc_end = $loc;
		$loc = $loc_end + ($nb_del+$nb_ins);
	    } else {
		$loc_end = $loc+($nb_del+$nb_ins);
	    }
	    ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = crossWithRefSeq(\%expressedGenes,\%knownExons,\%startExons,$infoExons,$tagNumber,$loc,$loc_end,$chr,$chr,$strand, $strand, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,$bound_max);
	} 
#	elsif(($tagNumber, $chr, $strand, $loc, $chr_end, $strand_end, $loc_end) = /^(\d+)\s+\S+\s+(\S+)?\|(\S+)?,(\d+)?\s+(\S+)?\|(\S+)?,(\d+)\s+/){
	elsif(($tagNumber, $chr, $strand, $loc, $chr_end, $strand_end, $loc_end) = /^(\d+)\s+[a-z]+\s+(\S+)?\|(\S+)?,(\d+)?,\s+(\S+)?\|(\S+)?,(\d+)\s+/){
	    # chimera case
	    ($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = crossWithRefSeq(\%expressedGenes,\%knownExons,\%startExons,$infoExons,$tagNumber,$loc,$loc_end,$chr,$chr_end,$strand, $strand_end, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,$bound_max);
	} elsif (($chr, $loc, $loc_end, $junction_name, $score, $strand, $posShift, $gapLength) = $_ =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S*)\s+(\S+)\s+([+-])\s+[0-9]+\s+[0-9]+\s+[0-9,]+\s+[2-9][0-9]*\s+([0-9,]+)\s+([0-9,]+)$/){
	    # Tophat/MapSplice/GSNAP (bed format)
	    my $original_loc = $loc;
	    my $original_loc_end = $loc_end;
	    $chr =~ s/^chr//;
	    $tagNumber = $junction_name;
	    if ($strand eq '+'){
		$strand = 1;
	    }else{
		$strand = -1;
	    }
	    my @posShift = split /,/, $posShift;
	    my @gapLength = split /,/, $gapLength;
	    shift @gapLength;
	    while (@gapLength) {
		if ($stranded && $strand == -1) {
		    $loc_end = $loc;
		    $loc = $original_loc + $gapLength[0] - 1;
		    $loc_end += $posShift[0];
		} else {
		    $loc_end = $original_loc + $gapLength[0] - 1;
		    $loc += $posShift[0];
		}
# 		print "BED start=$loc, end = $loc_end\n";
		($overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive) = crossWithRefSeq(\%expressedGenes,\%knownExons,\%startExons,$infoExons,$tagNumber,$loc+1,$loc_end+1,$chr,$chr,$strand, $strand, $overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive,$bound_max);
		shift @gapLength;
		shift @posShift;
		if ($stranded && $strand == -1) {
		    $loc_end = $loc;
		} else {
		    $loc = $loc_end;
		}
	    }
	}
    }
    close(CRAC);
}
statSplicing(\%expressedGenes,$overlap,$splice,$newExon,$intra,$variant,$newTranscript,$falsePositive);    

# print "No seem junctions in bed file [Gene,strand]\t(startExon1,endExon1)->(startExon2,EndExon2) :\n";
# my $cmp = 0;
# foreach my $chr (keys %knownExons){
#     foreach my $exon (@{$knownExons{$chr}}){
# 	if ( $infoExons->{$chr.'@'.$exon}{'STATUT'} == 0){
# 	    if (($cmp % 4) == 0){
# 		print "\n[".$infoExons->{$chr.'@'.$exon}{'GENE'}."]\t($exon";
# 	    }elsif (($cmp % 2) == 0) {
# 		print "\t->\t ($exon";
# 	    }else{
# 		print ",$exon)";
# 	    }
# 	    $cmp++;
# 	}
#     }
# }
	
