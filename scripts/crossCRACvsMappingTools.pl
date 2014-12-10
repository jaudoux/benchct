#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;

if (@ARGV < 5) {
    print STDERR "Usage: $0 [<file_bowtie> || <file_bwa> || <file_soap>] <file_single_crac> <file_duplicate_crac> <file_multiple_crac> <file_none_crac> [other_file_crac]+\n";
    exit(1);
}

my $filename = basename($ARGV[0]);
my ($libId) = $filename =~ /\-(.*)\.[a-zA-Z]+$/;

# -------------------  get all otherTools mapped reads -------------------
open(OTHER,$ARGV[0]) or die(impossible to open $!);
my $totalSingleOther++;
my $totalMultipleOther++;
my %hash = ();
my $oldid = 0;
my ($id,$strand,$tag,$type,$tagtmp,$flag);
my $otherTools;
while(<OTHER>){
    chomp;
    if (($id,$strand, $tag) = /^(.*)\s+([\-\+])\s+\S+\s+\d+\s+([CATG]+)/){
	# bowtie
	if (!defined $otherTools){
	    $otherTools = "Bowtie";
	}
        #print "$idtag,$strand,$tagbowtie\n";
	if ($strand eq '-'){
	    $tagtmp = scalar reverse $tag;
	    $tag = $tagtmp; 
	    $tag =~ tr/ACTG/TGAC/;
	}
	#print "$tagbowtie\n";
	if ($oldid eq $id){
	    $totalSingleOther--;
	    $totalMultipleOther++;
	    $hash{$tag}{'multiple'} = $_;
	    $hash{$tag}{'single'} = undef;
	}else{
	    $totalSingleOther++;
	    $hash{$tag}{'single'} = $_;
	}
	$oldid = $id;
    }elsif (($flag,$tag,$type) = /^\S+\s+(\S+).*?([CATG]+)\s+.*\s+XT:A:(\S)/){
	# bwa
	if (!defined $otherTools){
	    $otherTools = "BWA";
	}
	if ($flag & 16){
	    $tagtmp = scalar reverse $tag;
	    $tag = $tagtmp; 
	    $tag =~ tr/ACTG/TGAC/;
	}
	if ($type eq 'R'){
	    $hash{$tag}{'multiple'} = $_;
	    $totalMultipleOther++;
	}elsif ($type eq 'U'){
	    $hash{$tag}{'single'} = $_;
	    $totalSingleOther++;
	}
    }elsif (($tag,$type,$strand) = /^.*?\s+([CATG]+)\s+\S+\s+(\d)\s+.*?\s+([\-\+])\s+/){
        #soap2
	if (!defined $otherTools){
	    $otherTools = "SOAP2";
	}	
	if ($strand eq '-'){
	    $tagtmp = scalar reverse $tag;
	    $tag = $tagtmp; 
	    $tag =~ tr/ACTG/TGAC/;
	}
	if ($type == 1){
	    $hash{$tag}{'single'} = $_;
	    $totalSingleOther++;
	}else{
	    $hash{$tag}{'multiple'} = $_;
	    $totalMultipleOther++;
	}
    }
}
close(OTHER);

if (! defined $otherTools) {
    print STDERR "Sorry the files you gave me have not been recognised\n";
    exit 3;
}

# -------------------- cross with CRAC mapping step -------------------------
print "first step : what is the distribution of $otherTools single mapped reads in CRAC ?\n\n";
my (@totalsuff, @commons, @onlycracs, @suffFile, @hl_suffFile);
my ($percentCommonSingleSuffCRAC, $percentCommonSingleSuffOtherTools, $percentCommonMultipleSuffCRAC, $percentCommonMultipleSuffOtherTools);  

for (my $i =1 ; $i <= $#ARGV ; $i++){
    ($suffFile[$i]) = $ARGV[$i] =~ /^.*\.(\S+)$/;
    open($hl_suffFile[$i],">CRAC.inter$otherTools-$libId\_$suffFile[$i]") or die("impossible to open CRAC.inter$otherTools\_$suffFile[$i]");
    print {$hl_suffFile[$i]} $suffFile[$i]."type of CRAC classification\tinfo(CRAC)\tsingle or multiple dans $otherTools(S/M)\tinfo($otherTools)\n";
}
open(OUT,">CRAC.no$otherTools-$libId\_mapped") or die("impossible to open CRAC.no$otherTools\_mapped");
printf("Type\tTotal\tSingle_$otherTools\tMultiple_$otherTools\t%%single_$otherTools\t%%multiple_$otherTools\n");
for (my $i =1 ; $i <= $#ARGV ; $i++){
    open(CRAC,$ARGV[$i]) or die("impossible to open $ARGV[$i] : $!");
    $totalsuff[$i] = 0;
    $totalsuff[$i] = 0;
    $commons[$i]{'single'} = 0;
    $commons[$i]{'multiple'} = 0;
    $onlycracs[$i] = 0;
    ($percentCommonSingleSuffCRAC, $percentCommonSingleSuffOtherTools, $percentCommonMultipleSuffCRAC, $percentCommonMultipleSuffOtherTools)
	= (0,0,0,0);  
    while(<CRAC>){
	if (my ($tagcrac) = $_ =~ /\s+([CATG]+)\s+\d+/){
	    $totalsuff[$i]++;
	    if (defined $hash{$tagcrac}{'single'}){
		$commons[$i]{'single'}++;
		print {$hl_suffFile[$i]} $hash{$tagcrac}{'single'}."\tS\t$_";
	    }elsif (defined $hash{$tagcrac}{'multiple'}){
		$commons[$i]{'multiple'}++;
		print {$hl_suffFile[$i]} $hash{$tagcrac}{'multiple'}."\tM\t$_";
	    }else{
		if ($i <= 3){
		    print OUT "$tagcrac\n";
		#}elsif ($suffFile[$i] eq 'normal'){
		 #   print "$tagcrac\n";
		}
		$onlycracs[$i]++;
	    }
	}
    }
    close(CRAC);
    if ($totalsuff[$i] > 0){
	$percentCommonSingleSuffOtherTools = $commons[$i]{'single'}*100/$totalSingleOther;
	$percentCommonSingleSuffCRAC = $commons[$i]{'single'}*100/$totalsuff[$i];
	$percentCommonMultipleSuffOtherTools = $commons[$i]{'multiple'}*100/$totalMultipleOther;
	$percentCommonMultipleSuffCRAC = $commons[$i]{'multiple'}*100/$totalsuff[$i];
    }
    printf("%s\t%d\t%d\t%d\t%.2f\t%.2f\n",$suffFile[$i],$totalsuff[$i],$commons[$i]{'single'},$commons[$i]{'multiple'},
	   $percentCommonSingleSuffOtherTools, $percentCommonMultipleSuffOtherTools);
#     printf("commons single $otherTools and %s CRAC : %d (%.2f %% of single $otherTools and %.2f %% of %s CRAC)\n",$suffFile[$i],$commons[$i]{'single'},$percentCommonSingleSuffOtherTools,$percentCommonSingleSuffCRAC,$suffFile[$i]);
#     printf("commons multiple $otherTools and %s CRAC : %d (%.2f %% of multiple $otherTools and %.2f %% of %s CRAC)\n",$suffFile[$i],$commons[$i]{'multiple'},$percentCommonMultipleSuffOtherTools,$percentCommonMultipleSuffCRAC,$suffFile[$i]);
}
close(OUT);
for (my $i =1 ; $i <= $#ARGV ; $i++){
    close($hl_suffFile[$i]);
}

# -------------------- outside the intersection-------------------------
print "--------------------------------------------\n\n";
print "second step : mapped reads outside the intersection ?\n\n";
my $percentOnlyCRAC;
for (my $i = 1 ; $i <= 3 ; $i++){
    $percentOnlyCRAC = 0;    
    if ($totalsuff[$i] > 0){
	$percentOnlyCRAC = $onlycracs[$i]*100/$totalsuff[$i];
    }
    printf("only mapped CRAC : $onlycracs[$i] (%.2f %% of %s mapped CRAC)\n",$percentOnlyCRAC,$suffFile[$i]);
}

print "--------------------------------------------\n\n";
print "third step : what is the classification of the supplemental CRAC mapped reads ?\n\n"; 

my %hashcrac = ();
open(CRAC,"CRAC.no$otherTools-$libId\_mapped") or die("impossible to open CRAC.no$otherTools\_mapped : $!");
my $countMapped = 0;
while (<CRAC>){
    chomp;
    if ($_ =~ /^[CATG]+/){
	$hashcrac{$_} = 1;
	$countMapped++;
    }
}
close(CRAC);

my ($countSuff, $countTotal);
my ($percentMapped,$percentSuff);
for (my $i=1; $i<=$#ARGV ; $i++){
    open(CRAC,$ARGV[$i]) or die(impossible to open $!);
    $countTotal = 0;
    $countSuff = 0;
    $percentSuff = 0;
    $percentMapped = 0;
   while(<CRAC>){
	if (my ($tagcrac) = $_ =~ /\s+([CATG]+)\s+\d+/){
	    if (defined $hashcrac{$tagcrac}){
		$countSuff++;
	    }
	    $countTotal++;
	}
    }
    close(CRAC);
    if ($countTotal > 0 && $countMapped > 0){
	$percentMapped = $countSuff*100/$countMapped;
	$percentSuff = $countSuff*100/$countTotal;
    }
    printf("number of CRAC %s : %d (%.2f %% of supplementary mapped CRAC and %.2f %% of %s CRAC) \n", $suffFile[$i], $countSuff, $percentMapped, $percentSuff, $suffFile[$i]);
}
