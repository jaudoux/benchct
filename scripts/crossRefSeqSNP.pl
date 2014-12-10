#!/usr/bin/perl -w
use strict;
use Getopt::Long;

use constant DEFAULT_THRESHOLD => 0;

sub usage {
    print STDERR "Usage: $0 [-t <nb>] [--false-positives file] [--true-positives file] <soft SNPs> <refSeq SNPs>

For both input files, the format must be the following one:
* first column: chromosome name (same names for both, of course!)
  the name of the chromosome may be prefixed by chr in one of the
  file and not in the other, I can deal with that;
* second column: position on the chromosome.

Files must be sorted on the chromosomes (lexicographic ascending order)
and on the positions (ascending order).

On STDOUT, the script displays the number of SNPs from the first
file that are seen in the second file.

-t is the tolerance for the locations on the chr (default ".DEFAULT_THRESHOLD.")\n";
    exit 1;
}

our ($help,$threshold, $falsePos, $truePos) = (0,DEFAULT_THRESHOLD, undef, undef);
GetOptions('help|?' => \$help,
	   'false-positives=s' =>\$falsePos,
	   'true-positives=s' => \$truePos,
	   't=i' => \$threshold) or exit 1;

if ($help) {
    usage;
}

if (@ARGV < 2) {
    usage;
}

my ($softFile, $refFile, $tpFile, $fpFile);
open($softFile, $ARGV[0]) or die ("Unable to open $ARGV[0]: $!\n");
open($refFile, $ARGV[1]) or die("Unable to open $ARGV[1]: $!\n");

open($tpFile, ">".$truePos) or die("Unable to write in $truePos: $!\n")
    if (defined ($truePos));
open($fpFile, ">".$falsePos) or die("Unable to write in $falsePos: $!\n")
    if (defined($falsePos));

my $softLine;
my $refLine;
my $countTP=0;
my $total = 0;
my $wasTP = 1; 			# Was the last SNP a TP?

my ($softChr, $softPos, $refChr, $refPos);

while (defined ($softLine = <$softFile>) && (! $wasTP || defined($refLine = <$refFile>))) {
    ($softChr, $softPos) = split /\s+/, $softLine;
    $softChr =~ s/^chr//;
#    print "** $softLine";

    if ($wasTP) {
	($refChr, $refPos) = split /\s+/, $refLine;
	$refChr =~ s/^chr//;
#	print "$refLine";
	$wasTP = 0;
    }

    while (defined $refLine && 
	   (($refChr eq $softChr && $refPos < $softPos - $threshold)
	   || ($refChr lt $softChr))){
	$refLine = <$refFile>;
	if (defined($refLine)) {
	    ($refChr, $refPos) = split /\s+/, $refLine;
	    $refChr =~ s/^chr//;
	}
#   	print "\t$refLine";
    }

    $total++;
#   print "$refChr == $softChr: ".($refChr eq $softChr)." abs($refPos - $softPos) = ".abs($refPos - $softPos)."<= $threshold\n";
    if ($refChr eq $softChr && abs($refPos - $softPos) <= $threshold) {
	$countTP++;
	print $tpFile "$softChr\t$softPos\n"
	    if (defined($truePos));
	$wasTP = 1;
    } elsif(defined $falsePos) {
#	print "FP: $softLine";
	print $fpFile "$softChr\t$softPos\n";
    }
}
#print "Done (almost)\n";
while (<$softFile>) {
    $total++;
    ($softChr, $softPos) = split /\s+/;
    $softChr =~ s/^chr//;
    if(defined $falsePos) {
	print $fpFile "$softChr\t$softPos\n";
    }
}

close($softFile);
close($refFile);
close($tpFile)
    if (defined($truePos));
close($fpFile)
    if(defined($falsePos));

printf "TP:\t%d (%.2f%%)\nTotal:\t%d\n", $countTP, $countTP*100./$total,$total;
