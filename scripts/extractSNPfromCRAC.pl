#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use constant DEFAULT_SCORE => 0;
sub usage {
    print STDERR "Usage: $0 [--score <float>] [--duplicate] [--report-once] [--forbidden <crac IDs>]+ [--mandatory <crac IDs>]+ <file of crac> 

The input Crac file, must be the following one:
* the SNP file of CRAC 

--forbidden is a list of ID that won't be output as SNP
--mandatory is a list of ID where all the output SNPs will be

The IDs are the read numbers in the input file.

On STDOUT, the script displays \"chr pos\" for each reads considered.
Beware! Positions start at 0.

By default only the unique causes are considered but with --duplicate option the duplicate causes are kept.

--score (necessary <= 0) is the minimal score given by Crac to consider a SNP
(default ".DEFAULT_SCORE.")

--report-once: when the option is specified, a given SNP is reported just
    once. A third column is output, containing the coverage of the SNP.\n";
    exit 1;

}

our ($help,$score, $duplicate, $reportOnce) = (0,DEFAULT_SCORE, 0, 0);
my %report;
my @forbidden_ids;
my @mandatory_ids;
GetOptions('help|?' => \$help,
	   'duplicate' =>\$duplicate,
	   'score=i' => \$score,
           'report-once' => \$reportOnce,
	   'forbidden=s' => \@forbidden_ids,
	   'mandatory=s' => \@mandatory_ids) or exit 1;

if ($help || @ARGV < 1) {
    usage;
}

my $crac_snp=$ARGV[0];


my %forbidden_ids;
my %mandatory_ids;
    
my $handle;
foreach my $file (@forbidden_ids) {
    open($handle, $file) or die("Unable to open $file: $!\n");
    my $id;
    while(<$handle>) {
	if (($id) = /^([0-9]+)\s/) {
	    $forbidden_ids{$id}=1;
	}
    }
    close($handle);
}
foreach my $file (@mandatory_ids) {
    open($handle, $file) or die("Unable to open file: $!\n");
    my $id;
    while(<$handle>) {
	if (($id) = /^([0-9]+)\s/) {
	    $mandatory_ids{$id}=1;
	}
    }
    close($handle);
}


my @field;

open(CRAC_IN,$crac_snp) or die("Unable to open $crac_snp: $!\n");

while(<CRAC_IN>){
    # format without duplication: 10 single 19|-1,110856 pos_SNV=28 G->A score=-17.5553 15|1,102463096 pos_location=41 GACTGAAGACAGAGACAGATCAATGAGTAAGAGGTTGGCTAGCAGGAAGTACATGGGAGAGTGAAGGTGGGAGTC 
    # format with duplication: 10 duplicate 19|-1,110856 pos_SNV=28 G->A
    # score=-17.5553 15|1,102463096 pos_location=41
    # GACTGAAGACAGAGACAGATCAATGAGTAAGAGGTTGGCTAGCAGGAAGTACATGGGAGAGTGAAGGTGGGAGTC 
    if (/^#/) {
      next;
    }
    @field = split(/\s+/, $_);
    if ((! %mandatory_ids || defined $mandatory_ids{$field[0]})
	&& ! defined($forbidden_ids{$field[0]})) {
	my ($s,$single,$chr,$loc,$pos_single) = ($field[5], $field[1], $field[2], $field[2]);
	if ($s =~ /score=/) {
	    $s =~ s/score=//;
	} 
	$chr =~ s/\|.*//;
	$loc =~ s/.*,//;
	if (! defined $s || $s <= $score){
	    if ($single eq 'single' || $duplicate){
		if ($reportOnce) {
                  $report{$chr}->{$loc}++;
                } else {
		  print "$chr $loc\n";
                }
	    }
	}
    }
}

if ($reportOnce) {
    foreach my $chr (sort keys %report) {
        foreach my $pos (sort {$a <=> $b} keys %{$report{$chr}}) {
            print "$chr\t$pos\t".$report{$chr}->{$pos}."\n";
        }
    }
}
close(CRAC_IN);
