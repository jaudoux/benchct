#! /usr/bin/perl
#
use strict;
use warnings;

use Pod::Usage;
use CracTools::Utils;
use CracTools::BenchCT::Utils;

=head1 SYNOPSIS
  
  extractSamAlignmentsFromSpliceId.pl splices.bed alignments.bed mapping.sam id

=head1 AUTHOR

Jérôme Audoux

=cut

pod2usage() unless @ARGV >= 4;

my $splices = shift;
my $bed = shift;
my $sam = shift;
my $id = shift;
my $nb_lines = shift;

print STDERR "Id: $id\n";

my $i;
my $j;
my @reads_id;
#my $chim;
#my @reads_name;

my $bed_it = CracTools::Utils::bedFileIterator($splices);

# First we get the ids of the reads we are interested in
$i = 0;
while(my $bed_line = $bed_it->()) {
  if($i == $id) {
    #$splice = $bed_line;
    #print STDERR "Looking at chimera: ".join("\t",$chim->{chr1},$chim->{pos1},$chim->{strand1},$chim->{chr2},$chim->{pos2},$chim->{strand2})."\n";
    push @reads_id, split(":",$bed_line->{name});
    last;
  }
  $i++;
}

@reads_id = sort { $a <=> $b} @reads_id;

print STDERR "Found reads: ".join(",",@reads_id)."\n";
print STDERR "Now looking for SAM alignments...\n";

# Second we get the name of those reads
#my $bed_it = CracTools::Utils::getFileIterator(file => $bed,
#  parsing_method => \&CracTools::BenchCT::Utils::parseGSBedLineLite,
#);
my $bed_fh = CracTools::Utils::getReadingFileHandle($bed);
$i = 0;
$j = 0;
#while(my $bed_line = $bed_it->()) {
while(<$bed_fh>) {
    last if $j >= @reads_id;
    last if $j >= $nb_lines;
  if($i == $reads_id[$j]) {
      my $bed_line = CracTools::BenchCT::Utils::parseGSBedLineLite($_);
      print $_;
    print STDERR "Requesting: ".$bed_line->{chr}.":".$bed_line->{start}."\n";
    my $sam_it = CracTools::Utils::bamFileIterator($sam,$bed_line->{chr}.":".$bed_line->{start}."-".$bed_line->{end});
    #my $sam_it = CracTools::Utils::bamFileIterator($sam,$bed_line->{chr}.":".$bed_line->{start});
    while(my $sam_line = $sam_it->()) {
      if($sam_line =~ $bed_line->{name}) {
	  print $sam_line;
	  last;
      }
    }
    # $sam_it = CracTools::Utils::bamFileIterator($sam,$chim->{chr2}.":".$chim->{pos2}."-".$chim->{pos2});
    # while(my $sam_line = $sam_it->()) {
    #   if($sam_line =~ $bed_line->{name}) {
    #     print $sam_line;
    #     last;
    #   }
    # }
    #$reads_name[$j] = $bed_line->{name};

    $j++;
  }
  $i++;
}


