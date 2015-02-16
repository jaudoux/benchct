#! /usr/bin/perl
#
use strict;
use warnings;

use Pod::Usage;
use CracTools::Utils;
use CracTools::BenchCT::Utils;

=head1 SYNOPSIS
  
  extractSamAlignmentsFromChimId.pl chimeras.tsv alignments.bed mapping.sam id

=head1 AUTHOR

Jérôme Audoux

=cut

pod2usage() unless @ARGV >= 4;

my $chimeras = shift;
my $bed = shift;
my $sam = shift;
my $id = shift;

print STDERR "Id: $id\n";

my $i;
my $j;
my @reads_id;
my $chim;
#my @reads_name;

my $chim_it = CracTools::Utils::getFileIterator(file => $chimeras,
  parsing_method => \&CracTools::BenchCT::Utils::parseChimeraLine,
);

# First we get the ids of the reads we are interested in
$i = 0;
while(my $chim_line = $chim_it->()) {
  if($i == $id) {
    $chim = $chim_line;
    print STDERR "Looking at chimera: ".join("\t",$chim->{chr1},$chim->{pos1},$chim->{strand1},$chim->{chr2},$chim->{pos2},$chim->{strand2})."\n";
    push @reads_id, @{$chim_line->{read_ids}};
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
  if($i == $reads_id[$j]) {
    my $bed_line = CracTools::BenchCT::Utils::parseGSBedLineLite($_);
    print STDERR "Requesting: ".$bed_line->{chr}.":".$bed_line->{start}."\n";
    my $sam_it = CracTools::Utils::bamFileIterator($sam,$chim->{chr1}.":".$chim->{pos1}."-".$chim->{pos1});
    #my $sam_it = CracTools::Utils::bamFileIterator($sam,$bed_line->{chr}.":".$bed_line->{start});
    while(my $sam_line = $sam_it->()) {
      if($sam_line =~ $bed_line->{name}) {
        print $sam_line;
        last;
      }
    }
    $sam_it = CracTools::Utils::bamFileIterator($sam,$chim->{chr2}.":".$chim->{pos2}."-".$chim->{pos2});
    while(my $sam_line = $sam_it->()) {
      if($sam_line =~ $bed_line->{name}) {
        print $sam_line;
        last;
      }
    }
    #$reads_name[$j] = $bed_line->{name};

    $j++;
  }
  $i++;
}


