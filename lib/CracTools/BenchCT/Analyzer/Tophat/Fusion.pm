use strict;
use warnings;
package CracTools::BenchCT::Analyzer::Tophat::Fusion;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;
#use Data::Dumper;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'chimera') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $chimera_file = $args{file};
  my $chim_it = CracTools::Utils::getFileIterator(file => $chimera_file,
    parsing_method => \&parseTophatFusionLine,
  );
  while (my $chimera = $chim_it->()) {
    my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
    
    if($true_chimera) {
      $self->getStats('chimera')->addTruePositive($true_chimera);
    } else {
      #print STDERR Dumper($chimera);
      $self->getStats('chimera')->addFalsePositive();
    }
  }
}

=head2 parseTophatFusionLine

  chr20-chr17 49411707  59445685  ff  106 116 167 0 37  36  0.569598  @ 
  11 25 38 49 63 @ 
  CAGCGGGGCGCGCGAGCTCGCGCTCTTCCTGACCCCCGAGCCTGGGGCCG AGGTAGGGGACGGGGCTGTGGAGTTGGAGGAGAGGGTTCTCGCGGTTAGG @ 
  CCTGCTCCCTGAAGGTGTGGACTCAACGTCAGATGTCCCGTGTGTGCCAC AGGTACCTTTGACAGGAGCGTGACCCTGCTGGAGGTGTGCGGGAGCTGGC @ 
  106 106 106 106 106 106 106 106 106 106 106 106 106 106 98 90 84 79 79 78 74 68 65 64 63 63 59 59 56 55 52 51 20 15 12 10 8 0 0 0 0 0 0 0 0 0 0 0 0 0 @ 
  106 106 106 106 106 106 106 106 106 106 106 106 106 98 96 94 91 86 55 54 51 50 47 47 43 43 42 41 38 32 28 27 27 22 16 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 @ 
  -6:1 11:0 16:-3 18:1 14:6 14:6 15:7 23:0 5:21 31:5 18:19 36:-1 ... 


=cut

sub parseTophatFusionLine {
  my $line = shift;
  my ($chr1,$chr2,$pos1,$pos2,$strand1,$strand2,$spanning_junction,$spanning_PE) = $line =~ /(\S+)?-(\S+)\s+(\d+)\s+(\d+)\s+(\S)(\S)\s+(\d+)/;
  my %convert_hash = ('f' => 1, 'r' => -1);
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => $convert_hash{$strand1},
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => $convert_hash{$strand2},
    spanning_junction => $spanning_junction,
    spanning_PE => $spanning_PE,
  };
}

1;
