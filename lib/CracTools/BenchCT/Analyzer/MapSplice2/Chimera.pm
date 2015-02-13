use strict;
use warnings;
package CracTools::BenchCT::Analyzer::MapSplice2::Chimera;
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
    parsing_method => \&parseMapSplice2ChimeraLine,
    header_regex => '^#'
  );
  while (my $chimera = $chim_it->()) {
    my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
    
    if($true_chimera) {
      $self->getStats('chimera')->addTruePositive(id => $true_chimera);
    } else {
      #print STDERR Dumper($chimera);
      $self->getStats('chimera')->addFalsePositive();
    }
  }
}

=head2 parseMapSplice2ChimeraLine

  chr5 60039689 60039689 JUNC_22189 638 + 60039689 60039689 255,0,0 2 30,44, 0,45, 0.126481 0 1 1 0 4 0.805643 51 587 1 0 1 0 1 637 0 C-637,T-1,

=cut

sub parseMapSplice2ChimeraLine {
  my $line = shift;
  my($chr1,$chr2,$pos1,$pos2,$cover,$strand1,$strand2) = $line =~ /(\S+)?-(\S+)\s+(\d+)\s+(\d+)\s+\S+(\d+)\s+(\S)(\S)/;
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => CracTools::Utils::convertStrand($strand1),
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => CracTools::Utils::convertStrand($strand2),
    cover => $cover,
  };
}

1;
