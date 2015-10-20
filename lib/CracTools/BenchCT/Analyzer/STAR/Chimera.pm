use strict;
use warnings;
package CracTools::BenchCT::Analyzer::STAR::Chimera;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
 
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
    parsing_method => \&parseSTARChimeraLine,
    header_regex => '^#'
  );
  while (my $chimera = $chim_it->()) {
    my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
    
    if($true_chimera) {
      $self->getStats('chimera')->addTruePositive(id => $true_chimera);
    } else {
      #print STDERR Dumper($chimera);
      $self->getStats('chimera')->addFalsePositive(id => $chimera->{chr1}.";".$chimera->{pos1}.";".$chimera->{strand1}.";".$chimera->{chr2}.";".$chimera->{pos2}.";".$chimera->{strand2});
    }
  }
}

=head2 parseSTARChimeraLine

=cut

sub parseSTARChimeraLine {
  my $line = shift;
  my($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = split("\t",$line);
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => CracTools::Utils::convertStrand($strand1),
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => CracTools::Utils::convertStrand($strand2),
  };
}

1;
