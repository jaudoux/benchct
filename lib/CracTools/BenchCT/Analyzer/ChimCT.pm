use strict;
use warnings;
package CracTools::BenchCT::Analyzer::ChimCT;
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
  my $min_crac_score = $args{options}->{min_crac_score};
  my $chim_it = CracTools::Utils::chimCTFileIterator($chimera_file);
  while (my $chimera = $chim_it->()) {
    my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
    
    if(!defined $chimera->{comments}->{CRAC_score} || !defined $min_crac_score || $chimera->{comments}->{CRAC_score} >= $min_crac_score) {
      if($true_chimera) {
        $self->getStats('chimera')->addTruePositive($true_chimera);
      } else {
        #print STDERR Dumper($chimera);
        $self->getStats('chimera')->addFalsePositive($chimera->{chim_key});
      }
    }
  }
}

1;
