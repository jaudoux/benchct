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
  my $chim_it = CracTools::Utils::chimCTFileIterator($chimera_file);
  while (my $chim_line = $chim_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($chim_line);
  }
}

sub _processLine {
  my $self = shift;
  my $chimera = shift;
  my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
  
  if(!defined $chimera->{comments}->{CRAC_score} || $chimera->{comments}->{CRAC_score} > 70) {
    if($true_chimera) {
      $self->getStats('chimera')->addTruePositive($true_chimera);
    } else {
      #print STDERR Dumper($chimera);
      $self->getStats('chimera')->addFalsePositive();
    }
  }
}

1;
