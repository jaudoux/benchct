use strict;
use warnings;
package CracTools::BenchCT::Analyzer::Tophat::Insertion;
# ABSTRACT: Analyze deletion bed files (as tophat produce)
 
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'insertion') {
    return 1;
  }
  return 0;
}

sub _processLine {
  my $self = shift;
  my $bed_line = shift;

  my $true_insertion = $self->checker->isTrueInsertion($bed_line->{chr},$bed_line->{start},length($bed_line->{name}));
  if($true_insertion) {
    $self->getStats('insertion')->addTruePositive(id => $true_insertion);
  } else {
    $self->getStats('insertion')->addFalsePositive();
  }
}

1;
