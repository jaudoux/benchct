use strict;
use warnings;
package CracTools::BenchCT::Analyzer::BED::Insertion;
# ABSTRACT: Analyze deletion bed files (as tophat produce)
#
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;
#use Data::Dumper;

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

  if($self->checker->isTrueMutation('insertion',$bed_line->{chr},$bed_line->{start},length($bed_line->{name}))) {
    $self->getStats('insertion')->addTruePositive();
  } else {
    #print STDERR Dumper($bed_line);
    $self->getStats('insertion')->addFalsePositive();
  }
}

1;
