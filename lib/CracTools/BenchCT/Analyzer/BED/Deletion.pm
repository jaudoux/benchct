use strict;
use warnings;
package CracTools::BenchCT::Analyzer::BED::Deletion;
# ABSTRACT: Analyze deletion bed files (as tophat produce)
 
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;
#use Data::Dumper;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'deletion') {
    return 1;
  }
  return 0;
}

sub _processLine {
  my $self = shift;
  my $bed_line = shift;

  my $true_deletion = $self->checker->isTrueDeletion($bed_line->{chr},$bed_line->{start},$bed_line->{end}-$bed_line->{start});
  if($true_deletion) {
    $self->getStats('deletion')->addTruePositive(id => $true_deletion);
  } else {
    #print STDERR Dumper($bed_line);
    $self->getStats('deletion')->addFalsePositive();
  }
}

1;
