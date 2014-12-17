use strict;
use warnings;
package CracTools::BenchCT::Analyzer::BED::Junction;
# ABSTRACT: Analyze junction bed fils (as tophat produce)
#
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;
#use Data::Dumper;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'splice') {
    return 1;
  }
  return 0;
}

sub _processLine {
  my $self = shift;
  my $bed_line = shift;
  # Loop over blocks
  for(my $i=1; $i < @{$bed_line->{blocks}}; $i++) {

    # Extract splice coordinates
    my $chr = $bed_line->{chr};
    my $strand = CracTools::Utils::convertStrand($bed_line->{strand});
    my $start = $bed_line->{blocks}[$i-1]->{ref_end};
    my $end = $bed_line->{blocks}[$i]->{ref_start};
    my $length = $end - $start;

    if($self->checker->isTrueSplice($chr,$start,$length,$strand)) {
      $self->getStats('splice')->addTruePositive();
    } else {
      #print STDERR Dumper($bed_line);
      $self->getStats('splice')->addFalsePositive();
    }
  }
}

1;
