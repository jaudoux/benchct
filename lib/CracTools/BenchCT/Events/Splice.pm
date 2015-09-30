use strict;
use warnings;
package CracTools::BenchCT::Events::Splice;
# ABSTRACT: A collection of events of type 'splice' on the genome

use parent 'CracTools::BenchCT::Events::Interval';

use CracTools::Interval::Query;

sub addSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_; 
  return $self->addInterval($chr,$start,$start+$length-1,$strand);
}

sub isTrueSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_; 
  return $self->isTrueInterval($chr,$start,$start+$length-1,$strand,$self->threshold);
}

1;
