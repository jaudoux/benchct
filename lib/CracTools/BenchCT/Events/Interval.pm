use strict;
use warnings;
package CracTools::BenchCT::Events::Interval;
# ABSTRACT: A generic collection of events that are intervals over the genome

use parent 'CracTools::BenchCT::Events';

use CracTools::Interval::Query;
use List::Util qw(min max);

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

  $self->{interval_query} = CracTools::Interval::Query->new();
  
  return $self;
}

sub intervalQuery {
  my $self = shift;
  return $self->{interval_query};
}

sub addInterval {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_; 

  # First we look if this exon already exists using a 0 threshold
  my $id = $self->_foundInterval($chr,$start,$end,$strand,0);
  if($id) {
    return $id - 1;
  } else {
    $id = $self->addEvent("$start-$end"); # Increment the number of exons
    $self->intervalQuery->addInterval($chr,$start,$end,$strand,$id);
    return $id;
  }
}

sub isTrueInterval {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_; 
  return $self->_foundInterval($chr,$start,$end,$strand,$self->threshold);
}

sub _foundInterval {
  my $self = shift;
  my ($chr,$start,$end,$strand,$threshold) = @_; 

  # Check if we have already seen this exon
  # The 1 at the we means that we use a "windowed" query
  # TODO there is a bug in Set::IntervalTree windowed query, we are no longer using them
  my @matching_intervals = @{$self->intervalQuery->fetchByRegion($chr,$start-$threshold,$end+$threshold,$strand,0)};
  my $best_overlap = 0;
  my $best_interval;
  foreach my $interval (@matching_intervals) {
    my ($interval_start,$interval_end) = split('-',$self->getEvent($interval));
    my $overlap_length = min($interval_end,$end) - max($interval_start,$start) + 1;
    my $total_length   = max($interval_end,$end) - min($interval_start,$start) + 1;
    my $overlap_cover  = $overlap_length/$total_length;
    if(($total_length - $overlap_length) <= $threshold && $overlap_cover > $best_overlap) {
      $best_interval    = $interval;
      $best_overlap = $overlap_cover;
      last if $best_overlap == 1;
    }
  }
  
  if(defined $best_interval) {
    return $best_interval + 1;
  } else {
    return 0;
  }
}

1;
