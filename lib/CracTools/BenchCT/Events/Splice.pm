use strict;
use warnings;
package CracTools::BenchCT::Events::Splice;
# ABSTRACT: A collection of events of type 'mutation' (snp,indel) on the genome

use parent 'CracTools::BenchCT::Events';

use CracTools::Interval::Query;

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

sub addSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_; 
    
  # If we have not found the splice we add it to the IntervalQuert
  #if(!$self->_foundSplice($chr,$start,$length,$strand,1)) {
  $self->intervalQuery->addInterval($chr,$start,$start+$length,$strand,$length);
  $self->addEvent(); # Increment the number of splices
  #} else {
  #}
}

sub isTrueSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_; 
  return $self->_foundSplice($chr,$start,$length,$strand,$self->threshold);
  #return $self->_foundSplice($chr,$start,$length,1,$self->threshold) || $self->_foundSplice($chr,$start,$length,-1,$self->threshold) || 
  #  $self->_foundSplice($chr,$start-$length,$length,1,$self->threshold) || $self->_foundSplice($chr,$start-$length,$length,-1,$self->threshold) ||
  #  $self->_foundSplice($chr,$start+$length,$length,1,$self->threshold) || $self->_foundSplice($chr,$start+$length,$length,-1,$self->threshold);
}

sub _foundSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand,$threshold) = @_; 

  # Check if we have already seen this splice
  # The 1 at the we means that we use a "windowed" query
  # TODO there is a bug in Set::IntervalTree windowed query, we are no longer using them
  my @matching_splices = @{$self->intervalQuery->fetchByRegion($chr,$start-$threshold,$start+$length+$threshold,$strand,0)};
  my $found_splice = 0;

  foreach my $splice (@matching_splices) {
    if(abs($splice - $length) <= $threshold) {
      $found_splice = 1;
      last;
    }
  }

  return $found_splice;
}

1;
