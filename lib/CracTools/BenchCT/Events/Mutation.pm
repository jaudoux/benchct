use strict;
use warnings;
package CracTools::BenchCT::Events::Mutation;
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

=head2 addEvent

Add a new event to the collection

=cut

sub addMutation {
  my $self  = shift;
  my $info_line = shift;
  # Add the length of the mutation as the event value
  my $id = $self->addEvent($info_line->{length});
  $self->intervalQuery->addInterval($info_line->{chr},
    $info_line->{old_pos},
    $info_line->{old_pos},
    undef, # A mutation does not have a strand....
    $id,
  );
}

sub isTrueMutation {
  my $self = shift;
  my ($chr,$pos, $length) = @_;
  my @mutations = @{$self->intervalQuery->fetchByRegion($chr,$pos - $self->threshold,$pos + $self->threshold)};
  # We return true if we have found a matching mutation that have the same length
  foreach my $mutation (@mutations) {
    my $mutation_length = $self->getEvent($mutation);
    if(abs($mutation_length - $length) <= $self->threshold) {
      return $mutation + 1;
    }
  }
  return 0;
}

1;
