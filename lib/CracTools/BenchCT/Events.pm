package CracTools::BenchCT::Events;
# ABSTRACT: A generic collection of events on the genome

use Moo;
use Carp;

has threshold => (is  => 'rw', default => 0);
has verbose   => (is  => 'rw', default => 0);

has events    => (
  is      => 'ro',
  default => sub { return [] },
);

1;

=head2 new

=cut

=head2 addEvent

  my $event_id = $self->addEvent($event);

Add a new event to the collection and return its id

=cut

sub addEvent($) {
  my $self = shift;
  my $event = shift;
  # Add the new event
  push(@{$self->events},$event) if defined $event;
  # Return the last event ids
  return $self->getLastEventId;
}

=head2 printHeader
  
  $events->printHeader($fh)

Print Header line(s) in the output stream

=cut

sub printHeader($) {
  my $self = shift;
  my $fh = shift;
  print $fh "event_id\n";
}

=head2 printEvent

  $events->printEvent($id,$fh)

Print the event on the output stream

=cut

sub printEvent($$) {
  my $self = shift;
  my $fh = shift;
  my $event_id = shift;
  print $fh $event_id,"\n";
}

=head2 getEvent

  my $event = $self->getEvent($event_id);

Return the event associated to the id in parameter

=cut

sub getEvent($) {
  my $self = shift;
  my $id = shift;
  croak "Supplied id ($id) is not valid" unless $id < $self->nbEvents;
  return $self->events->[$id];
}

=head2 getLastEventId

Return the id of the last event added to the structure

=cut

sub getLastEventId {
  my $self = shift;
  return $self->nbEvents - 1;
}

=head2 nbEvents

Return the number of events added to the structure

=cut

sub nbEvents {
  my $self = shift;
  return scalar @{$self->events};
}

1;
