use strict;
use warnings;
package CracTools::BenchCT::Events;
# ABSTRACT: A generic collection of events on the genome

=head2 new

=cut

sub new {
  my $class = shift;

  # Get args
  my %args = @_;

  my $self = bless {
    nb_events => 0,
    threshold => $args{threshold},
    verbose   => defined $args{verbose}? $args{verbose} : 0,
    events    => [],
  }, $class;

  return $self;
}

sub threshold {
  my $self = shift;
  return $self->{threshold};
}

sub verbose {
  my $self = shift;
  return $self->{verbose};
}

=head2 addEvent

  my $event_id = $self->addEvent($event);

Add a new event to the collection and return its id

=cut

sub addEvent {
  my $self = shift;
  my $event = shift;
  # Add the new event
  push(@{$self->{events}},$event) if defined $event;
  # Increment the counter
  $self->{nb_events}++;
  # Return the last event ids
  return $self->getLastEventId;
}

=head2 printHeader
  
  $events->printHeader($fh)

Print Header line(s) in the output stream

=cut

sub printHeader {
  my $self = shift;
  my $fh = shift;
  print $fh "event_id\n";
}

=head2 printEvent

  $events->printEvent($id,$fh)

Print the event on the output stream

=cut

sub printEvent {
  my $self = shift;
  my $event_id = shift;
  my $fh = shift;
  print $fh $event_id,"\n";
}

=head2 getEvent

  my $event = $self->getEvent($event_id);

Return the event associated to the id in parameter

=cut

sub getEvent {
  my $self = shift;
  my $id = shift;
  return $self->{events}[$id];
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
  return $self->{nb_events};
}

1;
