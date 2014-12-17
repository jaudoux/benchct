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
    #interval_query => CracTools::Interval::Query->new(),
    nb_events => 0,
    threshold => $args{threshold},
    verbose   => defined $args{verbose}? $args{verbose} : 0,
    #event_reads => [],
    #genome_mask => CracTools::GenomeMask->new(),
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

Add a new event to the collection

=cut

sub addEvent {
  my $self = shift;
  $self->{nb_events}++;
  #my $info_line = shift;
  #push(@{$self->{event_reads}},$info_line->{read_ids});
}

sub getLastEventId {
  my $self = shift;
  return $self->nbEvents - 1;
}

#sub getEventReads {
#  my $self = shift;
#  my $i = shift;
#  return $self->{event_reads}->[$i];
#}

sub nbEvents {
  my $self = shift;
  return $self->{nb_events};
}

1;
