use strict;
use warnings;
package CracTools::BenchCT::Analyzer;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use CracTools::BenchCT::Stats;
use Carp;

our %check_events = (
  "insertion"  => 0,
  "mapping"    => 0,
  "error"      => 0,
  "snp"        => 0,
  "splice"     => 0,
  "deletion"   => 0,
  "insertion"  => 0,
  "chimera"    => 0,
  "exon"       => 0,
  "transcript" => 0,
);

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    checker => $args{checker},
    args => \%args,
  }, $class;

  my @check_types;

  my $false_positives_file = $args{false_positives_file};
  my $false_negatives_file = $args{false_negatives_file};
  my $true_positives_file = $args{true_positives_file};

  if($args{check} eq 'all') {
    push @check_types, @{$self->allCheckableEvents};
  } else {
    push @check_types, @{$args{check}};
  }

  foreach my $check_type (@check_types) {
    $self->addStats($check_type,
      CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents($check_type),
        false_positives_file => defined $false_positives_file? "$false_positives_file-$check_type.log" : undef,
        false_negatives_file => defined $false_negatives_file? "$false_negatives_file-$check_type.log" : undef,
        true_positives_file => defined $true_positives_file? "$true_positives_file-$check_type.log" : undef,
        print_header => sub { $self->checker->getEvents($check_type)->printHeader(@_) },
        print_element => sub { $self->checker->getEvents($check_type)->printEvent(@_) },
      )
    );
  }

  #$self->_init(@_);

  return $self;
}

sub run {
  my $self = shift;
  $self->_init(%{$self->{args}});
  # Now we close stat files outputs
  foreach my $stats (values %{$self->{stats}}) {
    $stats->closeOutputs();
  }
}

sub _init {
  my $self = shift;
}

sub checker {
  my $self = shift;
  return $self->{checker};
}

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if(defined $check_events{$event_type}) {
    return $check_events{$event_type};
  } else {
    carp "Unkown event: $event_type";
  }
}

sub allCheckableEvents {
  my $self = shift;
  my @checkable_events;
  foreach my $event_type (keys %check_events) {
    if($self->canCheck($event_type)) {
      push(@checkable_events,$event_type);
    }
  }
  return \@checkable_events;
}

sub allCheckedEvents {
  my $self = shift;
  my @checked_events;
  foreach my $event_type (keys %check_events) {
    push(@checked_events,$event_type) if defined $self->getStats($event_type);
  }
  return \@checked_events;
}

sub getStats {
  my $self = shift;
  my $event_type = shift;
  return $self->{stats}->{$event_type};
}

sub addStats {
  my $self = shift;
  my $event_type = shift;
  my $stats = shift;
  $self->{stats}->{$event_type} = $stats;
}

1;
