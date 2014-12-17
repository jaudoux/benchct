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
);

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    checker => $args{checker},
  }, $class;

  my @check_types;

  if($args{check} eq 'all') {
    push @check_types, @{$self->allCheckableEvents};
  } else {
    push @check_types, @{$args{check}};
  }

  foreach my $check_type (@check_types) {
    if($check_type eq 'mapping' && $self->canCheck($check_type)) {
      $self->addStats($check_type,
        CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbReads())
      );
    } elsif($check_type eq 'error' && $self->canCheck($check_type)) {
      $self->addStats($check_type,
        CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbErrors())
      );
    } elsif($check_type =~ /^(snp|splice|deletion|insertion)$/ && $self->canCheck($check_type)) {
      $self->addStats($check_type,
        CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents($check_type))
      );
    }
  }

  $self->_init(@_);

  return $self;
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
