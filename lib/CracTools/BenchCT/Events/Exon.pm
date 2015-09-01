use strict;
use warnings;
package CracTools::BenchCT::Events::Exon;
# ABSTRACT: A collection of events of type 'exon' on the genome

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

sub addExon {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_; 

  # First we look if this exon already exists using a 0 threshold
  my $id = $self->_foundExon($chr,$start,$end,$strand,0);
  if($id) {
    #print STDERR "No new exon added (found exon ",($id -1),")\n";
    return $id - 1;
  } else {
    $id = $self->addEvent($end-$start); # Increment the number of exons
    #print STDERR "New exon added (id: $id)\n";
    $self->intervalQuery->addInterval($chr,$start,$end,$strand,$id);
    return $id;
  }
}

sub isTrueExon {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_; 
  return $self->_foundExon($chr,$start,$end,$strand,$self->threshold);
}

sub _foundExon {
  my $self = shift;
  my ($chr,$start,$end,$strand,$threshold) = @_; 

  # Check if we have already seen this exon
  # The 1 at the we means that we use a "windowed" query
  # TODO there is a bug in Set::IntervalTree windowed query, we are no longer using them
  my @matching_exons = @{$self->intervalQuery->fetchByRegion($chr,$start-$threshold,$end+$threshold,$strand,0)};
  my $length = $end - $start;
  foreach my $exon (@matching_exons) {
    my $exon_length = $self->getEvent($exon);
    if(abs($exon_length - $length) <= $threshold) {
      return $exon + 1;
    }
  }

  return 0;
}

1;
