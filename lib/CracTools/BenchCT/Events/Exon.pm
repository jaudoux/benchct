use strict;
use warnings;
package CracTools::BenchCT::Events::Exon;
# ABSTRACT: A collection of events of type 'exon' on the genome

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

sub printEvent {
  my ($self,$fh,$id) = @_;
  my $exon = $self->getEvent($id);
  print $fh join("\t",$id,split('-',$exon)),"\n";
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
    return $id - 1;
  } else {
    $id = $self->addEvent("$start-$end"); # Increment the number of exons
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
  my $best_overlap = 0;
  my $best_exon;
  foreach my $exon (@matching_exons) {
    my ($exon_start,$exon_end) = split('-',$self->getEvent($exon));
    my $overlap_length = min($exon_end,$end) - max($exon_start,$start) + 1;
    my $total_length   = max($exon_end,$end) - min($exon_start,$start) + 1;
    my $overlap_cover  = $overlap_length/$total_length;
    if(($total_length - $overlap_length) <= $threshold && $overlap_cover > $best_overlap) {
      $best_exon    = $exon;
      $best_overlap = $overlap_cover;
      last if $best_overlap == 1;
    }
  }
  
  if(defined $best_exon) {
    return $best_exon + 1;
  } else {
    return 0;
  }
}

1;
