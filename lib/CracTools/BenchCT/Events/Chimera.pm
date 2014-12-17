use strict;
use warnings;
package CracTools::BenchCT::Events::Chimera;
# ABSTRACT: A collection of events of type 'mutation' (snp,indel) on the genome

use parent 'CracTools::BenchCT::Events';

use CracTools::ChimCT::Chimera;
use CracTools::ChimCT::OverlapStructure;

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

  $self->{overlap_structure} = CracTools::ChimCT::OverlapStructure->new(max_overlapping_distance => $self->threshold);
  
  return $self;
}

sub overlapStructure {
  my $self = shift;
  return $self->{overlap_structure};
}

sub addChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_; 
  my $chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  $self->overlapStructure->addChimera($chim);
  $self->addEvent(); # Increment the number of splices
}

sub isTrueChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_; 
  my $chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  return $self->overlapStructure->nbOverlappingChimeras($chim);
}

1;
