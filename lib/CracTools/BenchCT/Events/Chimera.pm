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
  $self->{chim_ids} = {};
  
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
  my $id = $self->addEvent(); # Increment the number of chimeras
  $self->{chim_ids}{$chim->getKey} = $id;
  $self->overlapStructure->addChimera($chim);
}

sub isTrueChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_; 
  my $chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  my @overlapping_chimeras = @{$self->overlapStructure->getOverlappingChimeras($chim)};
  foreach my $overlap_chim (@overlapping_chimeras) {
    if((abs($overlap_chim->pos1 - $pos1) + abs($overlap_chim->pos2 - $pos2)) <= $self->threshold ||
       (abs($overlap_chim->pos2 - $pos1) + abs($overlap_chim->pos1 - $pos2)) <= $self->threshold) {
        return $self->{chim_ids}{$overlap_chim->getKey} + 1;
     }
  }
  return 0;
}

1;
