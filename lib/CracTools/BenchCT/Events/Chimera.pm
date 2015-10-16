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

sub printHeader {
  my $self = shift;
  my $fh = shift;
  print $fh join("\t",qw(event_id chr1 pos1 strand1 chr2 pos2 strand2)),"\n";
}

sub printEvent {
  my $self = shift;
  my $fh = shift;
  my $event_id = shift;
  my $chim = $self->getEvent($event_id);
  print $fh join("\t",$event_id,
    $chim->chr1,
    $chim->pos1,
    $chim->strand1,
    $chim->chr2,
    $chim->pos2,
    $chim->strand2,
  ),"\n";
}

sub overlapStructure {
  my $self = shift;
  return $self->{overlap_structure};
}

sub addChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_; 
  my $chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  my $id = $self->addEvent($chim); # Increment the number of chimeras
  $self->{chim_ids}{$chim->getKey} = $id;
  $self->overlapStructure->addChimera($chim);
}

sub isTrueChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_; 
  my $chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  my $chim2 = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1*-1,$chr2,$pos2,$strand2*-1);
  my @overlapping_chimeras = @{$self->overlapStructure->getOverlappingChimeras($chim)};
  push @overlapping_chimeras, @{$self->overlapStructure->getOverlappingChimeras($chim2)};
  my $chim_id = -1;
  my $chim_overlap_dist = -1;
  foreach my $overlap_chim (@overlapping_chimeras) {
    my $dist1 = abs($overlap_chim->pos1 - $pos1) + abs($overlap_chim->pos2 - $pos2);
    my $dist2 = abs($overlap_chim->pos2 - $pos1) + abs($overlap_chim->pos1 - $pos2);
    if($dist1 <= $dist2) {
      if($dist1 <= $self->threshold && ($chim_overlap_dist == -1 || $dist1 < $chim_overlap_dist)) {
        $chim_overlap_dist = $dist1;
        $chim_id = $self->{chim_ids}{$overlap_chim->getKey};
      }
    } else {
      if($dist2 <= $self->threshold && ($chim_overlap_dist == -1 || $dist2 < $chim_overlap_dist)) {
        $chim_overlap_dist = $dist2;
        $chim_id = $self->{chim_ids}{$overlap_chim->getKey};
      }
    }
  }
  return $chim_id+1;
}

1;
