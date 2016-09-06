use strict;
use warnings;
package CracTools::BenchCT::Events::Chimera;
# ABSTRACT: A collection of events of type 'mutation' (snp,indel) on the genome

use parent 'CracTools::BenchCT::Events';

use Set::IntervalTree;

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

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
  my $chim_key = $self->getEvent($event_id);
  print $fh join("\t",$event_id, splitChimKey($chim_key)),"\n";
}

sub isTrueChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;

  my @overlapping_chimeras = @{$self->getOverlappingChimeras(
    $chr1,$pos1,$strand1,$chr2,$pos2,$strand2
  )};

  # push @overlapping_chimeras, @{$self->getOverlappingChimeras(
  #   $chr1,$pos1,$strand1*-1,$chr2,$pos2,$strand2*-1
  # )};
  push @overlapping_chimeras, @{$self->getOverlappingChimeras(
    $chr2,$pos2,$strand2*-1,$chr1,$pos1,$strand1*-1
  )};

  my $best_chim_id = -1;
  my $chim_overlap_dist = -1;
  foreach my $chim_id (@overlapping_chimeras) {
    my $chim_key = $self->getEvent($chim_id);

    my ($c_c1,$c_p1,$c_s1,$c_c2,$c_p2,$c_s2) = splitChimKey($chim_key);
    my $dist1 = overlapDistance($pos1,$pos2,$c_p1,$c_p2);
    my $dist2 = overlapDistance($pos2,$pos1,$c_p1,$c_p2);

    if($dist1 <= $dist2) {
      if($dist1 <= $self->threshold && ($chim_overlap_dist == -1 || $dist1 < $chim_overlap_dist)) {
        $chim_overlap_dist = $dist1;
        $best_chim_id = $chim_id;
      }
    } else {
      if($dist2 <= $self->threshold && ($chim_overlap_dist == -1 || $dist2 < $chim_overlap_dist)) {
        $chim_overlap_dist = $dist2;
        $best_chim_id = $chim_id;
      }
    }
  }
  return $best_chim_id + 1;
}

sub addChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;

  # FIXME We should check that the chimera does not already exists!
  my $chim_key = getChimKey($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
  my $chim_id  = $self->addEvent($chim_key); # Increment the number of chimeras
  my $class    = getClass($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);

  $self->{nb_chimeras}++;

  # CLASS 1 CHIMERAS, we use arrays as index
  if($class == 1) {
    my $key1 = $chr1."@".$strand1;
    my $key2 = $chr2."@".$strand2;
    if(!defined $self->{indexes}{$class}{$key1}{$key2}) {
      $self->{indexes}{$class}{$key1}{$key2} = [];
    }
    push(@{$self->{indexes}{$class}{$key1}{$key2}},$chim_id);
  # CLASS 2,3,4 CHIMERAS, we use intervalTrees as index
  } else {
    my ($pos_min,$pos_max) = getChimeraPosMinMax($pos1, $pos2);
    if(!defined $self->{indexes}{$class}{$chr1}{$strand1}) {
      $self->{indexes}{$class}{$chr1}{$strand1} = Set::IntervalTree->new;
    }
    $pos_max++ if $pos_min == $pos_max; # Avoid null width intervals
    $self->{indexes}{$class}{$chr1}{$strand1}->insert($chim_id);
  }

  return $chim_id;
}

sub getOverlappingChimeras {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;
  my $class    = getClass($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);

  my @overlapping_chimeras = ();

  # CLASS 1 CHIMERAS
  if($class == 1) {
    my $key1 = $chr1."@".$strand1;
    my $key2 = $chr2."@".$strand2;
    foreach my $chim_id (@{$self->{indexes}{$class}{$key1}{$key2}}) {
      push(@overlapping_chimeras,$chim_id);
    }
  # CLASS 2,3,4 CHIMERAS
  } else {
    my ($pos_min,$pos_max) = getChimeraPosMinMax($pos1, $pos2);
    $pos_max++ if $pos_min == $pos_max; # Avoid null width intervals

    my $hits;
    if(defined $self->{indexes}{$class}{$chr1}{$strand1}) {
      $hits = $self->{indexes}{$class}{$chr1}{$strand1}->fetch($pos_min,$pos_max);
      foreach my $chim_id (@{$hits}) {
          push(@overlapping_chimeras,$chim_id);
      }
    }
  }
  return \@overlapping_chimeras;
}

sub getChimeraPosMinMax {
  my ($pos1,$pos2) = @_;
  ($pos1,$pos2) = ($pos2,$pos1) if $pos2 < $pos1;
  return ($pos1,$pos2);
}

sub getClass {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;
  if ($chr1 eq $chr2){
    if ($strand1 == $strand2){
      if (($strand1 == 1 && $pos1 > $pos2) || ($strand1 == -1 && $pos2 > $pos1)) {
        return 3;
      } else {
        return 2;
      }
    } else {
      return 4;
    }
  } else {
    return 1;
  }
}

sub getChimKey {
  return join("@", @_);
}

sub splitChimKey {
  my $chim_key = shift;
  return split "@", $chim_key;
}

sub overlapDistance {
  my ($pos1_A,$pos2_A,$pos1_B,$pos2_B) = @_;
  return abs($pos1_A - $pos1_B) + abs($pos2_A - $pos2_B);
}

1;
