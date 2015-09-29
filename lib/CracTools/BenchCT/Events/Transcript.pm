use strict;
use warnings;
package CracTools::BenchCT::Events::Transcript;
# ABSTRACT: A collection of events of type 'transcript' on the genome

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

  $self->{transcripts} = {};
  
  return $self;
}

sub printEvent {
  my ($self,$fh,$id) = @_;
  my $exons = $self->getEvent($id);
  print $fh join("\t",$id,join(",",@{$exons})),"\n";
}

sub addTranscript {
  my $self = shift;
  my @exons = @_; 

  my ($nb_match,$transcript_id) = $self->isTrueTranscript(\@exons);

  # Check if this transcript already exists
  if($transcript_id && $nb_match == @exons) {
    return $transcript_id - 1;
  # If not we create a new entry
  } else {
    my $id = $self->addEvent(\@exons); # Increment the number of exons

    # We place the transcript_id in the exons hash to access the transcript
    # structure from any of its exons
    foreach my $exon (@exons) {
      push @{$self->{exons}->{$exon}}, $id;
    }

    $self->{transcripts}->{$id} = \@exons;
    return $id;
  }
}

sub isTrueTranscript {
  my $self = shift;
  my $exons = shift; 
  return 0 if !defined $exons || @{$exons} == 0;
  my $first_exon = $exons->[0];
  my $potential_transcripts = $self->{exons}->{$first_exon};
  my $best_match;
  my $nb_exons_best_match = 0;
  foreach my $transcript (@{$potential_transcripts}) {
    #next if @{$self->{transcripts}->{$transcript}} != @{$exons};
    my $nb_match_exons = 0;
    foreach my $exon (@{$exons}) {
      if($exon ~~ @{$self->{transcripts}->{$transcript}}) {
        $nb_match_exons++;
      } else {
        last;
      }
    }
    if($nb_match_exons == @{$exons} && $nb_match_exons > $nb_exons_best_match) {
      $best_match          = $transcript;
      $nb_exons_best_match = $nb_match_exons;
    }
  }
  if(defined $best_match) {
    return ($nb_exons_best_match,$best_match + 1);
  } else {
    return 0;
  }
}

1;
