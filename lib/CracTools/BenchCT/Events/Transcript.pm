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

sub addTranscript {
  my $self = shift;
  my @exons = @_; 

  my $id = $self->addEvent(\@exons); # Increment the number of exons

  # We place the transcript_id in the exons hash to access the transcript
  # structure from any of its exons
  foreach my $exon (@exons) {
    push @{$self->{exons}->{$exon}}, $id;
  }

  $self->{transcripts}->{$id} = \@exons;

  return $id;
}

sub isTrueTranscript {
  my $self = shift;
  my $exons = shift; 
  return 0 if !defined $exons || @{$exons} == 0;
  my $first_exon = $exons->[0];
  my $potential_transcripts = $self->{exons}->{$first_exon};
  foreach my $transcript (@{$potential_transcripts}) {
    #print STDERR "Try matching transcript $transcript: ";
    #print STDERR join ",",@{$self->{transcripts}->{$transcript}},"\n";
    #next if @{$self->{transcripts}->{$transcript}} != @{$exons};
    my $nb_match_exons = 0;
    foreach my $exon (@{$exons}) {
      if($exon ~~ @{$self->{transcripts}->{$transcript}}) {
        $nb_match_exons++;
      } else {
        last;
      }
    }
    if($nb_match_exons == @{$exons}) {
      #print STDERR "Transcript matched\n";
      return $transcript + 1;
    }
  }
  return 0;
}

1;
