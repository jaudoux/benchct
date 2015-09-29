use strict;
use warnings;
package CracTools::BenchCT::Analyzer::GTF;
# ABSTRACT: Analyze results from a transcript reconstruction over the simulated data handled by BenchCT

use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) || $event_type =~ /^(exon|transcript)$/) {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;

  my $gtf_file = $args{file};
  my $gtf_it = CracTools::Utils::gffFileIterator($gtf_file,'gtf');

  my $prev_transcript_id;
  my $transcript_id;
  my @transcript_exons;
  
  while (my $gtf_line = $gtf_it->()) {
    next if $gtf_line->{feature} ne 'exon';

    $transcript_id = $gtf_line->{attributes}->{transcript_id};
    next if !defined $transcript_id;

    if(!defined $prev_transcript_id || $prev_transcript_id ne $transcript_id) {
      if(defined $prev_transcript_id && defined $self->getStats('transcript')) {
        $self->_checkTranscript(\@transcript_exons,$prev_transcript_id);
      }
      @transcript_exons = ();
      $prev_transcript_id = $transcript_id;
    }

    # Remember that checker is 0-based
    my $true_exon = $self->checker->isTrueExon(
      $gtf_line->{chr},
      $gtf_line->{start}-1,
      $gtf_line->{end}-1,
      CracTools::Utils::convertStrand($gtf_line->{strand}),
    );

    if($true_exon) {
      $self->getStats('exon')->addTruePositive(id => $true_exon) if defined $self->getStats('exon');
      push @transcript_exons, $true_exon - 1;
    } else {
      $self->getStats('exon')->addFalsePositive(out_string => join("\t",
          $gtf_line->{chr},
          $gtf_line->{start},
          $gtf_line->{end},
          $gtf_line->{strand},
        ),
      ) if defined $self->getStats('exon');
    }
  }
  # Check the last transcript
  $self->_checkTranscript(\@transcript_exons,$transcript_id) if defined $transcript_id 
  && defined $self->getStats('transcript');
}

sub _checkTranscript {
  my $self = shift;
  my ($transcript_exons,$transcript_id) = @_;

  my $true_transcript = $self->checker->isTrueTranscript($transcript_exons);
  if($true_transcript) {
    $self->getStats('transcript')->addTruePositive(id => $true_transcript);
  } else {
    $self->getStats('transcript')->addFalsePositive();
  }
}

1;
