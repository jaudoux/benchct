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

  my $transcript_id;
  my @transcript_exons;
  #my $corrupted_transcript = 0;
  
  while (my $gtf_line = $gtf_it->()) {
    next if $gtf_line->{feature} ne 'exon';
    next if !defined $gtf_line->{attributes}->{transcript_id};

    if(!defined $transcript_id || $transcript_id ne $gtf_line->{attributes}->{transcript_id}) {
      $self->_checkTranscript(\@transcript_exons,$transcript_id) if defined $transcript_id;
      #$corrupted_transcript = 0;
      @transcript_exons = ();
      $transcript_id = $gtf_line->{attributes}->{transcript_id};
    }

    # Remember that checker is 0-based
    my $true_exon = $self->checker->isTrueExon(
      $gtf_line->{chr},
      $gtf_line->{start}-1,
      $gtf_line->{end}-1,
      CracTools::Utils::convertStrand($gtf_line->{strand}),
    );

    if($true_exon) {
      #print STDERR "Found exon: ", ($true_exon-1) ,"\n";
      $self->getStats('exon')->addTruePositive(id => $true_exon);
      push @transcript_exons, $true_exon - 1;
    } else {
      #print STDERR "Miss exon: ",join("\t",
      #    $gtf_line->{chr},
      #    $gtf_line->{start},
      #    $gtf_line->{end},
      #    $gtf_line->{strand}),"\n";
      #$corrupted_transcript = 1;
      $self->getStats('exon')->addFalsePositive(out_string => join("\t",
          $gtf_line->{chr},
          $gtf_line->{start},
          $gtf_line->{end},
          $gtf_line->{strand},
        ),
      );
    }
  }
  # Check the last transcript
  $self->_checkTranscript(\@transcript_exons,$transcript_id) if defined $transcript_id;
}

sub _checkTranscript {
  my $self = shift;
  my ($transcript_exons,$transcript_id) = @_;

  my $true_transcript = $self->checker->isTrueTranscript($transcript_exons);
  if($true_transcript) {
    $self->getStats('transcript')->addTruePositive(id => $true_transcript);
  } else {
    $self->getStats('transcript')->addFalsePositive(out_string => $transcript_id);
  }
}

1;
