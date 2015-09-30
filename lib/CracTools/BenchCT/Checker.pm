use strict;
use warnings;
package CracTools::BenchCT::Checker;
# ABSTRACT: Benchmarking CRAC results

use CracTools::Utils;
use CracTools::BenchCT::Const;
use CracTools::BenchCT::Utils;
use CracTools::BenchCT::Events::Alignment;
use CracTools::BenchCT::Events::Error;
use CracTools::BenchCT::Events::Mutation;
use CracTools::BenchCT::Events::SNP;
use CracTools::BenchCT::Events::Splice;
use CracTools::BenchCT::Events::Chimera;
use CracTools::BenchCT::Events::Exon;
use CracTools::BenchCT::Events::Transcript;
use Data::Dumper; #for debugging
use Carp;

=head1 DOCUMENTATION

CracTools::BenchCT::Checker is 0-based!!!
CracTools::BenchCT::Checker is not strand specific by default!

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $is_stranded = $args{is_stranded};
  $is_stranded = 0 if !defined $is_stranded;

  if(defined $args{err_file} && !defined $args{bed_file}) {
    # We need the bed file to check errors
    croak "Missing 'bed_file' to create a conversion table from read name to read id in order to check errors";
  }

  my $self = bless {
    is_stranded         => $is_stranded,
    vcf_file            => $args{vcf_file}, # contains mutations
    bed_file            => $args{bed_file}, # contains alignments
    err_file            => $args{err_file}, # contains errors
    junction_bed_file   => $args{junction_bed_file},
    chimera_tsv_file    => $args{chimera_tsv_file},
    gtf_file            => $args{gtf_file},
    read_ids            => {}, # We could avoid that somehow by placing in read names their mapping position and error informations
    verbose             => defined $args{verbose}? $args{verbose} : 0,
    events              => {
      mapping => CracTools::BenchCT::Events::Alignment->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_MAPPING,
      ),
      snp       => CracTools::BenchCT::Events::SNP->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_SNP,
        crac_index_conf => $args{crac_index_conf},
      ),
      insertion => CracTools::BenchCT::Events::Mutation->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_INS,
      ),
      deletion  => CracTools::BenchCT::Events::Mutation->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_DEL,
      ),
      splice    => CracTools::BenchCT::Events::Splice->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_SPLICE,
      ),
      chimera   => CracTools::BenchCT::Events::Chimera->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_CHIMERA,
      ),
      exon   => CracTools::BenchCT::Events::Exon->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_EXON,
      ),
      transcript   => CracTools::BenchCT::Events::Transcript->new(
        verbose => $args{verbose},
        threshold => $CracTools::BenchCT::Const::THRESHOLD_TRANSCRIPT,
      ),
    },
  }, $class;

  $self->_init();

  return $self;
}

sub _init {
  my $self = shift;

  my $max_read_length = 0;

  # Read bed file for alignements
  if(defined $self->{bed_file}) {
    print STDERR "[checker] Reading bed file\n" if $self->verbose;
    
    # First we read the bed files and we record seek positions of each read
    my $bed_it = CracTools::Utils::getFileIterator(file =>$self->{bed_file},
      parsing_method => \&CracTools::BenchCT::Utils::parseGSBedLine,
    );

    while(my $bed_line = $bed_it->()) {
      # Add the alignment
      my $id = $self->getEvents('mapping')->addAlignment($bed_line);
      # Create an entry to the conversion table
      $self->{read_ids}->{$bed_line->{name}} = $id;
      # Update the maximum read_length
      my $read_length = 0;
      map { $read_length += $_->{size} } @{$bed_line->{blocks}};
      $max_read_length = $read_length if $read_length > $max_read_length;
    }

    print STDERR "[checker] ".scalar $self->nbEvents('mapping')." alignment(s) read\n" if $self->verbose;

    # Read errors only if we have a bed file
    if(defined $self->{err_file}) {

      # We create the error library set
      $self->{events}->{error} = CracTools::BenchCT::Events::Error->new(
        nb_reads      => $self->nbEvents('mapping'),
        max_length    => $max_read_length,
        verbose       => $self->{verbose},
      );

      print STDERR "[checker] Reading error file\n" if $self->verbose;

      # Now we read the error file
      my $err_it = CracTools::Utils::getFileIterator(file =>$self->{err_file},
        parsing_method => \&CracTools::BenchCT::Utils::parseErrLine,
      );

      while(my $err_line = $err_it->()) {
        if($err_line->{read_id} >= $self->nbEvents('mapping')) {
          croak "There is more read in 'err_file' than in the 'bed_file' alignment file";
        }
        $self->getEvents('error')->addError($err_line->{read_id},$err_line->{pos});
      }

      print STDERR "[checker] ".$self->nbEvents('error')." error(s) read\n" if $self->verbose;
    }

  }

  # Read junction bed
  if(defined $self->{junction_bed_file}) {
    print STDERR "[checker] Reading junction bed file\n" if $self->verbose;
    # First we read the bed files and we record seek positions of each read
    my $bed_it = CracTools::Utils::bedFileIterator($self->{junction_bed_file});

    while(my $bed_line = $bed_it->()) {
      $self->getEvents('splice')->addSplice(
        $bed_line->{chr},
        $bed_line->{start},
        $bed_line->{end} - $bed_line->{start},
        $self->isStranded? CracTools::Utils::convertStrand($bed_line->{strand}) : 1,
      );
    }
    print STDERR "[checker] ".scalar $self->nbEvents('splice')." splice(s) read\n" if $self->verbose;
  }

  # Read chimera file
  if(defined $self->{chimera_tsv_file}) {
    print STDERR "[checker] Reading chimera file\n" if $self->verbose;
    # First we read the bed files and we record seek positions of each read
    my $chim_it = CracTools::Utils::getFileIterator(file =>$self->{chimera_tsv_file},
      parsing_method => \&CracTools::BenchCT::Utils::parseChimeraLine,
    );

    while(my $chim_line = $chim_it->()) {
      $self->getEvents('chimera')->addChimera($chim_line->{chr1},
        $chim_line->{pos1},
        $chim_line->{strand1},
        $chim_line->{chr2},
        $chim_line->{pos2},
        $chim_line->{strand2},
      );
    }
    print STDERR "[checker] ".$self->nbEvents('chimera')." chimera(s) read\n" if $self->verbose;
  }

  # Read vcf file (snps, indels)
  if(defined $self->{vcf_file}) {
    print STDERR "[checker] Reading vcf file\n" if $self->verbose;
    # Secondly we read the info file to register mutations
    my $vcf_it = CracTools::Utils::vcfFileIterator($self->{vcf_file});

    # store tag ids for each mutation in an Interval Tree (one for each type of mutations
    while(my $vcf_line = $vcf_it->()) {
      my $ref_length = length $vcf_line->{ref};
      foreach my $alt (@{$vcf_line->{alt}}) {
        my $alt_length = length $alt;
        my $mutation_type = $vcf_line->{type};
        next if $alt eq 'N' && $vcf_line->{ref} eq 'N';
        # It is a substitution
        if($ref_length == $alt_length) {
          # We shift the pos if the reference has more than one nucleotide
          my $id = $self->getEvents('snp')->addMutation($vcf_line->{chr},$vcf_line->{pos} + $ref_length - 2,$alt);
        # This is a deletion
        } elsif($ref_length > $alt_length) {
          $self->getEvents('deletion')->addMutation($vcf_line->{chr},$vcf_line->{pos} + $alt_length - 2,$ref_length - $alt_length);
        # This is an insertion
        } else {
          $self->getEvents('insertion')->addMutation($vcf_line->{chr},$vcf_line->{pos} + $ref_length - 2,$alt_length - $ref_length);
        }
      }
    }
    print STDERR "[checker] ".$self->nbEvents('snp')." SNP(s) read\n" if $self->verbose;
    print STDERR "[checker] ".$self->nbEvents('insertion')." insertion(s) read\n" if $self->verbose;
    print STDERR "[checker] ".$self->nbEvents('deletion')." deletion(s) read\n" if $self->verbose;
  }

  # Read GTF annotations
  if(defined $self->{gtf_file}) {
    print STDERR "[checker] Reading gtf file\n" if $self->verbose;
    my $gtf_it = CracTools::Utils::gffFileIterator($self->{gtf_file},'gtf');
    my $prev_transcript_id;
    my $transcript_id;
    my @transcript_exons;
    while(my $gtf_line = $gtf_it->()) {
      next if $gtf_line->{feature} ne 'exon';

      $transcript_id = $gtf_line->{attributes}->{transcript_id};
      next if !defined $transcript_id;
      
      # If this exon belong to a different transcript we add the previous one
      if(!defined $prev_transcript_id || $prev_transcript_id ne $transcript_id) {
        if(defined $prev_transcript_id) {
          my $id = $self->getEvents('transcript')->addTranscript(@transcript_exons);
        }
        @transcript_exons = ();
        $prev_transcript_id = $transcript_id;
      }
      # Remove chr prefix, in case of
      $gtf_line->{chr} =~ s/^chr//;
      # Now we add the current exon
      my $exon_id = $self->getEvents('exon')->addExon(
        $gtf_line->{chr},
        $gtf_line->{start}-1,
        $gtf_line->{end}-1,
        $self->isStranded? CracTools::Utils::convertStrand($gtf_line->{strand}) : 1,
      );
      # And push it to the current transcript
      push @transcript_exons, $exon_id;
    }
    # We add the last transcript
    if(@transcript_exons) {
      $self->getEvents('transcript')->addTranscript(@transcript_exons);
    }
    print STDERR "[checker] ".$self->nbEvents('exon')." Exons(s) read\n" if $self->verbose;
    print STDERR "[checker] ".$self->nbEvents('transcript')." Transcript(s) read\n" if $self->verbose;
  }
}

=head2 isStranded

return true is the CracTools::BenchCT::Checker used stranded simulated data, and therfore is
carreful about strand given by each tools for alignments or events.

=cut

sub isStranded {
  my $self = shift;
  return $self->{is_stranded};
}

=head2 verbose 

=cut

sub verbose {
  my $self = shift;
  return $self->{verbose};
}

=head2 nbErrors

Return the number of reads in the simulated data

=cut

sub nbErrors {
  my $self = shift;
  return $self->{nb_errors};
}

=head2 nbEvents($type)

Return the number of events of a specific type (snp, ins, del, splice, chimera)

=cut

sub nbEvents {
  my $self = shift;
  my $type = shift;
  my $events = $self->getEvents($type);
  if(defined $events) {
    return $events->nbEvents;
  } else {
    return 0;
  }
}

=head2 getReadId($read_name)

Given a read_name, this method return the read id (from 0 to nb_reads).

=cut

sub getReadId {
  my $self = shift;
  my $read_name = shift;
  my $read_id = $read_name =~ /^\d+$/? $read_name : $self->{read_ids}->{$read_name};
  return $read_id;
}

=head2 getEvents($type)

Return the CracTools::BenchCT::Events of a specific type (snp, ins, del, splice, chimera)

=cut

sub getEvents {
  my $self = shift;
  my $type = shift;
  return $self->{events}->{$type};
}

=head2 isTrueInsertion($chr,$pos,$length)

=cut

sub isTrueInsertion {
  my $self = shift;
  return $self->getEvents('insertion')->isTrueMutation(@_);
}

=head2 isTrueSNP($chr,$pos,$nuc)

=cut

sub isTrueSNP {
  my $self = shift;
  return $self->getEvents('snp')->isTrueMutation(@_);
}

=head2 isTrueDeletion$chr,$pos,$length)

=cut

sub isTrueDeletion {
  my $self = shift;
  return $self->getEvents('deletion')->isTrueMutation(@_);
}

=head2 isTrueMutation($type,$chr,$pos,$length)

Return true is there is a mutation at this position

=cut

sub isTrueMutation {
  my $self = shift;
  my ($type,$chr,$pos,$length) = @_;
  my $snp_events = $self->getEvents($type);
  #if($type eq 'snp') {
  #  print STDERR "NB BITS: ".$snp_events->isTrueMutation($chr,$pos,$length)."\n";
  #}
  return $snp_events->isTrueMutation($chr,$pos,$length);
}

=head2 isTrueSplice($chr,$start,$length,$strand)

Return true is there is a mutation at this position

=cut

sub isTrueSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_;
  $strand = 1 if !$self->isStranded();
  my $splice_events = $self->getEvents('splice');
  return $splice_events->isTrueSplice($chr,$start,$length,$strand);
}

=head2 isTrueExon($chr,$start,$end,$strand)

Return true if the specificed exon match an annotated exon

=cut

sub isTrueExon {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_;
  $strand = 1 if !$self->isStranded();
  my $exon_events = $self->getEvents('exon');
  return $exon_events->isTrueExon($chr,$start,$end,$strand);
}

=head2 isTrueTranscript($exons)

Return true is the specified transcripts matches a annotated transcript

=cut

sub isTrueTranscript {
  my $self = shift;
  my $exons = shift;
  my $transcript_events = $self->getEvents('transcript');
  return $transcript_events->isTrueTranscript($exons);
}

=head2 isTrueChimera($chr1,$pos1,$strand1,$chr2,$pos2,$strand2)

Return true if there is chimeric junction for these coordinates.

=cut

sub isTrueChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;
  my $chimera_events = $self->getEvents('chimera');
  return $chimera_events->isTrueChimera($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
}

=head2 isTrueError 

=cut

sub isTrueError {
  my $self = shift;
  my ($read_name, $pos) = @_;
  my $read_id = $self->getReadId($read_name);

  # Check if the read id exists
  return 0 if !defined $read_id;

  my $error_events = $self->getEvents('error');
  return $error_events->isTrueError($read_id,$pos);
}

=head2 isGoodAlignment

  [read_name]
  [ref_name]
  [pos_start] 0 based
  [ref_start] 0 based

Return true is the alignment of this read is good.

=cut

sub isGoodAlignment {
  my ($self,$read_name,$ref_name,$pos_start,$ref_start,$strand) = @_;

  my $read_id = $self->getReadId($read_name);

  # Check if the read id exists
  return 0 if !defined $read_id;

  my $alignment_events = $self->getEvents('mapping');
  return $alignment_events->isGoodAlignment($read_id,$ref_name,$pos_start,$ref_start,$strand);
}

1;
