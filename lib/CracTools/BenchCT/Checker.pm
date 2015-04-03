use strict;
use warnings;
package CracTools::BenchCT::Checker;
# ABSTRACT: Benchmarking CRAC results

use CracTools::Utils;
use CracTools::BenchCT::Const;
use CracTools::BenchCT::Utils;
use CracTools::BenchCT::Events::Mutation;
use CracTools::BenchCT::Events::SNP;
use CracTools::BenchCT::Events::Splice;
use CracTools::BenchCT::Events::Chimera;
use Data::Dumper; #for debugging
use Carp;

=head1 DOCUMENTATION

CracTools::BenchCT::Checker is 0-based!!!

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $is_stranded = $args{is_stranded};
  $is_stranded = 0 if !defined $is_stranded;

  my $self = bless {
    is_stranded         => $is_stranded,
    #info_file           => $args{info_file}, # contains mutations
    vcf_file            => $args{vcf_file}, # contains mutations
    bed_file            => $args{bed_file}, # contains alignements
    bed_fh              => defined $args{bed_file}? CracTools::Utils::getReadingFileHandle($args{bed_file}) : undef,
    err_file            => $args{err_file}, # contains errors
    err_fh              => defined $args{err_file}? CracTools::Utils::getReadingFileHandle($args{err_file}) : undef,
    junction_bed_file   => $args{junction_bed_file},
    chimera_tsv_file    => $args{chimera_tsv_file},
    read_ids            => {}, # TODO avoid this conversion hash and use a interval tree to seek for bed line that could correspond to a read alignement, or we could name using there read_id number...
    bed_seek_pos        => [],
    nb_reads            => 0,
    err_seek_pos        => [],
    verbose             => defined $args{verbose}? $args{verbose} : 0,
    events              => {
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
    },
  }, $class;

  $self->_init();

  return $self;
}

sub _init {
  my $self = shift;

  # Read bed file for alignements
  if(defined $self->{bed_file}) {
    print STDERR "[checker] Reading bed file\n" if $self->verbose;
    
    # First we read the bed files and we record seek positions of each read
    my $bed_it = CracTools::Utils::getFileIterator(file =>$self->{bed_file},
      parsing_method => \&CracTools::BenchCT::Utils::parseGSBedLineLite,
    );

    my $id = 0;
    while(my $bed_line = $bed_it->()) {
      # Set number to read ids
      $self->{read_ids}->{$bed_line->{name}} = $id;
      # Set seek pos of each bed alignement according to the read id
      $self->{bed_seek_pos}[$id] = $bed_line->{seek_pos};
      $id++;
    }

    $self->{nb_reads} = $id;
    print STDERR "[checker] ".scalar $self->nbReads." read(s) parsed\n" if $self->verbose;
  }

  # Read junction bed
  if(defined $self->{junction_bed_file}) {
    print STDERR "[checker] Reading junction bed file\n" if $self->verbose;
    # First we read the bed files and we record seek positions of each read
    my $bed_it = CracTools::Utils::bedFileIterator($self->{junction_bed_file});

    while(my $bed_line = $bed_it->()) {
      $self->getEvents('splice')->addSplice($bed_line->{chr},$bed_line->{start},$bed_line->{end} - $bed_line->{start},CracTools::Utils::convertStrand($bed_line->{strand}));
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
        # It is a substitution
        if($ref_length == $alt_length) {
          # We shift the pos if the reference has more than one nucleotide
          $self->getEvents('snp')->addMutation($vcf_line->{chr},$vcf_line->{pos} + $ref_length - 2,$alt);
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

  # Read errors
  if(defined $self->{err_file}) {
    print STDERR "[checker] Reading error file\n" if $self->verbose;
    # Now we read the error file
    my $err_it = CracTools::Utils::getFileIterator(file =>$self->{err_file},
      parsing_method => \&CracTools::BenchCT::Utils::parseErrLine,
    );

    # TODO there can be more than one error per read...
    my $read_id = -1;
    while(my $err_line = $err_it->()) {
      if($err_line->{read_id} != $read_id) {
        $self->{err_seek_pos}[$err_line->{read_id}] = $err_line->{seek_pos};
        $read_id = $err_line->{read_id};
      }
      #$self->{err_pos}[$err_line->{read_id}] = $err_line->{pos};
      $self->{nb_errors}++;
    }
    print STDERR "[checker] ".$self->nbErrors." error(s) read\n" if $self->verbose;
  }
}

=head2 isStranded

return true is the CracTools::BenchCT::Checker used stranded simulated data, and therfore is
carreful about strand given by each tools for alignements or events.

=cut

sub isStranded {
  my $self = shift;
  return $self->{is_stranded};
}

=head2 nbReads

Return the number of reads in the simulated data

=cut

sub nbReads {
  my $self = shift;
  return $self->{nb_reads};
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
  # convert read name into a read id
  # If the read_name is already the read_id we use it directly
  my $read_id = $self->getReadId($read_name);
  my $err_lines = $self->getErrLines($read_id);
  foreach my $err (@{$err_lines}) {
    if(abs($err->{pos} - $pos) <= $CracTools::BenchCT::Const::THRESHOLD_ERR) {
      return 1;
    }
  }
  return 0;
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

  # convert read name into a read id
  # If the read_name is already the read_id we use it directly
  my $read_id = $self->getReadId($read_name);

  my $bed_line = $self->getParsedBedLine($read_id);

  if(defined $bed_line) {
    my $block_cumulated_size = 0;

    # We try to find the first block where the read is located
    foreach my $block (@{$bed_line->{blocks}}) {
      if($block->{block_end} > $pos_start) {
        # Check if the read has been mapped to the right reference
        if($block->{chr} ne $ref_name) {
          #print STDERR "Wrong chr\n";
          next;
        }
         
        # Check if protocol is stranded and the read is mapped to the wrong strand
        if($self->isStranded && $block->{strand} != $strand) {
          #print STDERR "Wrong strand\n";
          next;
        }

        # Difference between the first pos of the block to the
        # pos where the alignement is started
        my $delta = $pos_start - $block->{block_start};
       

        # Check if the position if right within a THRESHOLD_MAPPING window
        if(abs($block->{ref_start} + $delta - $ref_start) <= $CracTools::BenchCT::Const::THRESHOLD_MAPPING) {
          return $read_id + 1;
        }
      }
    }
  } else {
    warn "There is a problem with the read $read_name";
    print STDERR Dumper($bed_line);
    return 0;
  }
  return 0;
}

=head2 getBedFileHandle

Return a filehandle on the bed file.

=cut

sub getBedFileHandle {
  my $self = shift;
  my $new_fh = shift;
  if(defined $new_fh) {
    $self->{bed_fh} = $new_fh;
  }
  return $self->{bed_fh};
}

=head2 getErrFileHandle

Return a filehandle on the bed file.

=cut

sub getErrFileHandle {
  my $self = shift;
  my $new_fh = shift;
  if(defined $new_fh) {
    $self->{err_fh} = $new_fh;
  }
  return $self->{err_fh};
}

=head2 getBedLine($read_name)

=cut

sub getBedLine {
  my $self = shift;
  my $read_id = shift;

  # Get the seek position of this read in the bed file

  if(defined $read_id) {
    my $seek_pos = $self->{bed_seek_pos}[$read_id];

    # Retrieve the whole line in the bed file
    return CracTools::Utils::getLineFromSeekPos($self->getBedFileHandle,$seek_pos);
  } else {
    carp "No seek_pos for read: $read_id in the bed file";
    return undef;
  }
}

=head2 getParsedBedLine($read_name)

Given a read_name, we return the associated bed line (already parsed

=cut

sub getParsedBedLine {
  my $self = shift;
  my $read_id = shift;

  my $bed_line = $self->getBedLine($read_id);
  if(defined $bed_line) {
    # Parse this line in order to have a nice little hash
    return CracTools::BenchCT::Utils::parseGSBedLine($bed_line);
  } else {
    return undef;
  }
}

=head2 getErrLines($read_name)

=cut

sub getErrLines {
  my $self = shift;
  my $read_id = shift;

  my @err_lines;

  # Get the seek position of this read in the err file
  my $seek_pos = $self->{err_seek_pos}[$read_id];

  if(defined $seek_pos) {

    my $fh = $self->getErrFileHandle;

    # Retrieve the whole line in the bed file
    my $line = CracTools::Utils::getLineFromSeekPos($fh,$seek_pos);

    # Parse this line in order to have a nice little hash
    push(@err_lines,CracTools::BenchCT::Utils::parseErrLine($line));

    # Look at next lines, if they belong to the same read
    my $found_error = 1;
    while($found_error) {
      my $next_line = <$fh>;
      last if not $next_line;
      my $err_line = CracTools::BenchCT::Utils::parseErrLine($next_line);
      if($err_line->{read_id} == $read_id) {
        push(@err_lines,$err_line);
      } else {
        $found_error = 0;
      }
    }
  } else {
    # Should we carp?
    # This just means that the read has no error....
    #carp "No seek_pos for read: $read_name in the err file";
  }

  return \@err_lines;
}

1;
