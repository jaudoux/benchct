use strict;
use warnings;
package CracTools::BenchCT::Checker;
# ABSTRACT: Benchmarking CRAC results

use CracTools::Utils;
use CracTools::BenchCT::Const;
use CracTools::BenchCT::Events::Mutation;
use CracTools::BenchCT::Events::Splice;
use Data::Dumper; #for debugging

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $is_stranded = $args{is_stranded};
  $is_stranded = 0 if !defined $is_stranded;

  my $self = bless {
    is_stranded => $is_stranded,
    info_file  => $args{info_file},
    bed_file   => $args{bed_file},
    bed_fh    => CracTools::Utils::getReadingFileHandle($args{bed_file}),
    error_file => $args{error_file},
    read_ids   => {}, # TODO avoid this conversion hash and use a interval tree to seek for bed line that could correspond to a read alignement
    bed_seek_pos => [],
    nb_reads   => 0,
    events => { snp         => CracTools::BenchCT::Events::Mutation->new( threshold => $CracTools::BenchCT::Const::THRESHOLD_SNP ), # Create a special CracTools::BenchCT::Events type for snp that would use a genomeMask
                insertion   => CracTools::BenchCT::Events::Mutation->new( threshold => $CracTools::BenchCT::Const::THRESHOLD_INS ),
                deleletion  => CracTools::BenchCT::Events::Mutation->new( threshold => $CracTools::BenchCT::Const::THRESHOLD_DEL ),
                splice      => CracTools::BenchCT::Events::Splice->new( threshold => $CracTools::BenchCT::Const::THRESHOLD_SPLICE ),
              },
  }, $class;

  $self->_init();

  return $self;
}

sub _init {
  my $self = shift;
  
  # First we read the bed files and we record seek positions of each read
  my $bed_it = CracTools::Utils::getFileIterator(file =>$self->{bed_file},
    parsing_method => \&parseGSBedLine,
  );

  my $id = 0;
  while(my $bed_line = $bed_it->()) {
    # Set number to read ids
    $self->{read_ids}->{$bed_line->{name}} = $id;
    # Set seek pos of each bed alignement according to the read id
    $self->{bed_seek_pos}[$id] = $bed_line->{seek_pos};
    # Loop over blocks to find splices and chimeras
    for(my $i=1; $i < @{$bed_line->{blocks}}; $i++) {

      # Check that the splice is not chimeric
      next if $bed_line->{blocks}[$i-1]->{chr} ne $bed_line->{blocks}[$i]->{chr};
      next if $bed_line->{blocks}[$i-1]->{strand} ne $bed_line->{blocks}[$i]->{strand};
      next if $bed_line->{blocks}[$i]->{strand} eq '+' && $bed_line->{blocks}[$i-1]->{ref_end} > $bed_line->{blocks}[$i]->{ref_start};
      next if $bed_line->{blocks}[$i]->{strand} eq '-' && $bed_line->{blocks}[$i-1]->{ref_end} < $bed_line->{blocks}[$i]->{ref_start};

      # Extract splice coordinates
      my $chr = $bed_line->{blocks}[$i]->{chr};
      my $start = $bed_line->{blocks}[$i-1]->{ref_end};
      my $end = $bed_line->{blocks}[$i]->{ref_start};
      my $length = $end - $start;
      my $strand = $bed_line->{blocks}[$i]->{strand};

      $self->getEvents('splice')->addSplice($chr,$start,$length,$strand);
    }
    $id++;
  }

  $self->{nb_reads} = $id;

  # Secondly we read the info file to register mutations
  my $info_it = CracTools::Utils::getFileIterator(file =>$self->{info_file},
    parsing_method => \&parseInfoLine, # we use our own parsing method
    skip => 1, # skip the first line
  );

  # store tag ids for each mutation in an Interval Tree (one for each type of mutations
  while(my $info_line = $info_it->()) {
    my $mutation_type = $info_line->{type};
    $mutation_type = 'snp' if $info_line->{type} eq 'sub'; # rename substition in snps
    if(defined $self->getEvents($mutation_type)) {
      $self->getEvents($mutation_type)->addMutation($info_line);
    }
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
  return $self->{read_ids}->{$read_name};
}

=head2 getEvents($type)

Return the CracTools::BenchCT::Events of a specific type (snp, ins, del, splice, chimera)

=cut

sub getEvents {
  my $self = shift;
  my $type = shift;
  return $self->{events}->{$type};
}

=head2 isTrueMutation($type,$chr,$pos,$length)

Return true is there is a mutation at this position

=cut

sub isTrueMutation {
  my $self = shift;
  my ($type,$chr,$pos,$length) = @_;
  my $snp_events = $self->getEvents($type);
  return $snp_events->isTrueMutation($chr,$pos,$length);
}

=head2 isTrueSplice($chr,$start,$length,$strand)

Return true is there is a mutation at this position

=cut

sub isTrueSplice {
  my $self = shift;
  my ($chr,$start,$length,$strand) = @_;
  my $splice_events = $self->getEvents('splice');
  return $splice_events->isTrueSplice($chr,$start,$length,$strand);
}

=head2 getBedFileHandle

Return a filehandle on the bed file.

=cut

sub getBedFileHandle {
  my $self = shift;
  return $self->{bed_fh};
}

=head2 getBedLine($read_name)

Given a read_name, we return the associated bed line (already parsed

=cut

sub getBedLine {
  my $self = shift;
  my $read_name = shift;

  # concert read name into a read id
  # If the read_name is already the read_id we use it directly
  my $read_id = $read_name =~ /^\d+$/? $read_name : $self->getReadId($read_name);

  # Get the seek position of this read in the bed file
  my $seek_pos = $self->{bed_seek_pos}[$read_id];

  # Retrieve the whole line in the bed file
  my $line = CracTools::Utils::getLineFromSeekPos($self->getBedFileHandle,$seek_pos);

  # Parse this line in order to have a nice little hash
  my $bed_parsed = parseGSBedLine($line);

  return $bed_parsed;
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
  my $bed_line = $self->getBedLine($read_name);
  my $first_mapped_block;
  my $block_cumulated_size = 0;

  # We try to find the first block where the read is located
  foreach my $block (@{$bed_line->{blocks}}) {
    if($block->{block_end} > $pos_start) {
      $first_mapped_block = $block;
      last;
    }
  }

  # There should be a corresponding block, otherwise there is a problem
  # between the BED file and the query
  if(defined $first_mapped_block) {

    # Check if the read has been mapped to the right reference
    if($first_mapped_block->{chr} ne $ref_name) {
      #print STDERR "Wrong chr\n";
      return 0;
    }
     
    # Check if protocol is stranded and the read is mapped to the wrong strand
    if($self->isStranded && CracTools::Utils::convertStrand($first_mapped_block->{strand}) != $strand) {
      #print STDERR "Wrong strand\n";
      return 0;
    }

    # Difference between the first pos of the block to the
    # pos where the alignement is started
    my $delta = $pos_start - $first_mapped_block->{block_start};
   

    # Check if the position if right within a THRESHOLD_MAPPING window
    if(abs($first_mapped_block->{ref_start} + $delta - $ref_start) <= $CracTools::BenchCT::Const::THRESHOLD_MAPPING &&
      $first_mapped_block) {
      return 1;
    } else {
      #print STDERR "Bad position : ".($first_mapped_block->{ref_start} + $delta)."\n";
      #print STDERR Dumper($bed_line);
      return 0;
    }
  } else {
    warn "There is a problem with the read $ref_name";
    return 0;
  }
}

=head2 parseGSBedLine

This is a special parsing method for Genome Simulator bed that can contain chimeras

=cut

sub parseGSBedLine {
  my $line = shift;
  my %args = @_;
  my($chr,$start,$end,$name,$score,$strand,$thick_start,$thick_end,$rgb,$block_count,$block_size,$block_starts) = split("\t",$line);
  
  # Remove an eventual "chr" prefix to the reference name
  $chr =~ s/^chr//;

  # Manage chimeric end pos
  my ($end_chr,$end_strand,$end_pos) = $end =~ /(\S+)@([+-])(\d+)/;
  $end_chr =~ s/^chr// if defined $end_chr;
  $end_chr = $chr if !defined $end_chr;
  $end_strand = $start if !defined $end_strand;

  # Manage blocks
  my @blocks;
  my @block_size = split(",",$block_size);
  my @block_starts = split(",",$block_starts);
  my $cumulated_block_size = 0;
  for(my $i = 0; $i < $block_count; $i++) {
    # manage chimeric blocks
    my ($block_chr,$block_strand,$block_start) = $block_starts[$i] =~ /(\S+)@([+-])(\d+)/;
    $block_starts[$i] = $block_start if defined $block_start;
    my $ref_start = defined $block_start? $block_start : $block_starts[$i] + $start;
    my $ref_end = defined $block_start? $block_start + $block_size[$i] : $block_starts[$i] + $start + $block_size[$i];
    $block_chr =~ s/^chr// if defined $block_chr;
    $block_chr = $chr if !defined $block_chr;
    $block_strand = $strand if !defined $block_strand;
    push(@blocks,{size        => $block_size[$i], 
                 start        => $block_starts[$i], 
                 end          => $block_starts[$i] + $block_size[$i],
                 block_start  => $cumulated_block_size,
                 block_end    => $cumulated_block_size + $block_size[$i],
                 ref_start    => $ref_start,
                 ref_end      => $ref_end,
                 chr          => $block_chr,
                 strand       => CracTools::Utils::convertStrand($block_strand), # Convert strand from '+/-' to '1/-1' format
               });
    $cumulated_block_size += $block_size[$i];
  }

  return { chr        => $chr,
    start       => $start, 
    end         => $end, 
    end_chr     => $end_chr,
    end_pos     => $end_pos,
    name        => $name,
    score       => $score, 
    strand      => CracTools::Utils::convertStrand($strand), # Convert strand from '+/-' to '1/-1' format
    thick_start => $thick_start,
    thick_end   => $thick_end,
    rgb         => $rgb,
    blocks      => \@blocks,
  };
}

=head2 parseInfoFile

This is a simple parsing method for an info line

=cut

sub parseInfoLine {
  my $line = shift;
  my %args = @_;
  my ($tag_ids,$chr,$old_pos,$new_pos,$type,$length,$mutation) = split("\t",$line);
  my @tag_ids = sort {$a <=> $b} split(":",$tag_ids);
  return {
    chr       => $chr,
    old_pos   => $old_pos,
    new_pos   => $new_pos,
    type      => $type,
    length    => $length,
    mutation  => $mutation,
    tag_ids   => \@tag_ids,
  };
}

1;
