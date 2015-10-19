use strict;
use warnings;
package CracTools::BenchCT::Utils;
# ABSTRACT: Usefull subroutines to BenchCT

=head1 PARSING SUBROUTINES

=head2 parseGSBedLineLite

=cut

sub parseGSBedLineLite {
  my $line = shift;
  my($chr,$start,$end,$name) = split("\t",$line,5);
  # Remove an eventual "chr" prefix to the reference name
  $chr =~ s/^chr//;
  return {
    chr => $chr,
    start => $start,
    end => $end,
    name => $name
  };
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
  my ($read_ids,$chr,$old_pos,$new_pos,$type,$length,$mutation) = split("\t",$line);
  my @read_ids = sort {$a <=> $b} split(":",$read_ids);
  return {
    chr       => $chr,
    old_pos   => $old_pos,
    new_pos   => $new_pos,
    type      => $type,
    length    => $length,
    mutation  => $mutation,
    read_ids   => \@read_ids,
  };
}

sub parseErrLine {
  my $line = shift;
  my %args = @_;
  my ($read_id,$pos) = split(/\s+/,$line);
  return {
    read_id => $read_id,
    pos     => $pos,
  };
}

=head2 parseChimeraLine

=cut

sub parseChimeraLine {
  my $line = shift;
  my %args = @_;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$read_ids,$nb_reads) = split("\t",$line);
  if($strand1 =~ /^[\+-]{1}$/) {
    $strand1 = CracTools::Utils::convertStrand($strand1);
    $strand2 = CracTools::Utils::convertStrand($strand2);
  }
  my @read_ids = sort {$a <=> $b} split(":",$read_ids);
  return {
    chr1      => $chr1,
    pos1      => $pos1,
    strand1   => $strand1,
    chr2      => $chr2,
    pos2      => $pos2,
    strand2   => $strand2,
    read_ids  => \@read_ids,
    nb_reads  => $nb_reads,
  };
}

1;

