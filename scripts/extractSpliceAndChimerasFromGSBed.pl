#! /usr/bin/perl
#
use strict;
use warnings;

use CracTools::Utils;

my %splices;
my %chimeras;
my $is_stranded = 0;

if(@ARGV < 3) {
  print STDERR "Usage: extractSplicesAndChimerasFromGSBed.pl file.bed file.splice file.chimeras\n";
  exit 0;
}

my $bed_file = shift;
my $splice_file = shift;
my $chimera_file = shift;

my $bed_it = CracTools::Utils::getFileIterator(file =>$bed_file,
  parsing_method => \&parseGSBedLine,
);

my $id = 0;
while(my $bed_line = $bed_it->()) {
  # Loop over blocks to find splices and chimeras
  for(my $i=1; $i < @{$bed_line->{blocks}}; $i++) {

    # Check that the splice is not chimeric
    if ($bed_line->{blocks}[$i-1]->{chr} ne $bed_line->{blocks}[$i]->{chr} ||
      $bed_line->{blocks}[$i-1]->{strand} ne $bed_line->{blocks}[$i]->{strand} ||
      $bed_line->{blocks}[$i-1]->{ref_end} >= $bed_line->{blocks}[$i]->{ref_start}) {

      #print STDERR "read: ".$bed_line->{name}."\n";
    #next if $bed_line->{blocks}[$i]->{strand} eq '+' && $bed_line->{blocks}[$i-1]->{ref_end} > $bed_line->{blocks}[$i]->{ref_start};
    #next if $bed_line->{blocks}[$i]->{strand} eq '-' && $bed_line->{blocks}[$i-1]->{ref_end} < $bed_line->{blocks}[$i]->{ref_start};
      my $chr1 = $bed_line->{blocks}[$i-1]->{chr};
      my $chr2 = $bed_line->{blocks}[$i]->{chr};
      my $strand1 = $bed_line->{blocks}[$i-1]->{strand};
      my $strand2 = $bed_line->{blocks}[$i]->{strand};
      my $pos1 = $bed_line->{blocks}[$i-1]->{ref_end};
      my $pos2 = $bed_line->{blocks}[$i]->{ref_start};
      my $key = $chr1.'@'.$strand1.'@'.$pos1.'@'.$chr2.'@'.$strand2.'@'.$pos2;
      my $reversed_key = $chr1.'@'.($strand1*-1).'@'.$pos1.'@'.$chr2.'@'.($strand2*-1).'@'.$pos2;
      #my $reversed_key = $chr2.'@'.($strand2*-1).'@'.$pos2.'@'.$chr1.'@'.($strand1*-1).'@'.$pos1;

      if(!$is_stranded && defined $chimeras{$reversed_key}) {
        $key = $reversed_key;
      }

      if(defined $chimeras{$key}) {
        push(@{$chimeras{$key}},$id);
      } else {
        $chimeras{$key} = [$id];
      }

    } else {
      #
      # Extract splice coordinates
      my $chr = $bed_line->{blocks}[$i]->{chr};
      my $start = $bed_line->{blocks}[$i-1]->{ref_end};
      my $end = $bed_line->{blocks}[$i]->{ref_start};
      my $length = $end - $start;
      my $strand = $is_stranded? $bed_line->{blocks}[$i]->{strand} : 1;
      my $key = $chr.'@'.$start.'@'.$length.'@'.$strand;

      if(defined $splices{$key}) {
        push(@{$splices{$key}},$id);
      } else {
        $splices{$key} = [$id];
      }
    }
  }
  $id++;
}

my $splice_fh = CracTools::Utils::getWritingFileHandle($splice_file);

foreach my $splice_key (keys %splices) {
  my ($chr,$start,$length,$strand) = split('@',$splice_key);
  print $splice_fh join("\t",($chr,$start,$start+$length,join(":",@{$splices{$splice_key}}),scalar @{$splices{$splice_key}},CracTools::Utils::convertStrand($strand))),"\n";
}

my $chim_fh = CracTools::Utils::getWritingFileHandle($chimera_file);

foreach my $chimera_key (keys %chimeras) {
  my ($chr1,$strand1,$pos1,$chr2,$strand2,$pos2) = split('@',$chimera_key);
  print $chim_fh join("\t",($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,join(":",@{$chimeras{$chimera_key}}),scalar @{$chimeras{$chimera_key}})),"\n";
}

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
