#!/usr/bin/perl -w
use strict;
no strict 'refs';
use Switch;
use BitVector;

use constant THRESHOLD_POS => 5;
use constant VERBOSE => 1;

my $verbose = VERBOSE;
my $threshold_pos = THRESHOLD_POS;

sub usage {
  print STDERR "Usage: $0 
  --bed <BED file> 
  --tool-name <CRAC|BWA|Bowtie|Bowtie2|SOAP2|GSNAP>
  --false-pos <file>
  --true-pos <file>
  <result>+

  Return the sensitivity of the specified tool.
  Result file must only contain results that have been mapped once.

  Options --bed, --tool-name are mandatory.
  --false-pos records the false positive found by the tool in the 
  specified file.
  --true-pos records the true positive found by the tool in the
  specified file. \n";
  exit 1;
}


sub verbose($) {
  my ($string) = @_;
  chomp $string;
  print STDERR "$string\n" if ($verbose);
}

sub openFile($) {
  my ($filename) = @_;
  open(my $handle, $filename) or die("Unable to open $filename: $!\n");
  return $handle;
}


if (@ARGV < 5) {
  usage;
}

my ($bed, $toolName);
my ($truePosFile, $falsePosFile);
my @results;

while (@ARGV > 0) {
  switch($ARGV[0]) {
    case "--bed" {$bed = openFile($ARGV[1]); shift}
    case "--true-pos" {$truePosFile = $ARGV[1]; shift}
    case "--false-pos" {$falsePosFile = $ARGV[1]; shift}
    case "--tool-name" {$toolName = $ARGV[1]; shift}
    else {
      if ($ARGV[0] =~ /^\-/) {
        print STDERR "Unknown option ",$ARGV[0],"\n";
      } else {
        push @results, openFile($ARGV[0]);
      }
    }
  }
  shift;
}

if (! defined($bed)
  || ! @results) {
  usage;
}

my $lineNb;
my $nb_tags = 0;

# Storing data from BED file (ie. splice information)
verbose "Counting tags";
$lineNb = 0;
while (<$bed>) {
  $lineNb++;
}
seek($bed, 0, 0) or die("Unable to seek in BED file: $!\n");

$nb_tags = $lineNb;

# Tools
verbose "Processing results";

my %tags;
my $tag_num;
my $tag_name;
my $nb_match = 0;
my $nb_total_seen = 0;
my $tag_not_found;
my @cols;
my $has_matched;
my $next_tag;
my $current_line;
my $is_multiple = 0;
foreach my $tool (@results) {
  $current_line = <$tool>;
  ($tag_name) = ($current_line =~ /^(\S+)\s+/);
  my $next_line = <$tool>;
  while ($current_line) {
    if (defined $next_line) {
      ($next_tag) = ($next_line =~ /^(\S+)\s+/);
    }	

    @cols = split /\s+/, $current_line;

    if (defined $tag_name 
      && ((! defined($next_tag) || $next_tag ne $tag_name)
        && ($toolName ne 'BWA' || $current_line =~ /XT:A:U/)
        && (($toolName ne 'BWASW' && $toolName ne 'GSNAP' && $toolName ne 'Bowtie2' && $toolName ne 'TopHat2' && $toolName ne 'MapSplice') || $cols[5] ne '*')
        && ($toolName ne 'SOAP2' || $cols[3] == 1)
        && (($toolName ne 'GASSST' && $toolName ne 'TopHat2') || $current_line =~ /NH:i:1\D/))) {
      if (! $is_multiple) {

        chomp $current_line;
        my ($chr, $pos, $strand);
        if ($toolName eq 'BWA' || $toolName eq 'BWASW' || $toolName eq 'GSNAP' 
          || $toolName eq 'Bowtie2' || $toolName eq 'TopHat2' || $toolName eq 'MapSplice') {
          $nb_total_seen ++;
          $chr = $cols[2];
          $pos = $cols[3];
          $strand = ($cols[1] & 16) ? '-' : '+';
        } elsif ($toolName eq 'Bowtie') {
          $nb_total_seen ++;
          $chr = $cols[2];
          $pos = $cols[3];
          $strand = $cols[1];
        } elsif ($toolName eq 'CRAC') {
          if (($chr, $strand, $pos) = ($cols[1] =~ /^(.*?)\|(\-?1),([0-9]+)$/)) {
            $nb_total_seen ++;
            $strand = ($strand == -1) ? '-' : '+';
          } else {
            #  SAM output
            if ($current_line =~ /XU:i:1/ && $cols[2] ne '*') {
              # ($chr,$strand,$pos) = $current_line =~ /\sXO:Z:(\S+)\|(-?1),([0-9]+)\s/;
              $chr = $cols[2];
              $pos = $cols[3];
              $strand = ($cols[1] & 16) ? '-' : '+';
              # $strand = ($strand == -1) ? '-' : '+';
              $nb_total_seen++;
            }
          }
        } elsif ($toolName eq 'SOAP2') {
          $nb_total_seen++;
          $chr = $cols[7];
          $strand = $cols[6];
          $pos = $cols[8];
        } elsif ($toolName eq 'GASSST') {
          # NH:i:1 -> unique location
          # NH:i:k -> k>1 multiple locations
          # NM:i:k -> k differences (edit distance)
          $nb_total_seen ++;
          $pos = $cols[3];
          $strand = ($cols[1] & 16) ? '-' : '+';
          $chr = $cols[2];
        } else {
          print STDERR "I don't know the tool $toolName. \nMaybe you mistyped it or you should modify my code so that I can manage it!\n";
          exit 11;
        }

        if (defined $chr) {
          if ($chr !~ /^chr/) {
            $chr = 'chr'.$chr;
          }
          $tags{chr}{$tag_name} = $chr;
          $tags{strand}{$tag_name} = $strand;
          $tags{pos}{$tag_name} = $pos;
        }
      } else {
        $is_multiple = 0;
      }
    } elsif ((defined($next_tag) && $next_tag eq $tag_name)
      || (($toolName eq 'GASSST' || $toolName eq 'TopHat2') && $current_line !~ /NH:i:1\D/)) {
      $is_multiple = 1;
    }

    $current_line = $next_line;
    $next_line = <$tool>;
    $tag_name = $next_tag;
  }
  close($tool);
}



my $truePosHandle; 
my $falsePosHandle;

if (defined($truePosFile)) {
  open($truePosHandle, ">".$truePosFile) or die("Unable to open file $truePosFile: $!\n");
}
if (defined($falsePosFile)) {
  open($falsePosHandle, ">".$falsePosFile) or die("Unable to open file $falsePosFile: $!\n");
}

verbose "Reading BED file";

$lineNb = 0;
my $tagID;
my @positions;
while (<$bed>) { 
  my @cols = split /\s+/;

  if (defined($tags{chr}{$lineNb})) {
    $tagID = $lineNb;
  } elsif(defined($tags{chr}{$cols[3]})) {
    $tagID = $cols[3];
  } else {
    $tagID = undef;
  }
  @positions = ();

  if (defined($tagID)) {
    my @fragPos = split /,/, $cols[10];
    my @fragLength = split /,/, $cols[11];
    my $posStart = $cols[1];
    my $readLength = 0;
    for (my $i = 0; $i < $cols[9]; $i++) {
      # Is it a chimera?
      if ($fragLength[$i] =~ /@/) {
        my ($chr, $strand, $pos) = ($fragLength[$i] =~ /^(.*?)@([+-]?)(\d+)$/);
        push @positions, {chr => $chr, 
          pos => $pos,
          strand => $strand,
          length => $fragPos[$i]};
      } else {
        push @positions, {chr => $cols[0], 
          pos => $posStart + $fragLength[$i],
          strand => $cols[5],
          length => $fragPos[$i]};
      }
      $readLength += $fragPos[$i];
    }
    @fragPos = ();
    @fragLength = ();

    # Test if the location is correct.
    my $found = 0;
    my $readLengthProcessed = 0;


    # NOTE If one chunk of the read is found, then it is consider as "matched"
    foreach my $position (@positions) {
      my $multiplier = ($position->{strand} eq '+') ? 1 : -1;
      # ``Normal'' match
      if ($tags{chr}{$tagID} eq $position->{chr}
        && $tags{strand}{$tagID} eq $position->{strand}
        && abs($tags{pos}{$tagID} - $position->{pos})
        <= ($position->{length} + $threshold_pos)) {
        $found = 1;
      }
      $readLengthProcessed += $position->{length};
      if ($found) {
        $nb_match++;

        if (defined($truePosFile)) {
          print $truePosHandle "Flux=",$position->{pos},'@',
          $position->{strand},$position->{chr};
          print $truePosHandle "\t",$toolName,
          '=',$tags{pos}{$tagID},'@',$tags{strand}{$tagID},
          $tags{chr}{$tagID},"\t$tagID\n";
        }

        last;
      }
    }	

    # Writing not found
    if (! $found && (defined($falsePosFile))) {
      my $nb_printed = 0;
      my $chimera_print;
      print $falsePosHandle "Flux=";
      foreach my $position (@positions) {
        my $multiplier = ($position->{strand} eq '+') ? 1 : -1;
        if ($nb_printed > 0) {
          print $falsePosHandle "|";
        }
        print $falsePosHandle $position->{pos},'@',
        $position->{strand},$position->{chr};
        $nb_printed++;
      }
      print $falsePosHandle "\t",$toolName,
      '=',$tags{pos}{$tagID},'@',$tags{strand}{$tagID},
      $tags{chr}{$tagID},"\t$tagID\n";
    }
  }

  delete @positions[0..$#positions];
  $lineNb++;
}
close($bed);

printf "nb_tags: %d\nnb_uniq_matched: %d (%.2f%%)\nnb_true: %d (%.2f%%) (%.2f%%)\n",
$nb_tags, $nb_total_seen, $nb_total_seen*100./$nb_tags, $nb_match, 
$nb_match*100./$nb_tags, $nb_match*100./$nb_total_seen;
