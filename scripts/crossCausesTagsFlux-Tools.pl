#!/usr/bin/perl -w
use strict;
no strict 'refs';
use Switch;
use BitVector;

use constant THRESHOLD_SNP => 5;
use constant THRESHOLD_INS => 5;
use constant THRESHOLD_DEL => 5;
use constant THRESHOLD_CHIMERA => 20;
use constant THRESHOLD_ERR => 5;
use constant THRESHOLD_SPLICE => 5;

use constant VERBOSE => 1;

# Append to key in infosPosChr to know if
# the cause has already been seen
use constant KEY_CHECK => 'check';

my $verbose = VERBOSE;
my $debugEntry = 0;
my $stranded = 0;
my %threshold;
$threshold{snp} = THRESHOLD_SNP;
$threshold{ins} = THRESHOLD_INS;
$threshold{del} = THRESHOLD_DEL;
$threshold{chimera} = THRESHOLD_CHIMERA;
$threshold{errors} = THRESHOLD_ERR;
$threshold{splice} = THRESHOLD_SPLICE;

sub usage {
  print STDERR "Usage: $0 
  --bed <BED file> --flux-err <err file> --mutations <info file>
  --chr-lengths <.conf file>
  --snp <CRAC snp> --chimera <CRAC chimera> 
  --indel <CRAC (tag|genome) indel> 
  --errors <CRAC err> --splicing <CRAC overlap|splice> 
  --single 
  (--tools-[multiple|single|all] <Tool's output file>)+
  --tool-name <CRAC|BWA|Bowtie|SOAP2|TopHat|TopHat2|MapSplice|TopHatFusion|GSNAP|SensitivityFlux>

  Counts the number of causes found, missing by analising Flux files
    and CRAC (or another tool)'s files.
  Percentages are output on STDIN.

  --single: takes into account only breaks that are labelled as single.
  That allows to compute accuracy without false positives
  due to duplications.

  You can have an additional output by specifying filenames
  --output-ok: saves the causes found correctly
  --output-missing: saves the causes not found
  --output-false: save the false positives
  --output-tag-missing: output the tag numbers for which a cause is missing.
  --output-bad-classification: output the tag numbers whose classification is 
  not correct.
  --output-expression: output the expression of the true postiive junctions
  (chimeric or splicing)

  You can change the threshold for each cause (ie. 
  tolerance between the actual position and the position found):
  --threshold-snp: threshold for SNPs (default ".THRESHOLD_SNP.")
  --threshold-ins: for insertions     (default ".THRESHOLD_INS.")
  --threshold-del: for deletions      (default ".THRESHOLD_DEL.")
  --threshold-chimera: for chimeras   (default ".THRESHOLD_CHIMERA.")
  --threshold-errors: for errors      (default ".THRESHOLD_ERR.")
  --threshold-splice: for splices     (default ".THRESHOLD_SPLICE.")
  --threshold: change all the thresholds to this value
  ";
  exit 1;
}

sub verbose($) {
  my ($string) = @_;
  chomp $string;
  print STDERR "$string\n" if ($verbose);
}

sub verbose_progress($) {
  my ($nb) = @_;
#     print STDERR "\r$nb " 
# 	if ($nb % 1000 == 0);
}

sub openFile($) {
  my ($filename) = @_;
  open(my $handle, $filename) or die("Unable to open $filename: $!\n");
  return $handle;
}

sub writeFile($) {
  my ($filename) = @_;
  open(my $handle, ">".$filename) or die("Unable to open $filename: $!\n");
  return $handle;
}

sub read_file_line {
  my ($fh) = shift;

  if ($fh and my $line = <$fh>) {
    return $line;
  }
  return;
}

# Increment statistics for a given tag
# $tag: a tag number
# $type: a mutation type
# $pos: an array, one position per tag that tells us whether the tag is located
#       once or more.
# $stats: ref to a hash with statistics on the localisations.
sub incStatTag {
  my ($tag, $type, $pos, $stats, $output, $inc) = @_;
  $inc = 1
  if (! defined $inc);
  $stats->{$type}->{total}+=$inc;
  if (defined $pos->[$tag]) {
    print {$output->{ok}} "$type $tag\n"
    if (defined $output->{ok});
    if ($pos->[$tag] == 1) {
      $stats->{$type}->{unique}+=$inc;
    } else {
      $stats->{$type}->{multiple}+=$inc;
    }
  } else {
    print {$output->{missing}} "$type $tag\n"
    if (defined($output->{missing}));
  }
}	    

sub count_false_positives($$$$$) {
  my ($allMasks, $stat, $type, $output, $infosPos) = @_;
  my $mask = $allMasks->{$type};
  my $falsePos = 0;
  my $sumSet = 0;
  if (defined $mask) {
    foreach my $chr (keys %$mask) {
      if ($mask->{$chr}->nb_set > 0) {
        $sumSet += $mask->{$chr}->nb_set;
        my $handle = $output->{false};
        my $pos = 0;
        while (($pos = $mask->{$chr}->succ($pos)) != -1) {
          if (defined $infosPos) {
            foreach my $key (keys %{$infosPos->{$pos.'@'.$chr}}) {
              next if ($key eq 'notSeenYet');
              my $value = $infosPos->{$pos.'@'.$chr}->{$key};
              if (! defined $infosPos->{KEY_CHECK}->{$pos.'@'.$chr}->{$key}) {
                $falsePos++;
                print_cause($type, $output->{false}, $chr, $pos,
                  ($value =~ /^0$/) ? $key : $value.','.$key);
              }
            }
          } else {
            $falsePos++;
            if (defined $output->{false}){
              print $handle "$type $chr,$pos\n";
            }
          }
          $pos++;
        }
      }
    }
    $stat->{falsePositives} = $falsePos;
  }
}

sub print_cause($$$$$) {
  my ($type, $output, $chr, $pos, $key) = @_;

  if (defined $output) {
    print $output "$type $chr,$pos";
    print $output " $key"
    if (defined $key);
    print $output "\n";
  }
}

# Output the read numbers that have not been correctly classified
sub output_bad_classified($$$) {
  my ($reads, $current, $output) = @_;
  my $pos = 0;
  while (($pos = $reads->{$current}->succ($pos)) != -1) {
    print {$output->{'bad-classification'}} "$current $pos\n";
    $pos++;
  }
}

# Count the number of different causes found
sub count_causes($) {
  my ($mask) = @_;
  my $causes = 0;
  verbose "Counting causes";
  foreach my $chr (keys %{$mask}) {
    $causes += $mask->{$chr}->nb_set;
  }
  return $causes;
}

sub giveExonicPosition($$$$) {
  # In fact totalLength seems to be useless..
  my ($relStart, $totalLength, $currentFragLength, $shift) = @_;
  my $result;
  if ($shift =~ /@/) {
    my $temp;
    ($result->{chr}, $temp) = split /@/, $shift;
    ($result->{strand}, $result->{pos}) = split /[+-]/, $temp;
    $result->{pos} += $currentFragLength;
    $result->{chimera}=1;
  } else {
    $result->{strand} = $relStart->{strand};
    $result->{chr} = $relStart->{chr};
    $result->{pos} = $relStart->{pos} + $currentFragLength + $shift;
  }
  $result->{chr} =~ s/^chr//;
  return $result;
}

# Initialises the given bit vector for the specified type.
# $mask: ref of hash of bit vectors
# $type: a reference to an array containing hash keys.
# $chrLengths: handle of the file to be read for retrieving the chr lengths.
sub init_bit_vectors($$$) {
  my ($mask, $type, $chrLengths) = @_;
  verbose "Initialising bit vectors";
# Reading Conf file, and initialise bit vectors.
#Skip the first line (nb. of chr)
  my $junk = <$chrLengths>;
  while (<$chrLengths>) {
    my $length = <$chrLengths>;
    chomp;
    chomp $length;
    foreach my $elem (@$type) {
      $mask->{$elem}{$_} = BitVector->new($length);
    }
  }
  seek $chrLengths, 0, 0;
}

# @param ref to the src hash containing bit vectors
# @return a reference to a copy of the hash containing copies of bit vectors
sub copy_bit_vectors($) {
  my ($src) = @_;
  my $dest = {};
  foreach my $key (keys %$src) {
    $dest->{$key} = $src->{$key}->copy;
  }
  return $dest;
}

# Try to find at pos $pos on chr $chr a mutation of type $type.
# $infosPosChr: is a pointer on an array that may contain additional information
#               (position and chr.) to be matched for a given type to be found.
# $resPosChr: it is what $infosPosChr must match to.
# $mask records the mutations found.
# $causesFound: a hash (by type) of hash of all the causes already found
# $notInFP: a hash similar to causesFound that telle what causes are not FP.
# $countTags: records the number of tags which have seen the mutation.
# $threshold defines the tolerance for finding the mutation in terms of
#   number of base pairs.
# $stats defines statistics on the number of elements correctly found, not found, ...
# $output contains some handle to output informations about what we found and what we 
#   didn't.
# return 1, if there are still things to be found at that position.

sub findType($$$$$$$$$$$$) {
  my ($type, $chr, $pos, $infosPosChr, $resPosChr, $mask, $causesFound, $notInFP, $countTags, $threshold, $stats, $output) = @_;

  $stats->{$type}->{total}++;

  my $additionalChr = '';
  my @posBV;
  my $current_threshold = $threshold->{$type};
  my $posBV;
  my $posBVFound;

  if ($debugEntry) {
    print "finding type...\n";
  }

  my $tmpPos = $pos;
  # Check all the possibilities because of alternative splicing
  # and because we may allow a threshold
  while ($current_threshold >= 0 
    && ($posBV = $mask->{$type}{$chr}->succ($tmpPos, $current_threshold)) != -1) {
    $current_threshold -= $posBV - $tmpPos + 1;
    $tmpPos = $posBV + 1;
    push @posBV, $posBV;
  }
  $current_threshold = $threshold->{$type}-1;
  $tmpPos = $pos-1;
  while ($current_threshold >= 0 
    && ($posBV = $mask->{$type}{$chr}->prev($tmpPos, $current_threshold)) != -1) {
    $current_threshold -= $tmpPos - $posBV + 1;
    $tmpPos = $posBV - 1;
    push @posBV, $posBV;
  }

  if ($debugEntry) {
    print "We are here! ",scalar(@posBV),"\n";
    foreach $posBV (@posBV) {
      print "\t$posBV\n";
    }
    print "pos = $pos, chr=$chr, otherPos = ".((defined $resPosChr->[0])?$resPosChr->[0]:"")."\n";
  }

  my @found;

  foreach $posBV (@posBV) {
    if ($debugEntry) {
      print "posBV = $posBV\n";
    }
    # If the additional informations are defined, we also test them
    # to see if we match the corresponding stuff.
    if (defined($infosPosChr)) {
      if (defined($infosPosChr->{$posBV.'@'.$chr})) {
        foreach my $elem (keys %{$infosPosChr->{$posBV.'@'.$chr}}) {
          next if ($elem eq 'notSeenYet');
          if ($debugEntry) {
            print "elem = $elem  ($resPosChr->[0])\n";
          }
          if (abs($elem - $resPosChr->[0]) <= $threshold->{$type}) {
            if (! defined($resPosChr->[1])
              || $infosPosChr->{$posBV.'@'.$chr}->{$elem} eq $resPosChr->[1]) {
              if (defined($resPosChr->[1])) {
                $additionalChr = $resPosChr->[1].',';
              }
              if ($debugEntry) {
                print "\tfound \n";
              }
              push @found, {'otherPos' => $elem, 'bv' => $posBV};
            }
          }
        }
      }
    } else {
      # Everything but splice and chimera
      push @found, {'bv' => $posBV}
    }
  }

  # Sort the entries by 'quality' the closest to what is required are put first.
  @found = sort {abs($a->{bv} - $pos) <=> abs($b->{bv} - $pos) 
  || ! defined $a->{otherPos} 
  || abs($a->{otherPos} - $resPosChr->[0]) <=> abs($b->{otherPos} - $resPosChr->[0])
  } @found;

  if (@found) {
    my $nb_found = 0;
    foreach my $occurrence (@found) {
      my $additionalPos = '';
      my $additionalKey = '';
      $posBV = $occurrence->{bv};

      if (defined($occurrence->{otherPos})) {
        $additionalKey = '-'.$occurrence->{otherPos};
        $additionalPos = $occurrence->{otherPos};
      }

      if (defined $countTags && defined $countTags->{$type}{$chr}{$posBV.$additionalKey}) {
        # We account for all the reads sharing the same mutation
        $stats->{$type}->{tagsOK} += $countTags->{$type}{$chr}{$posBV.$additionalKey};
        # And we delete the entry so that we don't count it many times
        delete $countTags->{$type}{$chr}{$posBV.$additionalKey};
      }
      print_cause($type, $output->{ok}, $chr, $posBV, $additionalChr.$additionalPos
        .(defined($resPosChr->[2]) ? "\t".$resPosChr->[2]:""));

      # Just make sure that we didn't see that cause so that we don't count
      # it twice.
      if (! defined $causesFound->{$type}->{$chr}->{$posBV.$additionalKey}) {
        $stats->{$type}->{ok}++;
        $causesFound->{$type}->{$chr}->{$posBV.$additionalKey}=1;
      }
      $nb_found++;

      # Unset it so that we don't count this one
      # in the false postives!
      if (defined $infosPosChr->{$posBV.'@'.$chr}->{$additionalPos}) {
        if (!defined($infosPosChr->{KEY_CHECK}->{$posBV.'@'.$chr}->{$additionalPos})) {
          $infosPosChr->{KEY_CHECK}->{$posBV.'@'.$chr}->{$additionalPos} = -1;
          $infosPosChr->{$posBV.'@'.$chr}->{notSeenYet}--;

          if ($infosPosChr->{$posBV.'@'.$chr}->{notSeenYet}==0) {
            $notInFP->{$type}->{$chr}->{$posBV}=1;
          }
        }
      } else {
        $notInFP->{$type}->{$chr}->{$posBV}=1;
      }
    }
    return scalar @found;
  } else {
    $stats->{$type}->{notFound}++;
    my $additionalInfos = '';
    if (defined ($resPosChr->[0])) {
      $additionalInfos = $resPosChr->[0];
      if (defined ($resPosChr->[1])) {
        $additionalInfos .= '@'.$resPosChr->[1];
      }
    }
    print_cause($type, $output->{missing}, $chr, $pos, $additionalInfos);
    return 0;
  }
}

sub delete_in_mask($$) {
  my ($found, $mask) = @_;

  foreach my $type (keys %$found) {
    foreach my $chr (keys %{$found->{$type}}) {
      foreach my $pos (keys %{$found->{$type}->{$chr}}) {
        $mask->{$type}{$chr}->unset($pos);
      }
    }
  }
}


# Processing indels and SNPs.
# $typesAllowed: hash ref whose keys are allowed types (types that can be searched)
# $mutations: handle to the INFO file
# $mask : ref to the hash of bit vectors
# $countTags: ref to the hash of tag count per chromosomic position.
# $reads : ref to the hash that contains a bit vector for each cause
#          identifying which read contains a given cause
# $threshold: ref to the hash of threshold
# $stats: ref to the hash of statistics
# $output: ref to the hash of handles for the output.
sub read_INFO_file($$$$$$$$) {
  my ($typesAllowed, $mutations, $mask, $countTags, $reads, $threshold, $stats, $output) = @_;

  foreach my $current (keys %$typesAllowed) {
    $stats->{$current}->{count} = count_causes(\%{$mask->{$current}});
  }


  my %causesFound;
  my %notInFP;
  verbose "Reading INFO file";
# Storing data from INFO file (ie. mutations)
  while (<$mutations>) {
    my @line = split /\s+/;
    my @num_tags = split /:/, $line[0];
    my @type = split /,/, $line[4];
    my $type = $type[0];
    if ($line[5] =~ /^\d+$/) {
      if ($line[5] == 1) {
        $type = 'snp';
      }
      # We do not treat chimeras here since the information in .INFO is not sufficient.
      # We need to know the position of the first exon in the chimera which is not given
      # by the .info.
      if (defined $typesAllowed->{$type}) {
        my ($chr, $pos) = ($line[1], $line[2]);

        if (! findType($type, $chr, $pos, undef, undef, 
            $mask, \%causesFound, \%notInFP, $countTags, $threshold, $stats, $output)
          && defined $output->{"tag-missing"}) {
          print {$output->{"tag-missing"}} "$type\t$line[0]\n";
        } 
        # We unset the bits in the BV corresponding to the given cause
        foreach my $read_nb (@num_tags) {
          $reads->{$type}->unset($read_nb);
        }
      }
    }
  }
  seek $mutations, 0, 0;

  delete_in_mask(\%notInFP, $mask);

  verbose "Counting false positives";
  foreach my $current (keys %$typesAllowed) {
    if (defined $mask->{$current}) {
      count_false_positives($mask, $stats->{$current}, $current
        , $output, undef);
    }
    $stats->{$current}->{tagsBadClassified} = $reads->{$current}->nb_set;
    if (defined($output->{'bad-classification'})) {
      output_bad_classified($reads, $current, $output);
    }
    delete $mask->{$current};
    delete $reads->{$current};
  }
}

if (@ARGV < 4) {
  usage;
}

my ($bed, $fluxErr, $mutations);
my ($snp, $chimera, @indel, $none, %tools, $toolName,
  $errors, @splicing, $chrLengths, $single);
my %output;

while (@ARGV > 0) {
  switch($ARGV[0]) {
    case "--bed" {$bed = openFile($ARGV[1]); shift}
    case "--flux-err" {$fluxErr = openFile($ARGV[1]); shift}
    case "--mutations" {$mutations = openFile($ARGV[1]); shift}
    case "--snp" {$snp = openFile($ARGV[1]); shift}
    case "--chimera" {$chimera = openFile($ARGV[1]); shift}
    case "--chr-lengths" {$chrLengths = openFile($ARGV[1]); shift}
    case "--indel" {push @indel, openFile($ARGV[1]); shift}
    case "--errors" {$errors = openFile($ARGV[1]); shift}
    case "--splicing" {push @splicing, openFile($ARGV[1]); shift}
    case "--single" {$single=1;}
    case (/^\-\-tools/) {
      my ($type) = ($ARGV[0] =~ /^\-\-tools-([a-z]+)$/);
      if (! defined $type || ($type ne 'multiple'
          && $type ne 'single'
          && $type ne 'all')) {
        print STDERR "Please give --tools-multiple, --tools-single "
        ."or --tools-all (instead of $type)\n";
        usage;
      }
      push @{$tools{$type}}, openFile($ARGV[1]); shift}
    case "--tool-name" {$toolName = $ARGV[1]; shift}
    case "--threshold" {
      foreach my $type (keys %threshold) {
        $threshold{$type} = $ARGV[1];
      }
      shift}
    case (/^\-\-threshold\-/) { my ($type) = ($ARGV[0] =~ /\-\-threshold\-(.*)$/);

      if (! defined $threshold{$type}) {
        print STDERR "$type is an unknown type for the threshold\n";
        usage;
      }
      $threshold{$type} = $ARGV[1]; 
      shift}
    case "--output-ok" {$output{ok} = writeFile($ARGV[1]); shift}
    case "--output-missing" {$output{missing} = writeFile($ARGV[1]); shift}
    case "--output-false" {$output{false} = writeFile($ARGV[1]); shift}
    case "--output-tag-missing" {$output{"tag-missing"} = writeFile($ARGV[1]); shift}
    case "--output-bad-classification" {$output{"bad-classification"} = writeFile($ARGV[1]); shift}
    case "--output-expression" {$output{expression} = writeFile($ARGV[1]); shift}
    else {print STDERR "Unknown option ",$ARGV[0],"\n";}
  }
  shift;
}

my $missing;
if ((@splicing || defined $chimera) && ! defined($bed)) {
  $missing = "--bed";
}
if (defined($errors) && ! defined($fluxErr)) {
  $missing = "--flux-err";
}
if ((defined $snp || @indel)
  && (! defined($mutations) || ! defined $chrLengths)) {
  $missing = "--mutations or --chr-lengths";
}
if (@splicing && ! defined $chrLengths) {
  $missing = "--chr-lengths";
}
if  (defined($toolName) && (! defined($mutations) || ! defined($bed)
    || ! defined $fluxErr)) {
  $missing = "--mutations or --bed or --flux-err";
}
if (defined $missing) {
  print STDERR "Option(s) $missing missing\n";
  usage;
}

my @tags;
my @errors;
my @splices;
my @snp = ();
my @ins = ();
my @del = ();
my @chimera = ();
my $stats = {};
my $lineNb;
my $currentOutput;
my $hOutput;

my (%mask, %countTags, %totalTags);
my $tagNb = 0;

# Process basic tools
# Record which reads have been located (once or more) and then
# look which of them have mutations (indels, snps, chimeras, splices)
# or sequencing errors.
if (defined $toolName) {
  verbose "Processing tool";

  my @pos = ();
  my @exact;
  my @cols; 
  my $nb_reads = 0;
  my $nb_exact;

  my %readName = ();
  verbose "Reading BED file";
  # Read BED file for retrieving read names.
  while (<$bed>) {
    my ($name) = /^\S+\s+\S+\s+\S+\s+(\S+)\s+/;
    $readName{$name} = $tagNb;
    $tagNb++;
  }
  $nb_reads = $tagNb;
  $nb_exact = $nb_reads;

  foreach (my $i=0; $i < $nb_exact; $i++) {
    $exact[$i]=0;
  }

  seek $bed, 0, 0;


  my $lastTag='';
  # NOTE Loop over q(multiple single all)
  # It looks like we should not consider causes that fall into
  # a read that is located "multiple"
  foreach my $unicity (sort keys %tools) {
    foreach my $tool (@{$tools{$unicity}}) {
      verbose "Reading results";
      while (<$tool>) {
        @cols = split /\s+/;
        # Just keep reads which are located
        if (($toolName ne 'BWA' || /XT:A/)
          && (($toolName ne 'BWASW' && $toolName ne 'GSNAP' && $toolName ne 'Bowtie2')
            || $cols[2] ne '*')) {
          my $currentReadName = $cols[0];
          if ($toolName eq 'SensitivityFlux') {
            $currentReadName = $cols[2];
          }
          if ($currentReadName !~ /^[0-9]+$/ 
            && ! defined $readName{$currentReadName}) {
            print STDERR "Unknown read $currentReadName\n";
          } else { 
            my $readName;
            if ($currentReadName =~ /^[0-9]+$/) {
              $readName = $currentReadName;
            } else {
              # NOTE We should never go in there??
              $readName = $readName{$currentReadName};
            }
            # Unique or multiple
            if ($unicity eq 'single' ||
              ($unicity eq 'all' && $currentReadName ne $lastTag
                && ($toolName ne 'BWA' || /XT:A:U/)
                && ($toolName ne 'SOAP2' || $cols[3] == 1)
                && ($toolName ne 'GASSST' || /NH:i:1/))) {
              $pos[$readName] = 1; # NOTE Uniq tag localisation
            } else {
              $pos[$readName] = 2; # NOTE Multiple tag localisation
            }
            $lastTag = $currentReadName;
          }
        }
      }
      close($tool);
    }
  }
  undef %readName;

  verbose "Reading INFO file";
  # Read INFO file 
  # NOTE info fields :
  # [0] tagNumber
  # [1] chromosome
  # [2] oldPos
  # [3] newPos
  # [4] type of mutation
  # [5] length of mutation
  # [6] mutation
  # Fill info about SNPs and indels
  my %typeSeen;
  while (<$mutations>) {
    @cols = split /\s+/;
    # NOTE type can be q(ins del sub)
    my $type = $cols[4];
    if ("$type" eq  "sub") {
      $type = 'snp';
    }
    if (defined($threshold{$type})) {
      # NOTE if length of mutation is 1, the it is a snp
      # even if the type is 'ins' or 'del'
      # This line is redundant with previous test...?
      if ($cols[5] == 1) {
        $type = 'snp';
      }
      if ($type !~ /^chimera/) {
        my @tags = split /:/, $cols[0];
        # A same mutation can affect several tags, whose id are separated
        # with semicolons (:)
        foreach my $tag (@tags) {
          # NOTE It looks like a tag can only have one mutation of each type???
          if (! defined($typeSeen{$type}[$tag])) {
            $typeSeen{$type}[$tag]=1;
            $exact[$tag] = undef;
            # NOTE Not sure what we do here....????
            incStatTag($tag, $type, \@pos, $stats, \%output);
          }
        }
      }
    }
  }

  verbose "Reading BED file";
  # Read BED file to register Splice and Chimeras
  # [0] chrom
  # [1] chromStart
  # [2] chromEnd
  # [3] name
  # [4] score
  # [5] strand
  # [6] thickStart
  # [7] thickEnd
  # [8] itemRgb
  # [9] blockCount
  # [10] blockSizes
  # [11] blockStarts
  $tagNb = 0;
  while (<$bed>) {
    @cols = split /\s+/;
    if ($cols[9] >= 2) {
      my @fragLength = split /,/, $cols[11];
      my ($spliceDone, $chimeraDone) = (0, 0);
      for (my $i = 0; $i < $#fragLength; $i++) {
        # Check how many @ we have between two consecutive fragment lengths
        my $nbChimera = () = (($fragLength[$i].$fragLength[$i+1]) =~ /@/);
        # NOTE Why, nbChimera == 2?
        # If both consequitive frag contain an @, that means this is a splicing
        # on a chimeric gene!
        if ($nbChimera == 0 || $nbChimera == 2) {
          # Splice 
          if (! $spliceDone) {
            incStatTag($tagNb, 'splice', \@pos, $stats, \%output);
            $exact[$tagNb] = undef;
            $spliceDone = 1;
          }
        } elsif (! $chimeraDone){
          incStatTag($tagNb, 'chimera',\@pos, $stats, \%output);
          $exact[$tagNb] = undef;
          $chimeraDone = 1;
        }
      }
    }
    $tagNb++;
  }
  close($bed);

  verbose "Reading ERR file";
  # Read ERR file
  # Fill info about errors.
  my $previousTagNb = -1;
  while (<$fluxErr>) {
    my ($tagNb, $pos) = split /\s+/;
    if ($tagNb != $previousTagNb) {
      incStatTag($tagNb, 'errors', \@pos, $stats, \%output);
      $exact[$tagNb] = undef;
      $previousTagNb = $tagNb;
    }
  }
  close($fluxErr);

  verbose "Treating exact reads";
  for (my $i = 0; $i <= $#exact; $i++) {
    if (defined $exact[$i]) {
      incStatTag($i, 'exact', \@pos, $stats, \%output);
    }
  }
}

my %reads;
my $readCount = 0;

# NOTE Why are we counting reads now since we have open
# the bed file twice????
if (! %tools) {
  verbose "Counting reads";
# Read BED file for knowing the number of reads.
  while (<$bed>) {
    $readCount++;
  }
  seek $bed, 0, 0;
}

if (! %tools && @indel) {
  # NOTE These are bitvectors on chromosomes
  init_bit_vectors(\%mask, ['ins', 'del'], $chrLengths);

  verbose "Reading indels";

  # NOTE These are bitvectors of reads
  $reads{ins} = BitVector->new($readCount);
  $reads{del} = BitVector->new($readCount);


  my ($tag_num, $is_single, $pos, $ins, $del, $chr, $chrPos, $letter);
  foreach my $handle (@indel) {
    while (<$handle>) {
      # NOTE Here, we handle CRAC homemade outputs for indels
      if (($tag_num, $is_single, $chr, $chrPos, $pos, $ins, $del) = /^([0-9]+)\s+(\S+)\s+([^\|]+)\|\-?1,(\d+)\s+pos_indel=([0-9]+) nb_ins=(\d+) nb_del=(\d+)\s+/) {
        # NOTE this could be simplified, these both blocks are just duplicates
        if ($ins > 0) {
          if (! defined $single ||  $is_single eq 'single') {
            $countTags{ins}{$chr}{$chrPos}++;
            $totalTags{ins}++;
            $mask{ins}{$chr}->set($chrPos);
            $reads{ins}->set($tag_num);
          }
        }
        if ($del > 0) {
          if (! defined $single || $is_single eq 'single') {
            $countTags{del}{$chr}{$chrPos}++;
            $totalTags{del}++;
            $mask{del}{$chr}->set($chrPos);
            $reads{del}->set($tag_num);
          }
        }
      } elsif ( ($chr, $chrPos, $letter) = /^(\S+)\s+(\d+)\s+\*\s+([A-Z\/+-]+)\s+\d+\s+\d+/) {
        # Samtools pileup 
        # Indel
        my ($left, $right) = split /\//, $letter;
        # If it looks like +A/+A or -A/-A or */+A
        if (($left ne '*' && length $left > 2)
          || ($right ne '*' && length $right > 2)) {
          my $type;
          if ($letter =~ /\+/) {
            $type = 'ins';
          } else {
            $type = 'del';
          }
          # NOTE, here we only set the mask and do not increment others variable
          # as it is the case when parsing CRAC output
          $mask{$type}{$chr}->set($chrPos);	    
        }
      }
    }
    close($handle);
  }
  read_INFO_file({'ins' =>0, 'del' =>0}, $mutations, \%mask, \%countTags,
    \%reads, \%threshold, $stats, \%output);
  undef %{$countTags{ins}};
  undef %{$countTags{del}};
}

if (! %tools && defined $snp) {
  init_bit_vectors(\%mask, ['snp'], $chrLengths);
  $reads{snp} = BitVector->new($readCount);

  verbose "Reading SNPs";

  my ($tag_num, $is_single, $pos, $chr, $chrPos, $letter, $origLetter);
  while (<$snp>) {
    if (($tag_num, $is_single, $chr, $chrPos, $pos) = /^([0-9]+)\s+(\S+)\s+([^\|]+)\|\-?1,(\d+)\s+.*pos_SNV=([0-9]+)/) {
      # CRAC
      if (! defined $single || $is_single eq 'single') {
        $countTags{snp}{$chr}{$chrPos}++;
        $totalTags{snp}++;
        $mask{snp}{$chr}->set($chrPos);
        $reads{snp}->set($tag_num);
      }
    } elsif ( ($chr, $chrPos, $origLetter, $letter) = /^(\S+)\s+(\d+)\s+([ATGC\*])\s+([A-Z\/+-]+)\s+\d+\s+\d+\s+\d+\s+\d+\s+[^+-]+\s+/) {
      # Samtools pileup 
      if (length $letter == 1 && $letter ne $origLetter) {
        #  Substitution
        $mask{snp}{$chr}->set($chrPos);	    
      } else {
        # Indel
        my ($left, $right) = split /\//, $letter;
        # If it looks like +A/+A or -A/-A or */+A
        if (($left ne '*' && length $left == 2)
          || ($right ne '*' && length $right == 2)) {
          $mask{snp}{$chr}->set($chrPos);	    
        }
      }
    }
  }
  close($snp);
  read_INFO_file({'snp' =>0}, $mutations, \%mask, \%countTags, \%reads, 
    \%threshold, $stats, \%output);
  undef %{$countTags{snp}};
}

my $chimeraicPos;
if (! %tools && defined $chimera) {
  init_bit_vectors(\%mask, ['chimera'], $chrLengths);
  $reads{chimera} = BitVector->new($readCount);

  verbose "Reading chimeras";

  my ($tag_num, $is_single, $chr1, $chrPos1, $strand1, $chr2, $chrPos2, $strand2);
  while (<$chimera>) {
    my $found = 0;
    # #for v1.0 to v1.6 
    # if (($tag_num, $is_single, $chr1, $strand1, $chrPos1, $chr2, $chrPos2) = /^([0-9]+)\s+(\S+)\s+([^\|]+)\|(\-?1),(\d+)\s+([^\|]+)\|\-?1,(\d+)\s+/) {
    #from v1.9.9
    if (($tag_num, $is_single, $chr1, $strand1, $chrPos1, $chr2, $chrPos2) = /^([0-9]+)\s+(\S+)\s+\S+\s+\S+\s+([^\|]+)\|(\-?1),(\d+)\s+([^\|]+)\|\-?1,(\d+)\s+/) {
      # CRAC
      if (! defined $single || $is_single eq 'single') {
        $found = 1;
        $reads{chimera}->set($tag_num);
      }
    } elsif (($chr1, $chr2, $chrPos1, $chrPos2) = /^([^_]+)_(\S+)\s+(\d+)\s+(\d+)\s+/) {
      # MapSplice
      $found = 1;
      $chr1 =~ s/chr//;
      $chr2 =~ s/chr//;
    } elsif (($chr1,$chr2,$chrPos1,$chrPos2, $strand1, $strand2) = /^(.*?)\-(.*?)\s+([0-9]+)\s+([0-9]+)\s+([fr])([fr])/) {
      # TopHat Fusion | TopHat2
      $found = 1;
      $chr1 =~ s/chr//;
      $chr2 =~ s/chr//;
      $strand1 = ($strand1 eq 'r')  ? -1 : 1;
    }
    if ($found) {
      if (! $stranded && defined $strand1 && ($strand1 == -1)) {
        # Canonical chimeras
        my $tmpChr = $chr1;
        $chr1 = $chr2;
        $chr2 = $tmpChr;
        $tmpChr = $chrPos1;
        $chrPos1 = $chrPos2;
        $chrPos2 = $tmpChr;
      }
      $countTags{chimera}{$chr1}{$chrPos1."-".$chrPos2}++;
      $totalTags{chimera}++;
      $mask{chimera}{$chr1}->set($chrPos1);
      $stats->{chimera}->{count}++
      if (! defined $chimeraicPos->{$chrPos1.'@'.$chr1}->{$chrPos2});
      $chimeraicPos->{$chrPos1.'@'.$chr1}->{$chrPos2} = $chr2;
      $chimeraicPos->{$chrPos1.'@'.$chr1}->{notSeenYet}++;
    }
  }
  close($chimera);
}

my $splicePos;
if (@splicing) {
  my $count_splice = 0;
  init_bit_vectors(\%mask, ['splice'], $chrLengths);
  $reads{splice} = BitVector->new($readCount);

  verbose "Reading splices";

  my ($tag_num, $is_single, $chr, $chrPos, $gap_length, $posShift, $gapLength, $matched, $nbExons);
  foreach my $splice (@splicing) {
    while (<$splice>) {
      $matched = 0;
      # CRAC
      if (($tag_num, $is_single, $chr, $chrPos, $gap_length) = /^([0-9]+)\s+(\S+)\s+([^\|]+)\|\-?1,(\d+)\s+pos_junction=[0-9]+\s+gap_length=(\d+)\s+/) {
        if (! defined $single || $is_single eq 'single') {
          $count_splice++;
          $matched = 1;
          $countTags{splice}{$chr}{$chrPos."-".($chrPos+$gap_length)}++;
          $reads{splice}->set($tag_num);
        }
      } elsif (($chr, $chrPos, $nbExons, $posShift, $gapLength) = /^(\S+)\s+([0-9]+)\s+[0-9]+\s+\S+\s+[0-9]+\s+[+-]\s+[0-9]+\s+[0-9]+\s+[0-9,]+\s+([0-9]+)\s+([0-9,]+)\s+([0-9,]+)$/) {
        # BED file
        if ($nbExons > 1) {
          my @posShift = split /,/, $posShift;
          my @gapLength = split /,/, $gapLength;
          my $original_pos = $chrPos;
          shift @gapLength;
          $chr =~ s/^chr//;
          $matched = 0;
          $totalTags{splice}++;

          while (@gapLength) {
            $chrPos += $posShift[0];
            $gap_length =  $gapLength[0] - ($chrPos - $original_pos) + 1;

            if (! defined $splicePos->{$chrPos.'@'.$chr}->{$chrPos + $gap_length}) {
              $mask{splice}{$chr}->set($chrPos);
              $stats->{splice}->{count}++;
              $splicePos->{$chrPos.'@'.$chr}->{$chrPos + $gap_length} = 0;
              $splicePos->{$chrPos.'@'.$chr}->{notSeenYet}++;
            }

            $chrPos += $gap_length;
            shift @gapLength;
            shift @posShift;
          }
        }
      }

      if ($matched) {
        $totalTags{splice}++;
        if (! defined $splicePos->{$chrPos.'@'.$chr}->{$chrPos + $gap_length}) {
          $mask{splice}{$chr}->set($chrPos);
          $stats->{splice}->{count}++;
          $splicePos->{$chrPos.'@'.$chr}->{$chrPos + $gap_length} = 0;
          $splicePos->{$chrPos.'@'.$chr}->{notSeenYet}++;
        }
      }
    }
    close($splice);    
  }
  print "nb splice = $count_splice\n";
}

if (@splicing || defined($chimera)) {
  verbose "Reading BED file";

  my %spliceToBeSeen ;
  my %expression;
  my %causesFound;
  my %notInFP;

  $lineNb = 0;
  while (<$bed>) {
    $debugEntry = 0;
    my @line = split /\s+/;
    verbose_progress $lineNb;

    # We got a splice (or a chimera or both)
    if ($line[9] >= 2) {
      my @splicesCurrent;
      my @fragLength = split /,/,$line[10];
      my @shifts = split /,/, $line[11];
      my $sumLength = 0;
      if ($line[5] eq '+') {
        for (my $i=0; $i < $line[9]-1; $i++) {
          $sumLength += $fragLength[$i];
          $splicesCurrent[$i] = $sumLength-1;
        }
      } else {
        my $start = $line[9]-1;
        for (my $i=$start; $i > 0; $i--) {
          $sumLength += $fragLength[$i];
          $splicesCurrent[$start-$i] = $sumLength-1;
        }
      }

      # Treating chimeras. Retrieve their positions
      my $relStart = {pos => $line[1], chr => $line[0], strand => $line[5]};
      my $totalLength = $fragLength[0];
      for (my $i = 0; $i < $line[9] - 1; $i++) {
        my ($otherPos, $otherChr);
        # Position of the exon end
        my $first = giveExonicPosition($relStart, $totalLength, $fragLength[$i], $shifts[$i]);
        # Position of the next exon start
        my $second = giveExonicPosition($relStart, 0, 0, $shifts[$i+1]);

        my ($pos, $chr);
        $chr = $first->{chr};
        $pos = $first->{pos};
        ($otherChr, $otherPos) = ($second->{chr},
          $second->{pos});

        $chr =~ s/^chr//;
        my $isSplice = (! defined($first->{chimera}) && !defined($second->{chimera}))
        || (defined($first->{chimera}) && defined($second->{chimera}));

        # canonical chimera
# 		if (! $isSplice && defined $chimera) {
# 		    if ($chr gt $otherChr 
#                         || ($chr eq $otherChr && $pos > $otherPos)) {
# 			# Canonical chimeras
# 			my $tmpChr = $chr;
# 			$chr = $otherChr;
# 			$otherChr = $tmpChr;
# 			$tmpChr = $pos;
# 			$pos = $otherPos;
# 			$otherPos = $tmpChr;
# 		    }
# 		}

        if ($isSplice && @splicing) {
          $reads{splice}->unset($lineNb);
        } elsif(! $isSplice && defined($chimera)) {
          $reads{chimera}->unset($lineNb);
        }

        if ($debugEntry) {
          print STDERR "chr = $chr, pos = $pos, otherChr = $otherChr, otherPos = $otherPos\n";
        }
        if (! defined($spliceToBeSeen{$chr.'@'.$pos.'-'.$otherPos})) {
          my $result;
          # If it is a simple splice
          if ( $isSplice) {
            if (@splicing) {
              $result = findType('splice', $chr, $pos, $splicePos, 
                [$otherPos, undef, $lineNb."@".$line[3]], 
                \%mask, \%causesFound, \%notInFP, \%countTags, \%threshold, $stats, 
                \%output);
            }
          } elsif (defined($chimera)) {
            # We just do it once for each different chimera.
            $result = findType('chimera', $chr, $pos, $chimeraicPos, 
              [$otherPos, $otherChr, $lineNb."@".$line[3]], 
              \%mask, \%causesFound, \%notInFP, \%countTags, \%threshold, $stats, 
              \%output);
          }
          if (defined $result) {
            $spliceToBeSeen{$chr.'@'.$pos.'-'.$otherPos}=$result;
          } else {
            $spliceToBeSeen{$chr.'@'.$pos.'-'.$otherPos}=0;
          }
        }
        if ($spliceToBeSeen{$chr.'@'.$pos.'-'.$otherPos} == 0
          && defined($output{"tag-missing"})) {
          my $type = !$isSplice ? 'chimera':'splice';
          print {$output{"tag-missing"}} "$type\t$lineNb\n";
        } elsif ($spliceToBeSeen{$chr.'@'.$pos.'-'.$otherPos} >= 1
          && defined($output{expression})) {
          # Measure expresion of this junction
          $expression{$chr.'@'.$pos.'-'.$otherPos}->{value}++;
          $expression{$chr.'@'.$pos.'-'.$otherPos}->{type} = ($isSplice ? 'splice' : 'chimera');
        }
        $totalLength += $fragLength[$i+1];
      }
    }
    $lineNb++;
  }

  delete_in_mask(\%notInFP, \%mask);

  # Count false positives for chimera.
  verbose "";
  verbose "Counting false positives";
  if (defined $chimera) {
    count_false_positives(\%mask, $stats->{chimera}, 
      'chimera', \%output, $chimeraicPos);
    delete $mask{chimera};
    undef %$chimeraicPos;
    undef %{$countTags{chimera}};
    $stats->{chimera}->{tagsBadClassified} = $reads{chimera}->nb_set;
    if (defined($output{'bad-classification'})) {
      output_bad_classified(\%reads, 'chimera', \%output);
    }
    delete $reads{chimera};
  }

  foreach my $key (keys %spliceToBeSeen) {
    delete $spliceToBeSeen{$key};
  }

  if (@splicing) {
    count_false_positives(\%mask, $stats->{splice}, 
      'splice', \%output, $splicePos);
    delete $mask{splice};
    undef %$splicePos;
    undef %{$countTags{splice}};
    $stats->{splice}->{tagsBadClassified} = $reads{splice}->nb_set;
    if (defined($output{'bad-classification'})) {
      output_bad_classified(\%reads, 'splice', \%output);
    }
    delete $reads{splice};
  }

  if (defined ($output{expression})) {
    foreach my $junction (keys %expression) {
      print {$output{expression}} $expression{$junction}->{type},"\t",$junction,"\t",$expression{$junction}->{value},"\n";
      delete $expression{$junction};
    }
    close($output{expression});
  }
}


if ($errors && !%tools) {
  verbose "Reading ERR file";


  my ($tagNum, $pos, $seq, $numErr, $posErr);
  my $readBoth;
  my @fluxErrors = ();

  while (<$fluxErr>) {
    if (($numErr, $posErr) = /^(\d+)\s+(\d+)$/) {
      if (! defined $fluxErrors[$numErr]) {
        $fluxErrors[$numErr] = [];
      }
      push @{$fluxErrors[$numErr]}, $posErr;
      $stats->{errors}->{total}++;
    }
  }
  close($fluxErr);

  verbose "Reading errors";

  while (<$errors>) {
    if (($tagNum, $pos, $seq) = /^([0-9]+).*\s+pos_error=([0-9]+).*\s+([ACGTN]+)\s+/) {
      $totalTags{errors}++;
      $stats->{errors}->{count}++;

      my $found = $threshold{errors}+1;
      my $indexFound = -1;
      # Search the error to see if we really had it in Flux.
      if (defined $fluxErrors[$tagNum]) {
        for (my $i = 0; $i <= $#{$fluxErrors[$tagNum]}; $i++) {
          if (defined ($fluxErrors[$tagNum]->[$i])
            && $found > abs($fluxErrors[$tagNum]->[$i] - $pos)) {
            $found = abs($fluxErrors[$tagNum]->[$i] - $pos);
            $indexFound = $i;
          }
        }
      }
      if ($indexFound == -1) {
        # Not found in Flux -> false positive
        $stats->{errors}->{falsePositives}++;
        print {$output{false}} "errors $tagNum,$pos\n"
        if (defined $output{false});
      } else {
        # Found -> correct
        $stats->{errors}->{ok}++;
        print {$output{ok}} "errors $tagNum,$pos\n"
        if (defined $output{ok});
        # Delete the value from fluxErrors for not counting it in false
        # positives.
        delete $fluxErrors[$tagNum]->[$indexFound];
      }
    }
  }
  close($errors);

  # Count the errors not found
  verbose "Counting not found errors";
  for (my $i = 0; $i <= $#fluxErrors; $i++) {
    if (defined $fluxErrors[$i]) {
      foreach my $pos (@{$fluxErrors[$i]}) {
        if (defined $pos) {
          $stats->{errors}->{notFound}++;
          if (defined $output{missing}) {
            print {$output{missing}} "errors $i,$pos\n";
          }
        }
      }
      @{$fluxErrors[$i]} = undef;
    }
  }
  @fluxErrors = undef;
}



my $hasTagStat;
foreach my $category (sort keys %$stats) {
  $hasTagStat = 0;
  print "$category\n";
  foreach my $subCat (sort keys %{$stats->{$category}}) {
    if ($subCat ne 'total' && $subCat !~ /^tags/) {
      if ($subCat eq 'falsePositives' && defined $stats->{$category}->{count}
        && $stats->{$category}->{count} > 0) {
        printf "\t%s\t%d\t(%.2f%%)\n",$subCat,$stats->{$category}->{$subCat},
        $stats->{$category}->{$subCat}*100./$stats->{$category}->{count};
      } else {
        printf "\t%s\t%d\t(%.2f%%)\n",$subCat,$stats->{$category}->{$subCat},
        $stats->{$category}->{$subCat}*100./$stats->{$category}->{total};
      }
    } elsif ($subCat =~ /^tags/ && defined($totalTags{$category})) {
      printf "\t%s\t%d\t(%.2f%%)\n",$subCat, $stats->{$category}->{$subCat}, $stats->{$category}->{$subCat}*100./$totalTags{$category};
      $hasTagStat = 1;
    }
  }
  if ($hasTagStat) {
    printf "\t%s\t%d\n","tagsTotal", $totalTags{$category};
  }
  printf "\t%s\t%d\n","Total",$stats->{$category}->{total};
}

close($chrLengths)
if defined($chrLengths);
close($mutations)
if defined($mutations);
close($output{ok})
if (defined $output{ok});
close($output{false})
if (defined $output{false});
close($output{missing})
if (defined $output{missing});
