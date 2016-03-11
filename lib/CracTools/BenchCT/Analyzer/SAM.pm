use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
 
use parent 'CracTools::BenchCT::Analyzer';

=head1 SYNOPSIS

Only Primary alignements are verified. Secondary and chimeric alignements are not verified.

=head1 CAN CHECK

=over 1

=item mapping

=item splice

=back

SAM analyzer can check mapping and splice from a SAM alignement. B<Alignement>
information is extracted from the first mapped base. B<Splice> information is
extracted from the cigar by searching C<N> characters.

=head1 OPTIONS

=over 1

=item max_hits = INT

If the C<max_hits> option is specified, we will only check reads alignements 
that have a number of hits (NH field in SAM) inferior to a given value.
Discarded reads will be counted as I<false-negative>.

=item sampling_rate = FLOAT

Fraction of reads to subsample

=back

=cut

use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;
use File::Temp qw/tempfile/;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type =~ /^(mapping|splice)$/) {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $sam_file = $args{file};
  my $options = $args{options};

  # If the subsampling option is activated
  if(defined $options->{subsample}) {
    if($sam_file =~ /\.bam$/ && -e "$sam_file.bai") {
      #my $nb_reads = _nbReadsInBam($sam_file);

      my $sub_sampling_rate = $options->{subsample};
      # Create a temp bam file wich contains a subsamples of the reads
      my ($subsample_fh, $subsample_filename) = tempfile(SUFFIX => ".bam");
      my $command_line = "samtools view -bh $sam_file -s $sub_sampling_rate > $subsample_filename && samtools index $subsample_filename";
      system($command_line);

      # Count the number of reads that have been subsampled
      my $nb_reads_subsampled = _nbReadsInBam($subsample_filename);

      # Set the new bam file for the analysis
      $self->getStats('mapping')->nbElements($nb_reads_subsampled);
      $sam_file = $subsample_filename;

      #print STDERR "Nb reads before sampling : $nb_reads, Nb reads after : $nb_reads_subsampled\n";
    } else {
      print STDERR "Cannot subsample $sam_file, only indexed BAM files can be subsampled\n";
    }
  }

  my $sam_reader = CracTools::SAMReader->new($sam_file);
  my $sam_it = $sam_reader->iterator();
  while (my $sam_line = $sam_it->()) {
    # ADD /1 or /2 to read name
    my $name = $sam_line->qname;
    if($name !~ /\/1$/ && $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT})) {
      $name .= "/1";
    } elsif($name !~ /\/2$/ && $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{LAST_SEGMENT})) {
      $name .= "/2";
    }
    $sam_line->qname($name);
    # This is a hook line for subclasses
    $self->_processLine($sam_line,$args{options});
  }
}

sub _checkMapping {
  my $self = shift;
  my $sam_line = shift;
  my $options = shift;

  if(!$sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{UNMAPPED}) # Unmapped reads counts a False Negative
    ) {

    # Remove an eventual "chr" prefix to the reference name
    my $chr = $sam_line->rname;
    $chr =~ s/^chr//;

    my ($pos_start) = $sam_line->cigar =~ /^(\d+)[HS]/;
    $pos_start = 0 if !defined $pos_start;

    # Get strand, this can be usefull if the simulated data are stranded
    my $strand = $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})? -1 : 1;

    # Set read type (0 => not PE, 1 => first pair, 2 => second pair)
    my $read_type = 0;
    $read_type = 1 if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT});
    $read_type = 2 if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{LAST_SEGMENT});

    my $good_alignment = $self->checker->isGoodAlignment($sam_line->qname,$chr,$pos_start,$sam_line->pos-1,$strand,$read_type);

    # -1 because checker is 0 based
    if($good_alignment) {
      $self->getStats('mapping')->addTruePositive(id => $good_alignment, out_string => $sam_line->line); 
    } else {
      # TODO should we add the read name as ID here?
      # IT could be interesting if the mapper yeld multiple alignements for one read and
      # do not raise the flag "secondary alignment"
      $self->getStats('mapping')->addFalsePositive(out_string => $sam_line->line);
    }
  } else {
    # this is a false negative
  }
}

sub _checkSplice {
  my $self = shift;
  my $sam_line = shift;
  my $cigar = $sam_line->cigar;
  my $pos = $sam_line->pos;
  my $chr = $sam_line->chr;
  my $last_splice_end_pos = $sam_line->pos;
  my $last_key = undef;
  my @ops = $sam_line->cigar =~ /(\d+\D)/g;
  foreach (@ops) {
    my ($nb,$op) = $_ =~ /(\d+)(\D)/;
    if($op =~ /(M|X|=|D)/) {
      $pos += $nb;
    } elsif($op eq 'N') {

      # Get strand, this can be usefull if the simulated data are stranded
      my $strand = $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})? -1 : 1;

      # TODO We should add a stranded test, when we will have stranded data
      # For now its not a good thing
      #my $strand = 1;
      #if(($sam_line->isFlagged(1) && ((!$sam_line->isFlagged(16) && $sam_line->isFlagged(64)) || ($sam_line->isFlagged(16) && $sam_line->isFlagged(128)))) # PE case
      #  || (!$sam_line->isFlagged(1) && $sam_line->isFlagged(16))) { # Single end case
      #  $strand = -1;
      #}

      # pos -1 because checker is 0-based
      my $true_splice = $self->checker->isTrueSplice($chr,$pos-1,$nb,$strand);
      
      if($true_splice) {
        $self->getStats('splice')->addTruePositive(id => $true_splice, out_string => $sam_line->line);
      } else {
        #print STDERR Dumper($bed_line);
        $self->getStats('splice')->addFalsePositive(id => "$chr;$pos", out_string => $sam_line->line);
      }
      $pos += $nb;
    }
  }
}

sub _processLine {
  my $self = shift;
  my $sam_line = shift;
  my $options = shift;

  # Only check primaru alignements,
  # skip chimeric and secondary alignemts
  # this is more fair.
  # this alignements will count as false negatives
  return 0 if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{SECONDARY_ALIGNMENT});
  return 0 if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{CHIMERIC_ALIGNMENT});

  # If the nb_hits options is activated
  # we skip aligments that have more hits than 
  # defined.
  # These alignements will count as false negatives
  my $nb_hits = $sam_line->getOptionalField("NH");
  return 0 if defined $options->{max_hits} && defined $nb_hits && $nb_hits > $options->{max_hits};

  # Now we check mapping and splice according to the
  # stats objects we have
  $self->_checkMapping($sam_line,$options) if defined $self->getStats('mapping');
  $self->_checkSplice($sam_line,$options) if defined $self->getStats('splice');
}

sub _nbReadsInBam($) {
  my $bam_file = shift;

  # Create a temp file and run samtools idxstats output
  my ($nb_reads_fh, $nb_reads_filename) = tempfile(SUFFIX => ".txt");
  my $command_line = "samtools idxstats $bam_file > $nb_reads_filename";
  system($command_line);

  # Sum the counts for each chromosome into a single integer
  my $nb_reads = 0;
  while(<$nb_reads_fh>) {
    chomp;
    my ($chr,$length,$mapped_reads,$unmapped_reads) = split "\t", $_;
    $nb_reads += $mapped_reads;
    $nb_reads += $unmapped_reads;
  }

  return $nb_reads;
}

1;
