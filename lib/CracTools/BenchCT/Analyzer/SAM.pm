use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;

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
  my $sam_reader = CracTools::SAMReader->new($sam_file);
  my $sam_it = $sam_reader->iterator();
  while (my $sam_line = $sam_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($sam_line);
  }
}

sub _checkMapping {
  my $self = shift;
  my $sam_line = shift;
  # Is this is a not secondary alignement or if there is no multiple hits for this read, we do
  # consider its alignment
  if(!$sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{SECONDARY_ALIGNMENT}) &&
    !$sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{CHIMERIC_ALIGNMENT})) {

    if(!$sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{UNMAPPED}) # Unmapped reads counts a False Negative
      # && $sam_line->getOptionalField("NH") == 1
      ) {
      
      # Remove an eventual "chr" prefix to the reference name
      my $chr = $sam_line->rname;
      $chr =~ s/^chr//;

      my ($pos_start) = $sam_line->cigar =~ /^(\d+)[HS]/;
      $pos_start = 0 if !defined $pos_start;

      # Get strand, this can be usefull if the simulated data are stranded
      my $strand = $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})? -1 : 1;

      # -1 because checker is 0 based
      if($self->checker->isGoodAlignment($sam_line->qname,$chr,$pos_start,$sam_line->pos-1,$strand)) {
        $self->getStats('mapping')->addTruePositive(); 
      } else {
        $self->getStats('mapping')->addFalsePositive(); 
        #print $sam_line->line,"\n";
      }
    } else {
      # this is a false negative
    }

  } else {
    #print $sam_line->line,"\n";
  }
}

sub _checkSplice {
  my $self = shift;
  my $line = shift;
  # Is this is a not secondary alignement or if there is no multiple hits for this read, we do
  # consider its alignment
  if(!$line->isFlagged($CracTools::SAMReader::SAMline::flags{SECONDARY_ALIGNMENT}) &&
    !$line->isFlagged($CracTools::SAMReader::SAMline::flags{CHIMERIC_ALIGNMENT})) {
    my $cigar = $line->cigar;
    my $pos = $line->pos;
    my $chr = $line->chr;
    my $last_splice_end_pos = $line->pos;
    my $last_key = undef;
    my @ops = $line->cigar =~ /(\d+\D)/g;
    foreach (@ops) {
      my ($nb,$op) = $_ =~ /(\d+)(\D)/;
      if($op =~ /(M|X|=|D)/) {
        $pos += $nb;
      } elsif($op eq 'N') {
        my $strand = 1;
        if(($line->isFlagged(1) && (!$line->isFlagged(16) && $line->isFlagged(64)) || ($line->isFlagged(16) && $line->isFlagged(128))) # PE case
          || (!$line->isFlagged(1) && $line->isFlagged(16))) { # Single end case
          $strand = -1;
        }        
        my $true_splice = $self->checker->isTrueSplice($chr,$pos,$nb,$strand);
        
        if($true_splice) {
          $self->getStats('splice')->addTruePositive($true_splice);
        } else {
          #print STDERR Dumper($bed_line);
          $self->getStats('splice')->addFalsePositive();
        }
        $pos += $nb;
      }
    }
  }
}

sub _processLine {
  my $self = shift;
  my $sam_line = shift;
  $self->_checkMapping($sam_line) if defined $self->getStats('mapping');
  $self->_checkSplice($sam_line) if defined $self->getStats('splice');
}

1;
