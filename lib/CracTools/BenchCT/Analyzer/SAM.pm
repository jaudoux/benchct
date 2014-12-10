use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;

sub canCheckMapping {
  my $self = shift;
  return 1;
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

      # +1 because checker is 0 based
      if($self->checker->isGoodAlignment($sam_line->qname,$chr,$pos_start,$sam_line->pos+1,$strand)) {
        $self->mappingStats->addTruePositive(); 
      } else {
        $self->mappingStats->addFalsePositive(); 
        #print $sam_line->line,"\n";
      }
    } else {
      # this is a false negative
    }
  } else {
    #print $sam_line->line,"\n";
  }
}

sub _processLine {
  my $self = shift;
  my $sam_line = shift;
  $self->_checkMapping($sam_line) if defined $self->mappingStats;
}

1;
