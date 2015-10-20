use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM::GSNAPTransloc;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
 
use parent 'CracTools::BenchCT::Analyzer::SAM';

use CracTools::Utils;

=head2 new

=cut

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'chimera') {
    return 1;
  }
  return 0;
}

sub _checkChimera {
  my $self = shift;
  my $sam_line = shift;

  my @fields = split(",",$sam_line->getOptionalField("XT"));
  my($strand1,$chr1,$pos1,$strand2,$chr2,$pos2) = $fields[3] =~ /([+-])(\S+)?:(\d+)?\.\.([+-])(\S+)?:(\d+)/;
  my ($nb_matched) = $sam_line->cigar =~ /(\d+)/;

  #print STDERR "strand1: $strand1\tchr1: $chr1\tpos1: $pos1\tstrand2; $strand2\tchr2: $chr2\tpos2: $pos2\n";

  my $true_chimera = $self->checker->isTrueChimera($chr1,$pos1+$nb_matched,CracTools::Utils::convertStrand($strand1),$chr2,$pos2,CracTools::Utils::convertStrand($strand2));
  
  if($true_chimera) {
    $self->getStats('chimera')->addTruePositive(id => $true_chimera);
  } else {
    #print STDERR Dumper($chimera);
    $self->getStats('chimera')->addFalsePositive(id => "$chr1;$pos1;$strand1;$chr2;$pos2;$strand2");
  }
}

sub _processLine {
  my $self = shift;
  # call parent method to process the sam_line
  $self->SUPER::_processLine(@_);
  my $sam_line = shift;
  # Add a special treatment for snps if we have a snpsStats object
  $self->_checkChimera($sam_line) if defined $self->getStats('chimera');
}

1;
