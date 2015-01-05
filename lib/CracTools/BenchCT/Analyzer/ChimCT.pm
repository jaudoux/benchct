use strict;
use warnings;
package CracTools::BenchCT::Analyzer::ChimCT;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'chimera') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $chimera_file = $args{file};
  my $chim_it = CracTools::Utils::getFileIterator(file => $chimera_file,
    parsing_method => \&parseChimCTLine,
    header_regex => '^#'
  );
  while (my $chim_line = $chim_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($chim_line);
  }
}

# We do nothing here... childs will.
sub _processLine {
  my $self = shift;
  my $chimera = shift;
  my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
  
  if($true_chimera) {
    $self->getStats('chimera')->addTruePositive($true_chimera);
  } else {
    #print STDERR Dumper($bed_line);
    $self->getStats('chimera')->addFalsePositive();
  }
}

=head2 parseChimCTLine

=cut

sub parseChimCTLine {
  my $line = shift;
  my($id,$name,$chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$chim_value,$spanning_junction,$spanning_PE,$class,$comments,$others) = split("\t",$line);
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => $strand1,
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => $strand2,
  };
}

1;
