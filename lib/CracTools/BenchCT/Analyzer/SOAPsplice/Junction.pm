use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SOAPsplice::Junction;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'splice') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $junction_file = $args{file};
  my $splice_it = CracTools::Utils::getFileIterator(file => $junction_file,
    parsing_method => \&parseSOAPSpliceJunctionLine,
    header_regex => '^#',
  );
  while (my $splice = $splice_it->()) {
    my $true_splice = $self->checker->isTrueSplice(
      $splice->{chr},
      $splice->{start},
      $splice->{end} - $splice->{start},
      $splice->{strand},
    );
    
    if($true_splice) {
      $self->getStats('splice')->addTruePositive(id => $true_splice);
    } else {
      $self->getStats('splice')->addFalsePositive();
    }
  }
}

=head2 parseSOAPSpliceJunctionLine

1) chr: id of reference sequence that the junction comes from;
2) site1: left site of the junction site, one previous the left bound of the intron;
3) site2: right site of the junction site, one after the right round of the intron;
4) direction: the chain that the intron is on, "fwd" means it's on the direct chain, while "rev" means it's on the reverse chain.
5) number of reads supportting this junction.

=cut

sub parseSOAPSpliceJunctionLine {
  my $line = shift;
  my($chr,$start,$end,$strand,$cover) = split("\t",$line);
  return {
    chr           => $chr,
    start         => $start,
    end           => $end,
    strand        => $strand eq 'fwd'? 1 : -1,
    cover         => $cover,
  };
}

1;
