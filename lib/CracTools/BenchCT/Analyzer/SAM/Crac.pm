use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM::Crac;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer::SAM';

=head2 new

=cut

sub canCheckSnps {
  my $self = shift;
  return 1;
}

sub _checkSnps {
  my $self = shift;
  my $sam_line = shift;

  my @snps = @{$sam_line->events('SNP')};
  foreach my $snp (@snps) {
    if($self->checker->isTrueSnp($snp->{loc}->{chr},$snp->{loc}->{pos})) {
      $self->snpsStats->addTruePositive(); 
    } else {
      $self->snpsStats->addFalsePositive(); 
    }
  }

}

sub _processLine {
  my $self = shift;
  # call parent method to process the sam_line
  $self->SUPER::_processLine(@_);
  my $sam_line = shift;
  # Add a special treatment for snps if we have a snpsStats object
  $self->_checkSnps($sam_line) if defined $self->snpsStats;
}

1;
