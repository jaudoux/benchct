use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM::Crac;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer::SAM';

=head2 new

=cut

sub canCheckErrors {
  my $self = shift;
  return 1;
}

sub _checkErrors {
  my $self = shift;
  my $sam_line = shift;

  my @errors = @{$sam_line->events('Error')};
  foreach my $err (@errors) {
    next if $err->{type} =~ /Sub/i;
    my $pos = $err->{pos};
    if($sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})) {
      $pos = length($sam_line->seq) - $pos;
    }
    if($self->checker->isTrueError($sam_line->qname,$pos)) {
      $self->errorsStats->addTruePositive(); 
    } else {
      $self->errorsStats->addFalsePositive(); 
    }
  }
}

sub _processLine {
  my $self = shift;
  # call parent method to process the sam_line
  $self->SUPER::_processLine(@_);
  my $sam_line = shift;
  # Add a special treatment for snps if we have a snpsStats object
  $self->_checkErrors($sam_line) if defined $self->errorsStats;
}

1;
