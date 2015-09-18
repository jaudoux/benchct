use strict;
use warnings;
package CracTools::BenchCT::Analyzer::SAM::Crac;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer::SAM';

=head1 SYNOPSIS

L<SAM::Crac> analyzer is based on L<SAM> and make use of some of CRAC's extended fields
either to add filters or to check other types of events.

=head1 CAN CHECK

=over 1

=item mapping [inherited]

=item splice [inherited]

=item error

=back

L<SAM::Crac> can check I<Errors> in addition of I<Mapping> and I<Splice> already check by
L<SAM> analyzer.

=head1 OPTIONS

The C<classified> option will check is a read has the following
classification before to verify its alignement. If the classification
does not match, the read will be counted as I<false-negative>.

  classified: [unique|multiple|duplicated|normal]
  

The C<not_classified> option will check is a read has NOT the following
classification before to verify its alignement. If the classification
DOES match, the read will be counted as I<false-negative>.

  not_classified: [unique|multiple|duplicated|normal]

=head2 new

=cut

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'error') {
    return 1;
  }
  return 0;
}

sub _checkErrors {
  my $self = shift;
  my $sam_line = shift;

  my @errors = @{$sam_line->events('Error')};
  foreach my $err (@errors) {
    next if $err->{type} =~ /Sub/i;
    my $pos = $err->{pos};
    if($sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})) {
      $pos = length($sam_line->seq) - $pos - 1;
    }
    if($self->checker->isTrueError($sam_line->qname,$pos)) {
      $self->getStats('error')->addTruePositive(); 
    } else {
      $self->getStats('error')->addFalsePositive(); 
    }
  }
}

sub _checkSplice {
  my $self = shift;
  my $sam_line = shift;
  my $options = shift;

  # If there is multiple alignements, we count it as a false negative
  return 0 if defined $options->{classified} && !$sam_line->isClassified($options->{classified});
  return 0 if defined $options->{not_classified} && $sam_line->isClassified($options->{not_classified});

  return $self->SUPER::_checkSplice($sam_line,$options);
}

sub _checkMapping {
  my $self = shift;
  my $sam_line = shift;
  my $options = shift;

  # If there is multiple alignements, we count it as a false negative
  return 0 if defined $options->{classified} && !$sam_line->isClassified($options->{classified});
  return 0 if defined $options->{not_classified} && $sam_line->isClassified($options->{not_classified});

  return $self->SUPER::_checkMapping($sam_line,$options);
}

sub _processLine {
  my $self = shift;
  # call parent method to process the sam_line
  $self->SUPER::_processLine(@_);
  my $sam_line = shift;
  # Add a special treatment for snps if we have a snpsStats object
  $self->_checkErrors($sam_line) if defined $self->getStats('error');
}

1;
