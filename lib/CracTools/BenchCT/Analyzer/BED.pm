use strict;
use warnings;
package CracTools::BenchCT::Analyzer::BED;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub _init {
  my $self = shift;
  my %args = @_;
  my $bed_file = $args{file};
  my $bed_it = CracTools::Utils::bedFileIterator($bed_file);
  while (my $bed_line = $bed_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($bed_line);
  }
}

# We do nothing here... childs will.
sub _processLine {
  my $self = shift;
}

1;
