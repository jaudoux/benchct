use strict;
use warnings;
package CracTools::BenchCT;
# ABSTRACT: Create benchmark on simulated data

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    checker => $args{checker},
    softwares => [],
  }, $class;

  return $self;
}

1;
