use strict;
use warnings;

use Test::More tests => 2;
use CracTools::BenchCT::Events::Error;

my $errors = CracTools::BenchCT::Events::Error->new(
  nb_reads      => 10,
  max_length    => 10,
  sampling_rate => 10,
);


$errors->addError(0,1);
$errors->addError(3,3);
$errors->addError(3,5);
$errors->addError(3,8);
$errors->addError(7,8);
$errors->addError(8,0);

is($errors->isTrueError(3,3),2);
is($errors->isTrueError(7,8),5);

