use strict;
use warnings;

use Test::More tests => 3;
use CracTools::BenchCT::Analyzer;

my $analyzer = CracTools::BenchCT::Analyzer->new();

is($analyzer->canCheck('mapping'), 0);
is($analyzer->canCheck('insertion'), 0);
is($analyzer->canCheck('chimera'), 0);
