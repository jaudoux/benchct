use strict;
use warnings;

use Test::More tests => 12;
use CracTools::BenchCT::Events::Chimera;

my $chim_checker = CracTools::BenchCT::Events::Chimera->new(threshold => 10);

# Class 1
my @chim1 = ("1",120,1,"2",300,-1);
my $chim1_id = $chim_checker->addChimera(@chim1);

is($chim_checker->isTrueChimera("2",300,1,"1",120,-1), $chim1_id + 1);
is($chim_checker->isTrueChimera(@chim1), $chim1_id + 1);
$chim1[1] += 5;
$chim1[4] -= 5;
is($chim_checker->isTrueChimera(@chim1), $chim1_id + 1);
$chim1[1] += 1;
ok(!$chim_checker->isTrueChimera(@chim1));

# Class 4 (same implementation as class 2,3)
my @chim2 = ("1",120,1,"1",300,-1);
my $chim2_id = $chim_checker->addChimera(@chim2);
is($chim_checker->isTrueChimera(@chim2), $chim2_id + 1);
$chim2[1] += 3;
$chim2[4] -= 7;
is($chim_checker->isTrueChimera(@chim2), $chim2_id + 1);
$chim2[1] += 1;
ok(!$chim_checker->isTrueChimera(@chim2));
# Reverse chimera is also good
is($chim_checker->isTrueChimera("1",300,1,"1",120,-1), $chim2_id + 1);
# Check false cases
is($chim_checker->isTrueChimera("1",120,-1,"1",300,1),  0);
is($chim_checker->isTrueChimera("1",300,-1,"1",120,1),  0);
is($chim_checker->isTrueChimera("1",120,-1,"1",300,-1), 0);
is($chim_checker->isTrueChimera("1",300,1,"1",220,-1),  0);
