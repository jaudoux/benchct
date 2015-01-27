use strict;
use warnings;

use Test::More tests => 1;
use CracTools::BenchCT::Checker;
use Inline::Files 0.68;
use File::Temp;

ok(1);

#{
#  my $bed_file = new File::Temp( SUFFIX => '.bed', UNLINK => 1);
#  while(<BED>) {print $bed_file $_;}
#  close $bed_file;
#  my $info_file = new File::Temp( SUFFIX => '.info', UNLINK => 1);
#  while(<INFO>) {print $info_file $_;}
#  close $info_file;
#  my $err_file = new File::Temp( SUFFIX => '.err', UNLINK => 1);
#  while(<ERR>) {print $err_file $_;}
#  close $err_file;
#  my $benchCT = CracTools::BenchCT::Checker->new(info_file => $info_file,
#                                        bed_file  => $bed_file,
#                                        err_file  => $err_file,
#                                        );
#
#  # Matching the first block
#  is($benchCT->isGoodAlignement('read1','chr2L',0,8460),1);
#  # Mathcin in the second block
#  is($benchCT->isGoodAlignement('read1','chr2L',140,8678),1);
#  # False alignement
#  is($benchCT->isGoodAlignement('read1','chr2L',140,7678),0);
#}
#
#__BED__
#chr22	1000	5000	cloneA	960	+	1000	5000	0	2	567,488,	0,3512
#chr2L	8460	8738	read1	0	-	.	.	0,0,0	2	136,64	0,214
#__ERR__
#__INFO__
