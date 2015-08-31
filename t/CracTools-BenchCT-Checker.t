use strict;
use warnings;

use Test::More tests => 1;
use CracTools::BenchCT::Checker;
use CracTools::BenchCT::Const;
use Inline::Files 0.68;
use File::Temp;

# BUG #2
{
  my $chim_file = new File::Temp( SUFFIX => '.bed', UNLINK => 1);
  while(<CHIMERA>) {print $chim_file $_;}
  close $chim_file;
  $CracTools::BenchCT::Const::THRESHOLD_CHIMERA = 10;
  my $benchCT = CracTools::BenchCT::Checker->new(chimera_tsv_file => $chim_file);
  ok($benchCT->isTrueChimera(10,48830884,-1,6,3118689,-1));
}

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

__CHIMERA__
6	3118691	-1	10	48830885	-1	8239815:8239817:8239819:8239821:8239823:8239825:8240348:8240350:8240352:8240354:8240356:8240358:8240702:8240704:8240706:8240708:8240710:8240712:8241859:8241861:8241863:8241865:8241867:8241869:8241871:8241948:8241950:8241952:8241954:8241956:8241958:8242208:8242210:8242212:8242214:8242216:8242218:8242220:8242695:8242697:8242699:8242701:8242703:8242705:8242707:8242709:8242711:8242713:8242715:8242717:8242719	51
