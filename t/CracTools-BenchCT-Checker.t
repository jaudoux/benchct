use strict;
use warnings;

use Test::More tests => 4;
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
  my $benchCT = CracTools::BenchCT::Checker->new(chimera_tsv_file => $chim_file, is_stranded => 1);
  ok($benchCT->isTrueChimera(10,48830884,1,6,3118689,1));
}

# GTF and EXON/TRANSCRIPT checking
{
  my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<GTF>) {print $gtf_file $_;}
  close $gtf_file;
  $CracTools::BenchCT::Const::THRESHOLD_EXON = 10;
  my $benchCT = CracTools::BenchCT::Checker->new(gtf_file => $gtf_file);
  my @transcript;
  my $true_exon = $benchCT->isTrueExon(1,1312789,1312948,-1);
  ok($true_exon, 'trueExon (2)');
  push @transcript, $true_exon-1;
  $true_exon = $benchCT->isTrueExon(1,1319295,1319335,-1);
  ok($true_exon, 'trueExon (1)');
  push @transcript, $true_exon-1;
  ok($benchCT->isTrueTranscript(\@transcript),'trueTranscript');
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
__GTF__
1	havana	transcript	1312787	1319346	.	-	.	gene_id "ENSG00000127054"; gene_version "16"; transcript_id "ENST00000531377"; transcript_version "3"; gene_name "CPSF3L"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CPSF3L-042"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	1319296	1319346	.	-	.	gene_id "ENSG00000127054"; gene_version "16"; transcript_id "ENST00000531377"; transcript_version "3"; exon_number "1"; gene_name "CPSF3L"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CPSF3L-042"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; exon_id "ENSE00002193286"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	1314296	1314346	.	-	.	gene_id "ENSG00000127054"; gene_version "16"; transcript_id "ENST00000531377"; transcript_version "3"; exon_number "1"; gene_name "CPSF3L"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CPSF3L-042"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; exon_id "ENSE00002193286"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	1312787	1312949	.	-	.	gene_id "ENSG00000127054"; gene_version "16"; transcript_id "ENST00000531377"; transcript_version "3"; exon_number "7"; gene_name "CPSF3L"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CPSF3L-042"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; exon_id "ENSE00003654424"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
