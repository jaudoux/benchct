use strict;
use warnings;

use Test::More tests => 1;
use CracTools::BenchCT::Checker;
use CracTools::BenchCT::Const;
use CracTools::BenchCT::Analyzer::GTF;
use Inline::Files 0.68;
use File::Temp;

{
  # First we create the "checker" with the original file
  my $original_gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<ORIGINAL_GTF>) {print $original_gtf_file $_;}
  close $original_gtf_file;
  $CracTools::BenchCT::Const::THRESHOLD_EXON = 10;
  my $checker = CracTools::BenchCT::Checker->new(gtf_file => $original_gtf_file);

  # Then we create the analyzer
  my $reconstructed_gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<RECONSTRUCTED_GTF>) {print $reconstructed_gtf_file $_;}
  close $reconstructed_gtf_file;
  my $analyzer = CracTools::BenchCT::Analyzer::GTF->new(
    file => $reconstructed_gtf_file,
    checker => $checker,
    check => 'all',
  );

  # Then we check the events
  $analyzer->run();
  is($analyzer->getStats('transcript')->nbTruePositives,1);
}

__RECONSTRUCTED_GTF__
1	StringTie	transcript	925063	935792	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	925063	925189	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	925922	926013	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "161.272552";
1	StringTie	exon	930155	930336	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "81.992920";
1	StringTie	exon	931039	931089	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "4"; cov "207.931503";
1	StringTie	exon	935772	935792	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "5"; cov "210.666672";
__ORIGINAL_GTF__
1	havana	transcript	925150	935793	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	925150	925189	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001481182"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	925922	926013	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "2"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001763717"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	925942	926013	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "2"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000393181"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	930155	930336	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00002727207"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	930155	930336	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000393181"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	931039	931089	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "4"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00002696520"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	931039	931089	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "4"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000393181"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	935772	935793	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001631320"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	935772	935793	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; exon_number "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000393181"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	UTR	925150	925189	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	UTR	925922	925941	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000437963"; transcript_version "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-003"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_end_NF"; tag "mRNA_end_NF";
2	havana	transcript	925738	944575	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
1	havana	exon	925738	925800	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00001864899"; exon_version "1";
1	havana	exon	925922	926013	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "2"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00001763717"; exon_version "1";
1	havana	CDS	925942	926013	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "2"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	start_codon	925942	925944	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "2"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
1	havana	exon	930155	930336	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002727207"; exon_version "1";
1	havana	CDS	930155	930336	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "3"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	931039	931089	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "4"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002696520"; exon_version "1";
1	havana	CDS	931039	931089	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "4"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	935772	935896	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002703998"; exon_version "1";
1	havana	CDS	935772	935896	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	939040	939129	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "6"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002686739"; exon_version "1";
1	havana	CDS	939040	939129	.	+	2	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "6"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	939275	939460	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "7"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002715021"; exon_version "1";
1	havana	CDS	939275	939460	.	+	2	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "7"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	941144	941306	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "8"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00003477353"; exon_version "1";
1	havana	CDS	941144	941306	.	+	2	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "8"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	942136	942251	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "9"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00003681266"; exon_version "1";
1	havana	CDS	942136	942251	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "9"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	942410	942488	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "10"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00003675531"; exon_version "1";
1	havana	CDS	942410	942488	.	+	2	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "10"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	942559	943058	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "11"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002728091"; exon_version "1";
1	havana	CDS	942559	943058	.	+	1	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "11"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	943253	943377	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "12"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002733131"; exon_version "1";
1	havana	CDS	943253	943377	.	+	2	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "12"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	943698	943808	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "13"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00002692620"; exon_version "1";
1	havana	CDS	943698	943808	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "13"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	exon	943908	944575	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "14"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; exon_id "ENSE00001804027"; exon_version "1";
1	havana	CDS	943908	944150	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "14"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2"; protein_id "ENSP00000342313"; protein_version "3";
1	havana	stop_codon	944151	944153	.	+	0	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; exon_number "14"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
1	havana	UTR	925738	925800	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
1	havana	UTR	925922	925941	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
1	havana	UTR	944154	944575	.	+	.	gene_id "ENSG00000187634"; gene_version "8"; transcript_id "ENST00000342066"; transcript_version "5"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-010"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS2";
