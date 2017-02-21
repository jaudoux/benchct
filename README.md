What is it?
==========
The CracTools-BenchCT software (or simply benchCT) is a powerful, yet simple,
and flexible plateform to benchmark results of softwares that aim to analyse
data produced by the Next generation sequencing technologies.

The philosophy behind benchCT is to assess NGS softwares results based on a
response that we already know. This means, we can evaluate softwares based on
*simulated data*, but also evaluate softwares on real data and use a database
of known variants and verify their presence in the dataset.

Therefore, we can produce statistics about a software's sensitivity and
accuracy to detect a given event.  BenchCT is able to check the validity of
various type of events: read mapping, splice alignment, SNPs and InDels calling
or even sequencing error detection.

In the future we also intend to make benchCT able to assess *differential gene
expression* or transcript/genome *reconstruction*.

Installation
============
Installing benchCT is very simple, download the tarball and install BenchCT
with the following command:

    cpanm CracTools-BenchCT.tar.gz

If you do not have admin rights, you can use the option `-l` to specify cpanm a
local directory to install benchCT.

BenchCT has a strong dependance to CracTools-core distribution wich is freely
available at CPAN. If you use `cpanm` they will be automatically installed.

Documentation
=============
BenchCT is very simple to run, you just need to run the following command (-v options if for *verbose* mode):

    benchCT -v benchmark.yaml > stats.tsv

The output will be printed in STDOUT, make sure to redirect it to a file if you want to save the results.

BenchCT works with a configuration file, written in
[YAML](http://www.yaml.org/) syntax that hold the whole configuration of the
benchmark you want to run.


This is an example of a benchCT configuration file:


    ---
    checker:
        files:
            mutations:    /data/reads/Flux/GRCh37-mutated-200-simulated-errors-48M.vcf
            mapping:      /data/reads/Flux/GRCh37-mutated-200-simulated-errors-48M.NotMutatedBed
            splices:      /data/reads/Flux/GRCh37-mutated-200-simulated-errors-48M-junctions.bed
            chimeras:    /data/reads/Flux/GRCh37-mutated-200-simulated-errors-48M-chimeras.tsv
            errors:       /data/reads/Flux/GRCh37-mutated-200-simulated-errors-48M.err
        thresholds:
            MAPPING:    5
            SNP:        5
            INSERTION:  5
            DELETION:   5
            CHIMERA:    20
            ERROR:      5
            SPLICE:     5
    output:
      statistics:
        - sensitivity
        - accuracy
        - true-positives
        - false-positives
        - nb-elements
        - gain
      nb_decimals: 4
    softwares:
        - name: soft_name
          files:
            - name: file.sam
              type: SAM
              check: 
                - mapping
              options:

As you can see the configuration file is divided in 3 parts: **`checker`**,
**`output`**, **`softwares`**.

## `checker` section

This first section contains the informations about the
events that will be checked during this benchmark. In the `file` subsection you
can define 5 differents files for each type of events to be checked.

### Files

- `mutations`: this file holds the mutations (SNPs and Indels) that will be
  checked. The file format is
  [VCF](http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41)
- `mapping`: this file holds the mapping information of each read. It is
  written in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format,
  where the read name is encoded in the *name* colum
- `splice`s: this file holds the splicing events that will be checked. The file
  format is BED, where the start and end fields define the first and last
  genomic positions of the intron
- `chimeras`: this file holds the chimeric splices that will be checked. This
  is not a standard format (since none exists). It is a "TAB separated value"
  file 
    - chr1
    - pos1
    - strand1
    - chr2
    - pos2
    - strand2
    - read_ids
    - nb_reads
- `errors`: this files holds the sequencing errors positions for each read. It
  is a TSV file with two columns, 1. the read number (1 is the read number of
  the first read of the FASTQ_1 file and 2 is the read number of the first read
  of the FASTQ_2 file, and so on...), 2 the position of the error in the read
  (starting at 1).

### Thresholds

There is a threshold value for each type of events that can be
checked with benchCT. This values are used to validate events even if they do
not have the exacts genomic positions.

## `output` section

In this section you will be able to specify the kind of statistics you want
benchCT to produce, and some output configuration like the number of decimals
for real numbers.

Output format is TSV, the first line will hold the column names. The column
names depends of the `statistis` that you have chose and the events types that
are being check.

For example, if you check `mapping` and `splice` events using a SAM file and
you have choose to output `sensitivity` and `accuracy`, you will have the
following columns in your output

    Software    mapping::sensitivity    mapping::accuracy    splice::sensitivity    splice::accuracy

Each line of the output will correspond to one of the benched softwares. If a
value is not available for a software, a `NA` will be places instead.

## `software` section

This last section hold the softwares output that you want
to compare with the files that you have define in the `checker` section.  Each
softwares is en entry to the software section. You can define its `name` (that
will appear in the output) and the `files` that are assossiated to this tools.
You cannot define two files that *check* the same type of events.

### `file` entry

For each file associated to a software you need to provide it
`name` wich is the path of the file. Its `type` and the events that will be
checked by this file (this will depend on the type). You can not check `splice`
events with a VCF file.

#### File `type`

This is the complete list of file types supported by benchCT
(at this time). If you want more information about a given type you can check its man page (`man CracTools::BenchCT::Analyzer::THE_TYPE`)

- BED::DELETION
- BED::JUNCTION
- ChimCT
- CRAC::Chimera
- MapSplice2::Chimera
- MapSplice2::InsertionA
- SAM
- SAM::Crac
- SAM::GsnapTransloc
- STAR::Chimera
- STAR::Junction
- Tophat::Fusion
- Tophat::Insertion
- VCF

If you can not find a `type` that matches your file you can either, create your
own (see `man CracTools::BenchCT::Analyzer`) or convert your file in a standard
format.

Authors
=======
Jérôme Audoux - jaudoux@cpan.org
Nicolas Philippe - nphilippe@cpan.org
Mikaël Salson - mikael.salson@univ-lille1.fr

LICENSE
=======
