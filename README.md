# What is it?

The CracTools-BenchCT software (or simply benchCT) is a powerful, yet simple,
and flexible plateform to benchmark results of softwares that aim to analyse
data produced by the Next generation sequencing technologies.

BenchCT is a qualitative evaluation tool which assess any pipeline results against
a simulated dataset produce with SimCT to obtain a clear understanding of its performance char-
acteristics in answering a particular biological question.

BenchCT is able to check the validity of various type of events: read mapping, splice junction, chimeric junction, SNV and Indel
or even sequencing error.

# Table of contents
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [What is it?](#what-is-it)
- [Table of contents](#table-of-contents)
- [Installation](#installation)
	- [Requirements](#requirements)
	- [Install from tarball](#install-from-tarball)
	- [Install from sources](#install-from-sources)
- [Documentation](#documentation)
	- [checker section](#checker-section)
		- [Files](#files)
			- [File examples](#file-examples)
		- [Thresholds](#thresholds)
	- [output section](#output-section)
	- [software section](#software-section)
		- [file entry](#file-entry)
			- [File type](#file-type)
- [Authors](#authors)

<!-- /TOC -->

# Installation

## Requirements

* Perl5 distribution
* [cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7043/lib/App/cpanminus.pm) (CPAN minus)
* [CracTools-core](https://metacpan.org/release/CracTools) perl package. It will be automaically installed by cpanm along with all other CPAN dependancies.

## Install from tarball

This is the simpliest way to install benchCT.

1. Go to the [release page](https://github.com/jaudoux/benchct/releases) of the github projet.
1. Download the latest tarball (tar.gz) release : `wget https://github.com/jaudoux/benchct/releases/download/$VERSION/CracTools-BenchCT-$VERSION.tar.gz`
2. Install the package with cpanm : `cpanm [-l local_dir] CracTools-BenchCT.tar.gz`

If you do not have admin rights, you can use the option `-l` to specify cpanm a local directory to install BenchCT.

## Install from sources

To install SimCT from the sources, you will need [Dist::Zilla]('http://dzil.org/') software, which is a Perl package manager.

1. Clone SimCT repository : `git clone https://github.com/jaudoux/benchct.git`
2. Build and install : `dzil install --install-command 'cpanm [-l LOCAL_INSTALLATION_DIRECTORY].'` (ommit the `-l` option if you want to install simCT to the system).

# Documentation

BenchCT is very simple to run, you just need to run the following command (-v options if for *verbose* mode):

    benchCT [-p nb-threads] -v benchmark.yaml > stats.tsv

The output will be printed in STDOUT, make sure to redirect it to a file if you want to save the results.

BenchCT works with a configuration file, written in
[YAML](http://www.yaml.org/) syntax that hold the whole configuration of the
benchmark you want to run.

This is an example of a benchCT configuration file:


    ---
    checker:
        files:
            infos:       GRCh38-mutated-200-48M.info
            mutations:   GRCh38-mutated-200-48M.vcf.gz
            splices:     GRCh39-mutated-200-48M-junctions.bed
            chimeras:    GRCh38-mutated-200-48M-chimeras.tsv
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

As you can see the configuration file is divided in 3 parts: **checker**,
**output**, **softwares**.

## checker section

This first section contains the informations about the
events that will be checked during this benchmark. In the `file` subsection you
can define 4 differents files for each type of events to be checked. If you use SimCT
to generate simulated dataset, these files are automatically produced.

### Files

The following table describes each of the file needed by benchCT and the event types associated with them.

file | check | Description
-----|-------|------------
infos | mapping and errors | This file holds some basic information about the simulated data (number of reads, number of errors, ...).
mutations | SNV, insertion and deletion | this file holds the mutations (SNPs and Indels) that will be checked. The file format is [VCF](http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41)
splices | splices | this file holds the splicing events that will be checked. The file format is BED, where the start and end fields define the first and last genomic positions of the intron
chimeras | chimeras | this file holds the chimeric splices that will be checked. This is not a standard format (since none exists). It is a "TAB separated value" file with columns : `chr1`, `pos1`, `strand1`, `chr2`, `pos2`, `strand2`, `read_ids`, `nb_reads`.

#### File examples

**info file**<br>

    # Date: Tue Dec 27 13:40:00 2016
    # Module: CracTools (v 1.25)
    # Script: simCT (v 0.01)
    nb_chimeras	113
    nb_errors	69914341
    nb_mutations	52514
    nb_reads	160004064
    nb_splices	132511

**chimeras file**<br>

    14      67887204        +       12      91883        +       20828040:20828041:2082804  3

### Thresholds

There is a threshold value for each type of events that can be
checked with benchCT. This values are used to validate events even if they do
not have the exacts genomic positions.

## output section

In this section you will be able to specify the kind of statistics you want
benchCT to produce, and some output configuration like the number of decimals
for real numbers.

Output format is TSV, the first line will hold the column names. The column
names depends of the `statistics` that you have chose and the events types that
are being check.

For example, if you check `mapping` and `splice` events using a SAM file and
you have choose to output `sensitivity` and `accuracy`, you will have the
following content in your output

    software    event    variable     value
    crac        mapping  Sensitivity  0.8812
    crac        mapping  Accuracy     0.9982
    crac        splice   Sensitivity  0.7812
    crac        splice   Accuracy     0.9965

Each line of the output will correspond to one of the benched softwares. If a
value is not available for a software, a `NA` will be places instead.

## software section

This last section hold the softwares output that you want
to compare with the files that you have define in the `checker` section.  Each
softwares is en entry to the software section. You can define its `name` (that
will appear in the output) and the `files` that are assossiated to this tools.
You cannot define two files that *check* the same type of events.

### file entry

For each file associated to a software you need to provide it
`name` wich is the path of the file. Its `type` and the events that will be
checked by this file (this will depend on the type). You can not check `splice`
events with a VCF file.

#### File type

This is the complete list of file types supported by benchCT
(at this time). If you want more information about a given type (including available options), you can check its man page (`man CracTools::BenchCT::Analyzer::THE_TYPE`)

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
own (see `man CracTools::BenchCT::Analyzer`) or convert your file in a
format already supported.

# Authors

Jérôme Audoux - jaudoux@cpan.org
Nicolas Philippe - nphilippe@cpan.org
Mikaël Salson - mikael.salson@univ-lille1.fr
