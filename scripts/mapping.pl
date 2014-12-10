#!/usr/bin/perl -w
use strict;
no strict 'refs';
use Getopt::Long;
use File::Basename;
use File::Path;

use constant TMP_DIR => $ENV{'TMPDIR'}.'crac/'.$$;

use constant BOWTIE_VERSION => "bowtie-latest";
use constant SOAP_VERSION => "soap-latest";
use constant FACTOR_VERSION => "factors.2.1";
use constant CRAC_VERSION => "crac-dev";
use constant BWA_VERSION => "bwa-latest";
use constant GASSST_VERSION => "gassst-latest";
use constant GSNAP_VERSION => "gmap-latest";
use constant TOPHAT_VERSION => 'tophat-latest';
use constant TOPHAT2_VERSION => 'tophat-latest';
use constant TOPHAT_FUSION_VERSION => 'tophatfusion-latest';
use constant MAPSPLICE_VERSION => "MapSplice-latest";
use constant GEM_VERSION => "GEM-latest";


use constant BASE => "/auto/nphilippe/transcriptome/data/indexes";
use constant GENOME_DIR => BASE."/genomes";
use constant GENOME_FASTA_DIR => "/auto/nphilippe/transcriptome/data/genome";
use constant BASERESULT => BASE."/resultat";
use constant MEMSIZEINFO => BASE.'/tools/memsizeinfo/memsizeinfo';

#use constant MOSRUN => ""; #"qsub -pe smp 12 -l exclusive=true -q smp.q -b y  "; #`which time` --format='%E' ";

sub process_param_array( $$$ );
sub process_output( $$ );			
sub strip_extension( $ );

# In the field parameters:
# - library
# - genome
# - size
# - threshold
# - nb_loc_max
# - others

my $tools = {SOAP2 => 
	     {genome_dir => 'soap',
	      genome_indexed => 1,
	      default_option => '',
	      result_dir => 'r-soap2',
	      exe => 'tools/soap/'.SOAP_VERSION.'/soap',
	      category => 'Localisation',
	      parameters => {library => '-a',
			     threads => '-p',
			     genome => ['-D','.fa.index']},
	      output => {options => 
			 {-u => 'none_located',
			  -o => 'located'},
			     err => 'summary',
			     out => undef},
		  post_processing => {script => 'statsSoap.sh',
				      args => '$pref_filename_result.located '
					  .'$pref_filename_result.none_located',
					  out_append => 'summary'}
	  },

	     Bowtie => 
	     {genome_dir => 'bowtie',
	      genome_indexed => 1,
	      default_option => '--best -n 3 -y -t -k 2',
	      result_dir => 'r-bowtie',
	      category => 'Localisation',
	      exe => 'bowtie/'.BOWTIE_VERSION.'/bowtie',
	      parameters => {size => '-l',
			     others => '$genome $library',
			     threads => '-p'},
	      options => {library => {'\.fasta$' => '-f',
				      '\.fastq$' => '-q'},
		      },
	      output => {out => 'located',
			 err => 'summary'},
	      post_processing => {script => 'statsBowtie.sh',
				  args => '$library '
				      .'$pref_filename_result.located',
				      out_append => 'summary'}
	  },

	     Bowtie2 => 
	     {genome_dir => 'bowtie2',
	      genome_indexed => 1,
	      default_option => '--end-to-end --very-sensitive -k 2',
	      result_dir => 'r-bowtie2',
	      category => 'Localisation',
	      exe => 'tools/bowtie2/'.BOWTIE_VERSION.'/bowtie2',
	      parameters => {library => '-U',
			     threads => '-p',
			     genome => '-x',
# 			     others => '-S'
			 },
	      options => {library => {'\.fasta$' => '-f',
				      '\.fastq$' => '-q'}},
              output => {'options' => 
                         {'-S' => 'located'},
                         out => 'summary',
			 err => 'summary'},
	      post_processing => {script => 'statsBowtie2.pl',
				  args => '$pref_filename_result.located',
				      out_append => 'summary'}
	  },

	     BWA => 
	     {genome_dir => 'bwa',
	      genome_indexed => 1,
	      result_dir => 'r-bwa',
	      category => 'Localisation',
	      exe => 'tools/bwa/bwaMapping.sh',
	      parameters => {
#		  threads => '-t',
		  others => '$genome.fa $library $nb_loc_max'},
	      output => {err => 'summary',
			 out => 'located'},
	      post_processing => {script => 'statsBWA.pl',
				  args => '$pref_filename_result.located',
				  out_append => 'summary'}

	  },

	     BWASW => 
	     {genome_dir => 'bwa',
	      genome_indexed => 1,
	      result_dir => 'r-bwa',
	      category => 'Localisation',
	      exe => 'tools/bwa/'.BWA_VERSION.'/bwa bwasw',
	      parameters => {threads => '-t',
			     others => '-b 5 -q2 -r1 -z10 $genome.fa $library'},
	      output => {err => 'summary',
			 out => 'located'},
	      post_processing => {script => 'statsBWASW.pl',
				  args => '$pref_filename_result.located',
				  out_append => 'summary'}

	  },

	     TopHat => 
	     {genome_dir => 'bowtie',
	      genome_indexed => 1,
	      result_dir => 'r-tophat',
	      category => 'Splice',
	      exe => 'tools/TopHat/'.TOPHAT_VERSION.'/bin/tophat',
	      parameters => {"others" => '-o $rep_result/$pref_library/".strip_extension(basename($library))." $genome $library',
			     threads => '-p',
			 }
	  },

	     TopHat2 => 
	     {genome_dir => 'bowtie2',
	      genome_indexed => 1,
	      result_dir => 'r-tophat2',
	      category => 'Splice',
	      exe => 'tools/TopHat2/'.TOPHAT2_VERSION.'/tophat',
	      parameters => {"others" => '--b2-very-sensitive --fusion-search -o $rep_result/$pref_library/".strip_extension(basename($library))." $genome $library',
			     threads => '-p',
			 }
	  },

	     TopHatFusion => 
	     {genome_dir => 'bowtie',
	      genome_indexed => 1,
	      result_dir => 'r-tophatfusion',
	      category => 'Chimera',
	      exe => 'tools/TopHatFusion/'.TOPHAT_FUSION_VERSION.'/tophat-fusion',
	      parameters => {"others" => '-o $rep_result/$pref_library/".strip_extension(basename($library))." $genome $library',
			     threads => '-p',
			     'threshold' => '--fusion-anchor-length',
			 }
	  },
	     
	     MapSplice =>
	     { genome_dir => 'mapsplice',
	       genome_indexed => 1,
	       result_dir => 'r-mapsplice',
	       category => 'Splice',
	       exe_launcher => 'python',
	       exe => 'tools/MapSplice/'.MAPSPLICE_VERSION.'/bin/mapsplice_segments.py',
	       parameters => {"others" => '--not-rerun-all -o $rep_result/$pref_library/".strip_extension(basename($library))." -c '.GENOME_FASTA_DIR.'/".basename($genome)." --fusion',
			      'genome' => '-B',
			      'size' => '-w',
			      'library' => {'\.fasta$' => '-Q fa -u',
					    '\.fa$' => '-Q fa -u',
					    '\.fastq$' => '-Q fq -u',
					    '\.fq$' => '-Q fq -u'},
			      'threshold' => '-L'
			      },
				  output => { err => 'log'}
	   },

	     GSNAP => 
	     {genome_dir => 'GSNAP',
	      genome_indexed => 1,
	      result_dir => 'r-gsnap',
	      category => 'Localisation',
	      exe => 'tools/GSNAP/'.GSNAP_VERSION.'/src/gsnap',
	      parameters => {others => "-N 1 --novel-doublesplices -A sam -D ".GENOME_DIR.'/$current_tool->{genome_dir} -d $initial_genome',
			     threads => '-t',
			     library => ' ',
			 },
	      output => {'options' => {'--split-output=' => ''}},
	      post_processing => {script => 'postProcessGSNAP.sh',
				  args => '$pref_filename_result'}
	      
	  },
	     
	     GASSST => 
	     {genome_dir => 'soap',
	      genome_indexed => 1,
	      result_dir => 'r-gassst',
	      category => 'Localisation',
	      exe => 'tools/GASSST/'.GASSST_VERSION.'/Gassst',
	      parameters => {others => '".get_gassst_options($size)."',
			     genome => ['-d','.fa'],
			     threads => '-n',
			     library => '-i',
			 },
	      output => {'options' => 
			 {'-o' => 'located'}},
	      post_processing => {script => 'postProcessGASSST.sh',
				  args => '$library $pref_filename_result.located $pref_filename_result.sam',
				  out_append => 'summary'}

	  },

	     GEM => 
	     {genome_dir => 'GEM',
	      genome_indexed => 1,
	      result_dir => 'r-gem',
	      category => 'Localisation',
	      exe => 'tools/Gem/'.GEM_VERSION.'/gem-mapper',
	      parameters => {others => '-d 2 ',
                             genome => ['-I', '.gem'],
			     threads => '-T',
			     library => '-i',
			 },
	      output => {'options' => 
			 {'-o' => 'located'}},
	  },
	     
	     CRAC =>
	     {genome_dir => 'SSA',
	      genome_indexed => 1,
	      default_option => '',
	      result_dir => 'r-crac',
	      category => 'Localisation',
	      exe => 'r-fact/'.CRAC_VERSION.'/src/crac',
	      # 	      read_input_type => {'raw'},
	      parameters => {'library' => '-r',
			     'genome' => '-i',
			     'size' => '-m',
			     'threshold' => '-k',
			     'nb_loc_max' => '-n',
			     'threads' => '--nb-threads'
			 },
	      'output' => {'options' => 
			   {'--sam' => 'sam',
			    '--snp' => 'snp',
			    '--indel' => 'indel',
			    #'--genome-indel' => 'genomeIndel',
			    '--errors' => 'seqErr',
			   # '--overlap' => 'overlap',
			    '--splice' => 'splice',
			    '--weak-splice' => 'weakSplice',
			    '--chimera' => 'chimera',
			    '--undetermined' => 'undetermined',
			    '--repeat' => 'repetition',
			    '--duplicate' => 'duplication',
			    '--nothing' => 'nothing',
			    '--normal' => 'normal',
			    '--almost-normal' => 'almostNormal',
			    '--multiple' => 'multiple',
			    '--none' => 'none',
			    '--biological' => 'bioUndetermined',
			    '--single' => 'single',
			},
			       'out' => 'summary'}
	  }
	 };

# Return options for GASSST depending on the read length
sub get_gassst_options($) {
    my ($length) = @_;
    if ($length < 75) {
	return '-w 13 -p 86 -s 3 -g 4';
    } elsif ($length < 125) {
	return '-w 18 -p 94 -s 3 -g 4';
	# return '-w 18 -p 89 -s 3 -g 6';
    } else {
	return '-w 18 -p 89 -s 3 -g 9';
    }

#     if ($length < 120) {
#         return '-w 17 -p 88 -s 3 -g 5';
#     } else {
#         return '-w 18 -p 90 -s 3 -g 9';
#     }

#     if ($length <= 60) {
# 	# For 50-bp reads
# 	return '-w 13 -p 92 -g 3 -s 3';
#     } elsif ($length <= 120) {
# 	#  For 100-bp reads
# 	return ' -w 18 -p 94 -g 4 -s 3';
#     } else {
# 	# For large reads
# 	return '-w 18 -p 95 -g 5 -s 3';
# #	return '-w 18 -p 90 -s 3 -g 9';
#     }
}

sub strip_extension($) {
    my ($filename) = @_;
    $filename =~ s/\.([a-zA-Z0-9]+)$//;
    return $filename;
}

sub usage() {
    
    my %available_genomes;
    
    my $max_len=0;
    foreach my $mappingTool (keys %$tools) {
	my $tool = $tools->{$mappingTool}->{'genome_dir'};
	if (defined($tool) && $tools->{$mappingTool}->{'genome_indexed'} == 1) {
	    opendir(my $dir, GENOME_DIR.'/'.$tool) 
		or die ("Unable to open ".GENOME_DIR.'/'.$tool.": $!\n");
	    my @files = readdir $dir;
	    foreach my $file (@files) {
		if (-r GENOME_DIR.'/'.$tool.'/'.$file 
		    && (($mappingTool ne 'GSNAP' && -f GENOME_DIR.'/'.$tool.'/'.$file)
			|| ($mappingTool eq 'GSNAP' && -d GENOME_DIR.'/'.$tool.'/'.$file))) {
		    $file =~ s/\..*$//;
		    if (length $file > 0) {
			$available_genomes{$file}{$mappingTool}=1;
			if (length $file > $max_len) {
			    $max_len = length $file;
			}
		    }
		}
	    }
	    closedir $dir;
	}
    }
    
    my $available_genomes='';
    
    foreach my $genome (sort {uc($a) cmp uc($b)} keys %available_genomes) {
	$available_genomes .= "\t\t$genome:";
	for (my $i=length $genome; $i <= $max_len; $i++) {
	    $available_genomes .= ' ';
	}
	
	my @tools;
	foreach my $tool (keys %{$available_genomes{$genome}}) {
	    push @tools, $tool;
	}
	@tools = sort @tools;
	$available_genomes .= join(',', @tools)."\n";
    }

    my $available_tools = join('; ',keys %$tools);
    
print STDERR <<USAGE;
    
    $0 -l library -r rep_result -s size -t threshold -n nb_loc_max -g genome -o option [-h threads] [--verbose]
        
	library is a file with one sequence (occurrence) by line. You have to enter the {relative,absolute} path 
	of the library (libraries repertory by default :"transcriptome/data/indexes/tags").
	size is the size of the tag.
	threshold is the threshold for the perfect match.
	genome is the genome for the matching~: 
	hg18,GRCh37 for the human genome,
	mm9 for the mouse genome,
	panTro2 for the chimpanzee genome,
	Pinot for the vine genome,
	etc...
	option = $available_tools
	threads is optional. It can be used to specify the number of threads, the program should be run on.
	verbose is optional. We can use this option to print on STDIN the detail of the procedure.

	output : the repertory of the results is rep_result if define 
	otherwise "transcriptome/data/indexes/resultat/Localisation/r-[method matching]"
	The summary files (summary.[method matching]-size.library) contain time, allocation and percentage of location.
	The results files (results.[method matching]-size.library.[type of location]) contain differents type of tag
	location.

	available genomes:
USAGE
	print STDERR 	       $available_genomes,"\n";
}

if (@ARGV < 2) {
    usage;
    exit 1;
}

# options processing
our ($help,$library,$rep_result,$size,$threshold,$nb_loc_max,$options,$genome,$verbose, $threads)= (0,'',BASERESULT,75,20,100,4,"GRCh37",'', 1);
# parse options and print usage if there is a syntax error.
GetOptions('help|?'   => \$help,
	   "l=s"       => \$library,
	   "r=s"       => \$rep_result,
	   "s=s"       => \$size,
	   "t=s"       => \$threshold,
	   "n=s"       => \$nb_loc_max,
	   'o:s'       => \$options,
	   "g=s"       => \$genome,
	   "h=s"       => \$threads,
	   "verbose"   => \$verbose
	   )
    or exit 1;

if ($help or $library eq "") {
    usage; 
    exit 1;
}

if (! defined($tools->{$options})) {
    print STDERR "Unknown tool $options\n";
    usage;
    exit 1;
}

my $command_line='';

my $current_tool = $tools->{$options};

my $pref_filename_result;
if ($rep_result ne BASERESULT) {
    $pref_filename_result = $rep_result;
} else {
    $pref_filename_result = $rep_result.'/'.$current_tool->{category}.'/'.$current_tool->{result_dir};
}
$rep_result = $pref_filename_result;

my $pref_library;
if (defined($library)) {		
    $pref_library = basename($library);
    if ($pref_library =~ /^(.*?)[_\.\-]/) {
	($pref_library) = $pref_library =~ /^(.*?)[_\.\-]/;
    }

    $pref_filename_result .= '/'.$pref_library;
}


if ($library !~ /^\//) {
    $library = BASE.'/'.$library;
}

$command_line = BASE.'/'.$current_tool->{exe};
if (defined($current_tool->{default_option})) {
    $command_line .= ' '.$current_tool->{default_option};
}

my $initial_genome = $genome;
$genome=GENOME_DIR.'/'.$current_tool->{genome_dir}.'/'.$genome;

# Processing parameters
my $current_param;
my $matched;
my $space;
foreach my $type (('options', 'parameters')) {
    foreach my $param (keys %{$current_tool->{$type}}) {
	if ($param ne 'others') {
	    $current_param = $current_tool->{$type}->{$param};
	    if (defined($current_param) 
		&& (defined($$param) || $type eq 'options')) {
		if (ref($current_param) eq 'HASH') {
		    # A hash contains regular expression to match against
		    # the variable
		    $matched = 0;
		    foreach my $exp (keys %$current_param) {
			if ($$param =~ /$exp/) {
			    $matched = 1;
			    if (ref($current_param->{$exp}) eq 'ARRAY') {
				$command_line .= ' '.process_param_array($current_param->{$exp}, $$param, $type);
			    } else {
				if ($type eq 'options') {
				    $command_line .= ' '.$current_param->{$exp};
				} else {
				    $space = ($current_param->{$exp} =~ /=$/ ? '' : ' ');
				    $command_line .= ' '.$current_param->{$exp}.$space.$$param;
				}
			    }
			}
		    }
		    if ($matched == 0) {
			print STDERR "Warning: $param $$param ignored\n";
		    }
		} elsif (ref($current_param) eq 'ARRAY') {
		    $command_line .= ' '.process_param_array($current_param, $$param, $type);
		} else {
		    if ($type eq 'options') {
			$command_line .= ' '.$current_param;
		    } else {
			$space = ($current_param =~ /=$/ ? '' : ' ');
			$command_line .= ' '.$current_param.$space.$$param;
		    }
		}
	    }
	}
    }
}

if (defined($current_tool->{parameters}->{others})) {
    my $code;
    $code = '$command_line .= " '.$current_tool->{parameters}->{others}.'"';
    eval $code;
    warn $@ if $@;
}

my $tmp_dir_result = TMP_DIR.'/'.$pref_filename_result;
my $results_dir = $pref_filename_result;
if (! -d $pref_filename_result) {
    mkpath($pref_filename_result); 
}
if (! -d $tmp_dir_result) {
    mkpath($tmp_dir_result); 
}

$genome = $initial_genome;
$pref_filename_result .= "/$options";
if (defined($library)) {
    $pref_library = basename($library);
}

foreach my $type (($genome, $pref_library, $size, $threshold)) {
    if (defined($type)) {
	$pref_filename_result .= '-'.$type;
    }
}

my $tmp_pref_filename_result = TMP_DIR.'/'.$pref_filename_result;
my $suffix;
if (defined($current_tool->{output})) {
    my $output = $current_tool->{output};
    
    if (defined($output->{options})) {
	foreach my $option (keys %{$output->{options}}) {
	    $space = ($option =~ /=$/ ? '' : ' ');
	    $suffix = length($output->{options}->{$option}) == 0 ? 
		'' : '.'.$output->{options}->{$option};
	    $command_line .= ' '.$option.$space
		.$tmp_pref_filename_result.$suffix;
	}
    }
    $command_line .= process_output($tmp_pref_filename_result, $output);
}

if (defined($current_tool->{exe_launcher})) {
    $command_line = $current_tool->{exe_launcher}.' '.$command_line;
}


if ($verbose) {
    print $command_line,"\n";
}
my $start_time = time();
# Execute the comman in another process
system($command_line);

if ($verbose) {
    print "Moving results from $tmp_dir_result to $results_dir\n";
}
system("mv ".$tmp_dir_result."/* ".$results_dir);

# my $pid = fork();
# if (! $pid) {
#     # Launch the command
#     my @command = split /\s+/, $command_line;
#     exec { $command[0] } @command or die("Unable to launch $command[0]: $!\n");
# }
# # Otherwise monitor the space usage
# my $meminfo;
# my $maxMem=0;
# my $currentMem;
# open($meminfo, MEMSIZEINFO." $pid 1000 |") or die("Unable to launch memsizeinfo: $!\n");
# while (<$meminfo>) {
#     chomp;
#     (undef, $currentMem) = split /\s+/;
#     if ($maxMem < $currentMem) {
# 	$maxMem = $currentMem;
#     }
# }
# close($meminfo);
print 'Wall-time: '.(time() - $start_time).''."\n";
# print 'Memory: '.$maxMem."\n";

$command_line='';
if (defined($current_tool->{post_processing})) {
    my $proc = $current_tool->{post_processing};
    $command_line = BASE.'/scripts/Localisation/'.$proc->{script};
    $command_line .= ' '.eval '"'.$proc->{args}.'"';
    $command_line .= process_output($pref_filename_result,$proc);


    if ($verbose) {
	print $command_line,"\n";
    }
    system($command_line);
}


sub process_param_array( $$$ ) {
    my ($array, $param, $type) = @_;
    
    if ($type eq 'options') {
	return ' '.$array->[0];
    } else {
	$space = ($array->[0] =~ /=$/ ? '' : ' ');
	return ' '.$array->[0].$space.$param.$array->[1];
    }
}

sub process_output( $$ ) {
    my ($pref, $output) = @_;
    my $result='';
    if (defined($output->{out})) {
	$result .= ' > '.$pref.'.'.$output->{out};
    }
    if (defined($output->{err})) {
	$result .= ' 2> '.$pref.'.'.$output->{err};
    }
    if (defined($output->{out_append})) {
	$result .= ' >> '.$pref.'.'.$output->{out_append};
    }
    return $result;
}

