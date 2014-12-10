#!/bin/bash
## crick ##
# BASE_DIR='/auto/nphilippe/transcriptome/data/indexes'
# SCRIPT_DIR=$BASE_DIR'/scripts'
# RESULT_DIR=$BASE_DIR'/resultat'
# REFSEQ_SPLICE_DIR=$RESULT_DIR'/Splice/RefSeq'
# GENOME_DIR=$BASE_DIR'/genomes'
# REFSEQ_SNP_DIR=$GENOME_DIR/SNP/

# CRAC_GENOME_DIR=$GENOME_DIR'/SSA'

# Bowtie_TO_SAM='/auto/nphilippe/transcriptome/data/indexes/tools/samtools/samtools-latest/misc/bowtie2sam.pl'
# SOAP2_TO_SAM='/auto/nphilippe/transcriptome/data/indexes/tools/samtools/samtools-latest/misc/soap2sam.pl'

## pfloyd ##
BASE_DIR='/data'
#SCRIPT_DIR=$BASE_DIR'/projects/scripts'
SCRIPT_DIR='./'
RESULT_DIR=$BASE_DIR'/storage'
INDEX_DIR=$BASE_DIR'/indexes'
GENOME_DIR=$BASE_DIR'/genomes'
REFSEQ_SPLICE_DIR=$GENOME_DIR/SPLICE/
REFSEQ_SNP_DIR=$GENOME_DIR/SNP/

CRAC_GENOME_DIR=$INDEX_DIR'/crac'

Bowtie_TO_SAM='/usr/share/samtools/bowtie2sam.pl'
SOAP2_TO_SAM='/usr/share/samtools/soap2sam.pl'


#####
SPLICE_TOOL=(CRAC TopHat TopHat2 MapSplice TopHatFusion)
SPLICE_TopHat_DIR='/Splice/r-tophat'
SPLICE_TopHat_OUT_DIR='/Splice/r-tophat'
SPLICE_TopHat_INPUT_FILES='junctions.bed'
#---
SPLICE_TopHat2_DIR='/Splice/r-tophat2'
SPLICE_TopHat2_OUT_DIR='/Splice/r-tophat2'
SPLICE_TopHat2_INPUT_FILES='junctions.bed'
#---
SPLICE_TopHatFusion_DIR='/Chimera/r-tophatfusion'
SPLICE_TopHatFusion_OUT_DIR='/Chimera/r-tophatfusion'
SPLICE_TopHatFusion_INPUT_FILES='junctions.bed'
#---
SPLICE_CRAC_DIR='/Mapping/r-crac'
SPLICE_CRAC_OUT_DIR='/Splice/r-crac'
SPLICE_CRAC_INPUT_FILES='CRAC*${LIB_NAME}*.{splice,coverlessSplice}'
#---
SPLICE_GSNAP_DIR='/Mapping/r-gsnap'
SPLICE_GSNAP_OUT_DIR='/Splice/r-gsnap'
SPLICE_GSNAP_INPUT_FILES='GSNAP*${LIB_NAME}*.bed'
#---
SPLICE_MapSplice_DIR='/Splice/r-mapsplice'
SPLICE_MapSplice_OUT_DIR='/Splice/r-mapsplice'
SPLICE_MapSplice_INPUT_FILES='best_junction.bed'
#---
SPLICE_SCRIPT=$SCRIPT_DIR'crossSpliceCRAC-Db.pl' # Only for real data
SPLICE_JUNCTION_UPPER_BOUNDS="1 3 5 7 10"
SPLICE_SCRIPT_OUTPUT='splices.'

######
MAPPING_TOOLS=(Bowtie Bowtie2 BWA BWASW SOAP2 GASSST GSNAP TopHat TopHat2 MapSplice TopHatFusion TopHatFusionPost)
TOOLS_Bowtie_DIR='/Mapping/r-bowtie'
TOOLS_Bowtie_BASENAME='Bowtie*${LIB_NAME}*.'
TOOLS_Bowtie_INPUT_FILES='Bowtie*${LIB_NAME}*-${THRESHOLD}*.located'
TOOLS_Bowtie_OUTPUT_DIR='/CrossMappingTools/CrossCracBowtie'
TOOLS_Bowtie_OUT_FILE='statsCRACvsBowtie'
TOOLS_Bowtie_GENOME_DIR=$GENOME_DIR/bowtie
#---
TOOLS_Bowtie2_DIR='/Mapping/r-bowtie2'
TOOLS_Bowtie2_BASENAME='Bowtie2*${LIB_NAME}*.'
TOOLS_Bowtie2_INPUT_FILES='Bowtie2*${LIB_NAME}*-${THRESHOLD}*.located'
TOOLS_Bowtie2_OUTPUT_DIR='/CrossMappingTools/CrossCracBowtie2'
TOOLS_Bowtie2_OUT_FILE='statsCRACvsBowtie2'
TOOLS_Bowtie2_GENOME_DIR=$GENOME_DIR/bowtie2
#---
TOOLS_BWA_DIR='/Mapping/r-bwa'
TOOLS_BWA_BASENAME='BWA-*${LIB_NAME}*.'
TOOLS_BWA_INPUT_FILES='BWA-*${LIB_NAME}*-${THRESHOLD}*.located'
TOOLS_BWA_OUTPUT_DIR='/CrossMappingTools/CrossCracBWA'
TOOLS_BWA_OUT_FILE='statsCRACvsBWA'
TOOLS_BWA_GENOME_DIR=$GENOME_DIR/bwa
#---
TOOLS_BWASW_DIR='/Mapping/r-bwa'
TOOLS_BWASW_BASENAME='BWASW*${LIB_NAME}*.'
TOOLS_BWASW_INPUT_FILES='BWASW*${LIB_NAME}*-${THRESHOLD}*.located'
TOOLS_BWASW_OUTPUT_DIR='/CrossMappingTools/CrossCracBWASW'
TOOLS_BWASW_OUT_FILE='statsCRACvsBWASW'
TOOLS_BWASW_GENOME_DIR=$GENOME_DIR/bwa
#---
TOOLS_SOAP2_DIR='/Mapping/r-soap2'
TOOLS_SOAP2_BASENAME5B='SOAP2*${LIB_NAME}*.'
TOOLS_SOAP2_INPUT_FILES='SOAP2*${LIB_NAME}*-${THRESHOLD}*.located'
TOOLS_SOAP2_OUTPUT_DIR='/CrossMappingTools/CrossCracSOAP2'
TOOLS_SOAP2_OUT_FILE='statsCRACvsSOAP2'
TOOLS_SOAP2_GENOME_DIR=$GENOME_DIR/soap
#---
TOOLS_CRAC_DIR='/Mapping/r-crac'
TOOLS_CRAC_BASENAME='CRAC*${LIB_NAME}*-${THRESHOLD}*.'
TOOLS_CRAC_INPUT_FILES='CRAC*${LIB_NAME}*-${THRESHOLD}*.{single,duplication,multiple,none,almostNormal,bioUndetermined,indel,splice,chimera,normal,nothing,repetition,error,snp,coverlessSplice,undetermined}'
TOOLS_CRAC_SINGLE_LOCATED='CRAC*${LIB_NAME}*-${THRESHOLD}*.sam'
TOOLS_CRAC_GENOME_DIR=$INDEX_DIR/crac
#TOOLS_CRAC_GENOME_DIR=$GENOME_DIR/SSA
#---
TOOLS_GASSST_DIR='/Mapping/r-gassst'
TOOLS_GASSST_INPUT_FILES='GASSST*${LIB_NAME}*.sam'
TOOLS_GASSST_BASENAME='GASSST*${LIB_NAME}*-${THRESHOLD}*.'
TOOLS_GASSST_OUTPUT_DIR='/CrossMappingTools/CrossCracGASSST'
TOOLS_GASSST_OUT_FILE='statsCRACvsGASSST'
TOOLS_GASSST_GENOME_DIR=$GENOME_DIR/soap
#---
TOOLS_GSNAP_DIR='/Mapping/r-gsnap'
TOOLS_GSNAP_INPUT_FILES='GSNAP*${LIB_NAME}*.unpaired_{mult,uniq}'
TOOLS_GSNAP_BASENAME='GSNAP*${LIB_NAME}*-${THRESHOLD}*.'
TOOLS_GSNAP_OUTPUT_DIR='/CrossMappingTools/CrossCracGSNAP'
TOOLS_GSNAP_OUT_FILE='statsCRACvsGSNAP'
TOOLS_GSNAP_SINGLE_LOCATED='GSNAP*${LIB_NAME}*-${THRESHOLD}*.unpaired_uniq'
TOOLS_GSNAP_GENOME_DIR=$GENOME_DIR/GSNAP

# These tools are not used with $TOOLS_SCRIPT but this information is necessary
# in the following.
TOOLS_TopHat_DIR=$SPLICE_TopHat_DIR
TOOLS_TopHat_BASENAME=''
TOOLS_TopHat_INPUT_FILES=$SPLICE_TopHat_INPUT_FILES
#---
TOOLS_TopHat2_DIR=$SPLICE_TopHat2_DIR
TOOLS_TopHat2_BASENAME=''
#TOOLS_TopHat2_INPUT_FILES="$SPLICE_TopHat2_INPUT_FILES fusions.out"
TOOLS_TopHat2_INPUT_FILES='${LIB_NAME}/accepted_hits.sam'
#---
TOOLS_TopHatFusion_DIR=$SPLICE_TopHatFusion_DIR
TOOLS_TopHatFusion_BASENAME=''
TOOLS_TopHatFusion_INPUT_FILES="$SPLICE_TopHatFusion_INPUT_FILES fusions.out"
#---
TOOLS_TopHatFusionPost_DIR=$SPLICE_TopHatFusion_DIR
TOOLS_TopHatFusionPost_BASENAME=''
TOOLS_TopHatFusionPost_INPUT_FILES="$SPLICE_TopHatFusion_INPUT_FILES fusions.result.out"
#---
TOOLS_MapSplice_DIR=$SPLICE_MapSplice_DIR
TOOLS_MapSplice_BASENAME=''
TOOLS_MapSplice_INPUT_FILES='${LIB_NAME}/alignments.sam'

# We do not care about this one!
TOOLS_SCRIPT=$SCRIPT_DIR'/crossCRACvsMappingTools.pl'

#####
SNP_CRAC_OUT_DIR='/SNP/r-crac'
#--
SNP_Bowtie_OUT_DIR='/SNP/r-bowtie'
#--
SNP_Bowtie2_OUT_DIR='/SNP/r-bowtie2'
#--
SNP_GASSST_OUT_DIR='/SNP/r-gassst'
#--
SNP_GSNAP_OUT_DIR='/SNP/r-gsnap'
#--
SNP_SOAP2_OUT_DIR='/SNP/r-soap2'
#--
SNP_BWA_OUT_DIR='/SNP/r-bwa'
#####
Splices_CROSSED_JUNCTION=3
SPLICES_CROSSED_Bowtie_FILES='CRAC.interBowtie*${LIB_NAME}*_{splice,coverlessSplice}'
SPLICES_CROSSED_Bowtie_OUT_FILE='${PREFIX}_splicesFound'
SPLICES_CROSSED_Bowtie_EXTRACT_PREFIX='%_*'
#---
Splices_CROSSED_JUNCTION=3
SPLICES_CROSSED_Bowtie2_FILES='CRAC.interBowtie2*${LIB_NAME}*_{splice,coverlessSplice}'
SPLICES_CROSSED_Bowtie2_OUT_FILE='${PREFIX}_splicesFound'
SPLICES_CROSSED_Bowtie2_EXTRACT_PREFIX='%_*'
#---
SPLICES_CROSSED_BWA_FILES='CRAC.interBWA*${LIB_NAME}*_{splice,coverlessSplice}'
SPLICES_CROSSED_BWA_OUT_FILE='${PREFIX}_splicesFound'
SPLICES_CROSSED_BWA_EXTRACT_PREFIX='%_*'
#---
SPLICES_CROSSED_SCRIPT=$SCRIPT_DIR'/crossRefseqBowtieSplice.pl' # We do not care
#####
SNP_CREATE_INDEX="samtools faidx "
SNP_SORTED_BAM_EXT='.sorted'
SNP_CREATE_SORTED_BAM='samtools view -buSt @INDEX_FAI @SAM_FILE | samtools sort - @SAM_SORTED_BAM'
SNP_CREATE_PILEUP='samtools mpileup -ugf @INDEX_FA @SAM_SORTED_BAM.bam | bcftools view -bvcg - > @SAM_FILE.bcf &&  bcftools view @SAM_FILE.bcf | vcfutils.pl varFilter -D 100 | grep -v ^# | sort -k1,1 -k2n,2 > @SAM_FILE.vcf'
#####
CROSS_SNP_REFSEQ_SCRIPT=$SCRIPT_DIR'/crossRefSeqSNP.pl' # For real data
CRAC_SNP_EXTRACTION_SCRIPT=$SCRIPT_DIR'extractSNPfromCRAC.pl'
DEFAULT_THRESHOLD_SNP=5
#####
CHIMERA_SCRIPT=$SCRIPT_DIR'postProcessChimera.pl'
CHIMERA_NB_BASES_TOLERANCE=50
CHIMERA_INPUT_FILES=$TOOLS_CRAC_BASENAME'chimera'
CHIMERA_BOWTIE_GENOME_DIR=$TOOLS_Bowtie_GENOME_DIR
CHIMERA_REFSEQ_OUTPUT_FILE='chimera.'
CHIMERA_REFSEQ_GOOD_CATEGORIES='(overlap)|(splice)|(chimera)'
CHIMERA_OUTPUT_DIR='/Chimera/r-crac'
CHIMERA_OUTPUT_FILE='${PREFIX}.filtered.chimera'
CHIMERA_OUTPUT_FILE_REFSEQ='${PREFIX}.filtered.refseq.chimera'
CHIMERA_DEFAULT_FLAG=15
#####

#####
FLUX_ERRORS_EXTENSION="err"
FLUX_MUTATION_EXTENSION="info"
FLUX_NO_MUTATION_BED_EXTENSION="NotMutatedBed"
FLUX_CROSS_OUTPUT_DIR='/CrossMappingTools/Cross${TOOL_NAME}Flux/'
FLUX_TAGS_NOT_FOUND_OUTPUT='tags.notFound'
FLUX_TAGS_STAT_OUTPUT='tags.stats'
FLUX_TAGS_FP_STAT_OUTPUT='tags.falsePositives.stats'
FLUX_CAUSES_NOT_FOUND_OUTPUT='causes.notFound'
FLUX_CAUSES_FP_OUTPUT='causes.falsePositives'
FLUX_CAUSES_TAG_MISSING_OUTPUT='causes.tagMissing'
FLUX_CAUSES_BAD_CLASSIFICATION_OUTPUT='causes.badClassif'
FLUX_CAUSES_STAT_OUTPUT='causes.stats'

FLUX_CAUSES="splicing snp chimera errors indel" # snp indel chimera errors"
#---
FLUX_CAUSES_CRAC=1
FLUX_CAUSES_CRAC_splicing="splice coverlessSplice"
FLUX_CAUSES_CRAC_snp="snp"
FLUX_CAUSES_CRAC_indel="indel"
FLUX_CAUSES_CRAC_chimera="chimera"
FLUX_CAUSES_CRAC_errors="error"
FLUX_CAUSES_CRAC_noDup="single"
#---
FLUX_CAUSES_TopHat=1
FLUX_CAUSES_TopHat_splicing=junctions.bed
#---
FLUX_CAUSES_TopHat2=1
FLUX_CAUSES_TopHat2_splicing=junctions.bed
FLUX_CAUSES_TopHat2_chimera=fusions.out
#---
FLUX_CAUSES_MapSplice=1
FLUX_CAUSES_MapSplice_splicing=best_junction.bed
FLUX_CAUSES_MapSplice_chimera=fusion.junction
#---
FLUX_CAUSES_TopHatFusion=1
FLUX_CAUSES_TopHatFusion_splicing=junctions.bed
FLUX_CAUSES_TopHatFusion_chimera=fusions.out
#---
FLUX_CAUSES_TopHatFusionPost=1
FLUX_CAUSES_TopHatFusionPost_chimera=fusions.result.out
#---
FLUX_CAUSES_BWA=1
FLUX_CAUSES_BWA_snp=pileup
FLUX_CAUSES_BWA_indel=pileup
#---
FLUX_CAUSES_BWASW=1
FLUX_CAUSES_BWASW_snp=pileup
FLUX_CAUSES_BWASW_indel=pileup
#---
FLUX_CAUSES_Bowtie=1
FLUX_CAUSES_Bowtie_snp=sam.pileup
FLUX_CAUSES_Bowtie_indel=sam.pileup
#---
FLUX_CAUSES_Bowtie2=1
FLUX_CAUSES_Bowtie2_snp=sam.pileup
FLUX_CAUSES_Bowtie2_indel=sam.pileup
#---
FLUX_CAUSES_GASSST=1
FLUX_CAUSES_GASSST_snp=pileup
FLUX_CAUSES_GASSST_indel=pileup
#---
FLUX_CAUSES_GSNAP=1
FLUX_CAUSES_GSNAP_splicing=bed
#---
FLUX_CAUSES_SOAP2=1
FLUX_CAUSES_SOAP2_snp=sam.pileup
FLUX_CAUSES_SOAP2_indel=sam.pileup
#---
FLUX_CAUSES_SCRIPT=$SCRIPT_DIR'/crossCausesTagsFlux-Tools.pl'
FLUX_TAGS_SCRIPT=$SCRIPT_DIR'/crossCausesTagsFlux-Tools.pl'
#####

#####
ACCURACY_SCRIPT=$SCRIPT_DIR'/sensitivity-Flux.pl'
ACCURACY_STATS_EXT=.accuracy.stats
ACCURACY_STATS_FALSE_POS=.accuracy.falsePositive
#ACCURACY_STATS_TRUE_POS=.accuracy.truePositive
#####


LIBRARY_DIR=$BASE_DIR'/reads'

# Gets the available genomes for the RefSeq
AVAILABLE_GENOMES='';
for i in $REFSEQ_SPLICE_DIR/*
do
  if [ -d "$i" ]
      then
      NB_FILES=0
      AVAILABLE_FILE=''
      for file in $i/*.gff; do
#	  echo $file
	  TMP=`basename $file`
	  AVAILABLE_FILE=$AVAILABLE_FILE${TMP%%.gff}" "
	  NB_FILES=$((NB_FILES+1))
      done

      if [ $NB_FILES -eq 1 ]; then
	  AVAILABLE_GENOMES=$AVAILABLE_GENOMES`basename $i`" "
      elif [ $NB_FILES -ge 2 ]; then 
	  AVAILABLE_GENOMES=$AVAILABLE_GENOMES`basename $i`" ("$AVAILABLE_FILE"\b)\n\t\t"
      fi
  fi
done


if [ $# -lt 3 ]
then
    echo -e "Usage: $0 -l library -g genome -t threshold
         --full-lib <name> --refseq <file>
         --accept-threshold <threshold>
         --library-dir <dir> --snp <tool> --indel CRAC
         --splices <tool> --tools <tool> --splices-within-mapped <tool>
         --chimera --flag <flag> --flux-causes <tool> --flux-tags <tool>
         --flux-tags-fp <tool> --accuracy <tool> --no-dup" >&2

    if [ "$1" == "-h" -o "$1" == "--help" ]
	then
	echo -e "

  Cross the results obtained by CRAC with other tools.

  library: Name of the tags that were used for the mapping.
           Ex: K562, DrosoSimulated, SRR034309

  genome: Name of the genome for crossing with splices (options --splices 
          or --splices-within-mapped or --chimera), SNPs (and several others)

  threshold: k-mer size used for CRAC. Mandatory parameter
             with the option --splices (among others).
  
  --full-lib: Specify the full name of the library when the option -l
             may lead to an ambiguity (if you have several librairy names
             sharing the same prefix).

  --library-dir: Specify the directory in which are stored the reads.
            The directory can be specified relatively to the current
            directory or relatively to $LIBRARY_DIR.
            This option is mandatory with --flux-causes.

  --refseq: Specify the name of the RefSeq file to use (whether relative to the
            specified genome or the absolute path to the file).
            The RefSeq is used to check splices against known exons.
              Available: $AVAILABLE_GENOMES
            It is also used to check SNPs against known SNPs.

  --snp: Computes the SNPs and indels for the specified tool (one which has a
         SAM output) using SamTools.
         If the --refseq option is specified, it will compute the number of SNPs
         found that are in the RefSeq data.

  --indel: Juste used with CRAC (for now?). It is used for crossing indels with the
           RefSeq data from dbSNP (or elsewhere). Hence the --refseq option must
           be specified with that option.

  --accept-threshold: Specify the threshold when crossing splices or SNPs with --flux-causes
            or --snp (for example) to determine until which tolerance threshold
            a given cause should be accepted as correct.

  --splices: Cross the results from the specified tool (${SPLICE_TOOL[*]}), for
             a given library, using the RefSeq of the specified genome.

  --tools: Cross the results from the specified tool (${MAPPING_TOOLS[*]}), with
           CRAC results for the specified library.

  --splices-within-mapped: Cross the results between splices discovered by CRAC
           and the cross processed using the --tools option.
           You must specify the tool that was specified with the --tools option.

  --chimera: Validate the chimera found by CRAC for the library given in parameter.
             A first postprocessing step consists in removing chimeras with short breaks,
             low support, duplications, ... You may specify which chimera you want to 
             keep by specifying the --flag option.
             An optional second postprocessing step consists in checking if filtered 
             chimeras are on known genes boundaries (for doing that step, you 
             must specify the --refseq option).

  --flag: Determines which chimera one wants to keep with the filtration done by 
          --chimera. For being very specific you should specify a value of 1
          (default: ".$CHIMERA_DEFAULT_FLAG.").
          For more details, see the help in $CHIMERA_SCRIPT.

  --flux-causes: Count the number of SNP, indels, chimera, splices
             found by the specified tool (CRAC, TopHat, TopHat2, TopHatFusion, TopHatFusionPost,
             GSNAP, MapSplice).
             This is mainly meaningful for CRAC, TopHat{,2} and MapSplice since
             the others do not classify the reads. For TopHat
             we just count the number of splices found.

  --flux-tags: 
             It counts the percentage of tags located
             with a given cause (eg. 5% of tags with SNP are located).
             Tools can be one of ${MAPPING_TOOLS[*]}.

  --flux-tags-fp:
             Same as the previous one but counts the percentage of false
             positives instead. The option --accuracy must have been
             launched before using the option --flux-tags-fp.

  --accuracy: Computes the accuray of the single located reads of the specified
             tool (CRAC ${MAPPING_TOOLS[*]}).

  --no-dup: Option that can be given with --flux-causes to avoir computing
            FP, TP for reads that are considered as duplicated.
            That will mechanically lower the percentage of TP, but
            (hopefully) the percentage of FP." >&2
    fi
    exit 1
fi


verbose() {
    echo -e $*
}
fileFound() {
    FOUND=0
    for file in  $*
      do
      if [ -f $file -a -r $file ]
	  then
	  FOUND=1
      fi
    done
    return $FOUND
}

LIBRARY=''
GENOME=''
SPLICES=''
TOOLS=''
SPLICES_CROSSED=''
THRESHOLD=''
LIB_NAME=''
CHIMERA=''
REFSEQ_FILE=''
SNP=''
LIB_DIR=''
CAUSES=''
FLUX_TAGS=''
ACCURACY=''
THRESHOLD_ACCEPT=''
CAUSE_SPECIF=''
NO_DUP=''
INDEL=''
FLUX_TAGS_FP=''
FLAG=$CHIMERA_DEFAULT_FLAG

while [ $# -ne 0 ];
do
  case "$1" in
      -l) LIBRARY=$2; shift 2;;
      -g) GENOME=$2; shift 2;;
      -t) THRESHOLD=$2; shift 2;;
      --full-lib) LIB_NAME=$2; shift 2;;
      --library-dir) LIB_DIR=$2; shift 2;;
      --accept-threshold) THRESHOLD_ACCEPT=$2; shift 2;;
      --splices) SPLICES=$2; shift 2;;
      --tools) TOOLS=$2; shift 2;;
      --snp) SNP=$2; shift  2;;
      --indel) INDEL=$2; shift 2;;
      --splices-within-mapped) SPLICES_CROSSED=$2; shift 2;;
      --chimera) CHIMERA=1; shift 1;; 
      --flag) FLAG=$2; shift 1;;
      --refseq) REFSEQ_FILE=$2; shift 2;;
      --flux-causes) CAUSES=$2; shift 2;;
      --flux-tags) FLUX_TAGS=$2; shift 2;;
      --flux-tags-fp) FLUX_TAGS_FP=$2; shift 2;;
      --accuracy) ACCURACY=$2; shift 2;;
      --cause-specificity) CAUSE_SPECIF=1; shift 1;;
      --no-dup) NO_DUP='.noDuplication'; shift 1;;
      *) echo "Unknown option $1, ignored." >&2; shift 1;;
  esac
done

# Check functions

record_result() {
    if [ $# -eq 2 -a ! -z "$1" ]
	then
	eval "$1=\"$2\""
    fi
}

# Check if the genome parameter has been given.
check_genome() {
    if [ -z "$GENOME" ]
	then
	echo -e "Please specify a genome ($AVAILABLE_GENOMES)" >&2
	exit 13
    fi
}

# Check if there exists a refseq for the given genome. 
check_genome_refseq() {
    check_genome
    TYPE=$1
    eval DIR='$'REFSEQ_$TYPE'_DIR/$GENOME'
    if [ ! -d "$DIR" ]
	then
	echo "Invalid genome '$GENOME': $DIR does not exist" >&2
	exit 2
    fi
    record_result "$2" "$DIR"
}

# Check if a given tool has been launched on the given library.
check_library() {
    eval local DIR=$1
    if [ -z "$LIBRARY" ] || eval [ ! -d $RESULT_DIR$DIR/$LIBRARY ]
	then
	eval echo "Invalid library '$LIBRARY': $RESULT_DIR$DIR/$LIBRARY does not exist" >&2
	exit 4
    fi
    record_result "$2" "$RESULT_DIR$DIR/$LIBRARY"
}

#Check if the threshold has been given.
check_threshold() {
    if [ -z "$THRESHOLD" ]
	then
	echo "Invalid threshold" >&2
	exit 5
    fi    
}

#Check if the first parameter is valid name for a refseq.
# If it is record the path in a variable.
check_refseq() {
    verbose "Searching the appropriate RefSeq file"

    local TYPE=$1
    eval local DIR='$REFSEQ_'$TYPE'_DIR/$GENOME'
    local REFSEQ=$2
    if [ -z "$2" ]
	then
	verbose "No --refseq specified, trying to get default one"
	REFSEQ=`echo $DIR/*.gff`
    elif [ -f "$DIR/$REFSEQ" ]
	then
	REFSEQ=$DIR/$REFSEQ
    fi
    if [ ! -f "$REFSEQ" -o ! -r "$REFSEQ" ]
	then
	echo "I didn't find the appropriate RefSeq file or I can't read it: $REFSEQ" >&2
	exit 6
    fi
    verbose "\t-> $REFSEQ"
    record_result "$3" $REFSEQ
}

# Check if an output dir already exists or create it otherwise.
check_outdir() {
    eval local DIR=$1
    if eval [ ! -d $DIR ]
	then
	verbose "Create directory $DIR"
	eval mkdir -p $DIR
    fi

}

# Check if the library has been launched on CRAC
# and record the result in the first parameter.
check_CRAC_library() {
    if [ ! -d $RESULT_DIR$TOOLS_CRAC_DIR/$LIBRARY ]
	then
	echo "CRAC has not been launched on $LIBRARY" >&2
	exit 9
    fi
    record_result "$1" $RESULT_DIR$TOOLS_CRAC_DIR/$LIBRARY
}

# Check if the library name given is a valid library name
# and record it in the variable whose name is the second parameter.
check_libdir() {
    if [ -z "$LIB_DIR" ]
	then
	echo "--library-dir has not been specified. Aborting." >&2
	exit 10
    fi
    if [ ! -d $LIB_DIR ]
	then
	if [ ! -d $LIBRARY_DIR/$LIB_DIR ]
	    then
	    echo "$LIB_DIR is not a valid library directory" >&2
	    exit 21
	else
	    tmp_dir=$LIBRARY_DIR/$LIB_DIR
	fi
    else
	tmp_dir=$LIB_DIR
    fi
    record_result "$1" $tmp_dir
}

# Check if the first parameter is a file and if it is, sets a variable
# whose name is the second parameter.
fileFoundSet() {
    eval eval fileFound $1
    if [ $? -eq 1 ]
	then
	local FILENAME=`eval eval echo $1`
	record_result "$2" "$FILENAME"
    else
	eval eval echo "File $1 not found" >&2
	exit 15
    fi
}

# Check if the first parameter is a directory and is readable.
check_dir() {
    if [ ! -d "$1" -o ! -r "$1" ]
	then
	echo "$1 cannot be accessed" >&2
	exit 16
    fi
    record_result "$2" "$1"
}

# Check if it is a valid tool
check_tool() {
    if eval [ -z '$'TOOLS_${1}_DIR ]
	then
	echo "Invalid tool $1" >&2
	exit 7
    fi
}

# Check if the results exist for a given tool and record the path
# in the variable given in parameter.
get_results() {
    local TOOL=$1

    if eval [ a'$'TOOLS_${TOOL}_SINGLE_LOCATED == "a" ]
	then
	local eval eval INPUT=$RESULT_DIR'$'TOOLS_${TOOL}_DIR/$LIBRARY/'$'TOOLS_${TOOL}_INPUT_FILES
    else
	local eval eval INPUT=$RESULT_DIR'$'TOOLS_${TOOL}_DIR/$LIBRARY/'$'TOOLS_${TOOL}_SINGLE_LOCATED
    fi
    fileFoundSet "$INPUT" "$2"
}

if [ ! -z "$SPLICES" ]
then
    if eval [ -z '$'SPLICE_${SPLICES}_DIR ]
	then
	echo "Invalid tool $SPLICES" >&2
	exit 3
    fi

    
    check_genome_refseq SPLICE CURR_GENOME
    check_library '$'SPLICE_${SPLICES}_DIR CURR_LIBRARY
    check_threshold
    check_refseq SPLICE "$REFSEQ_FILE" REFSEQ_FILE
    check_outdir $RESULT_DIR'$'SPLICE_${SPLICES}_OUT_DIR/$LIBRARY/$LIB_NAME

    INPUT=`eval eval echo $CURR_LIBRARY/$LIB_NAME/'$'SPLICE_${SPLICES}_INPUT_FILES`
    fileFound  $INPUT
    if [ $? -eq 0 ]
	then
	INPUT=`eval eval echo $CURR_LIBRARY/'$'SPLICE_${SPLICES}_INPUT_FILES`
    fi

    fileFound $INPUT
    if [ $? -eq 0 ]
	then
	echo "I cannot find the input files. Please check the library name (maybe you should specify the --library-dir option).
  Didn't find $INPUT " >&2
	exit 5
    fi
    
    verbose "Launching the script $SPLICE_SCRIPT"

    for i in $SPLICE_JUNCTION_UPPER_BOUNDS
      do
      eval OUTPUT=$RESULT_DIR'$'SPLICE_${SPLICES}_OUT_DIR/$LIBRARY/$LIB_NAME/$SPLICE_SCRIPT_OUTPUT$i
      verbose "\t$SPLICE_SCRIPT $REFSEQ_FILE $THRESHOLD $i $INPUT > $OUTPUT"
      eval $SPLICE_SCRIPT $REFSEQ_FILE $THRESHOLD $i $INPUT > $OUTPUT
    done
fi

if [ ! -z "$SNP" -o "$INDEL" == "CRAC" ]
then
    if [ "$INDEL" == "CRAC" ]; then
	SNP=$INDEL
	INDEL=1
    fi
    
    check_genome
    check_tool "$SNP"
    check_library '$'TOOLS_${SNP}_DIR 
    
    # Check if the results exist for this tool and record them in RESULT_FILE
    get_results "$SNP" RESULT_FILE

    FINAL_RESULT_FILE=$RESULT_FILE.vcf
    if [ "$SNP" == "Bowtie" -o "$SNP" == "SOAP2" ]
	then
	FINAL_RESULT_FILE=$RESULT_FILE.sam.vcf
    fi

    # We process the SNPs if the refseq file was not given or if the SNP file does not
    # exist yet
    if [ "$SNP" != "CRAC" ] && [ -z "$REFSEQ_FILE" -o ! -s "$FINAL_RESULT_FILE" ]; then

	CONVERTED=0
    # if Bowtie, convert the output
	if [ "$SNP" == "Bowtie" -o "$SNP" == "SOAP2" ]
	    then
	    verbose "Convert $SNP output to SAM"
	    eval SCRIPT='$'$SNP'_TO_SAM'
	    $SCRIPT "$RESULT_FILE" > "$RESULT_FILE.sam"
	    RESULT_FILE="$RESULT_FILE.sam"
	    CONVERTED=1
	fi

    # Check if the genome .fa exists
	eval GENOME_FA='$'TOOLS_${SNP}_GENOME_DIR/$GENOME.fa
	fileFound "$GENOME_FA"
	if [ $? -eq 0 ]; then
	    GENOME_FA="$TOOLS_SOAP2_GENOME_DIR/$GENOME.fa" 
	fi
	if [ ! -f "$GENOME_FA" -o ! -r "$GENOME_FA" ]
	    then
	    echo "I cannot find the FASTA file of the genome in $GENOME_FA" >&2
	    exit 20
	fi
	
    # Should we create the fai file?
	if [ ! -f "${GENOME_FA}.fai" ]
	    then
	# Create it
	    verbose "Creating index: $SNP_CREATE_INDEX $GENOME_FA"
	    $SNP_CREATE_INDEX $GENOME_FA
	fi

    # Are the output consistent on the chromosome name?
	if [[ $(cat $RESULT_FILE | head -1 | awk '{print $4}') =~ ^chr ]] && [[ ! $(head -1 "${GENOME_FA}.fai" | awk '{print $1}') =~ ^chr ]]
	    then
	    verbose "Removing chr in chromosome names"
	    sed -i 's/chr//' "$RESULT_FILE"
	fi

    # Creating BAM and sorting it (use the variable needed by SNP_CREATE_SORTED_BAM)
	COMMAND=${SNP_CREATE_SORTED_BAM//@INDEX_FAI/${GENOME_FA}.fai}
	COMMAND=${COMMAND//@SAM_FILE/$RESULT_FILE}
	COMMAND=${COMMAND//@SAM_SORTED_BAM/$RESULT_FILE$SNP_SORTED_BAM_EXT}
	verbose 'Creating BAM and sorting it: '$COMMAND
	eval $COMMAND

	COMMAND=${SNP_CREATE_PILEUP//@INDEX_FA/$GENOME_FA}
	COMMAND=${COMMAND//@SAM_SORTED_BAM/$RESULT_FILE$SNP_SORTED_BAM_EXT}
	COMMAND=${COMMAND//@SAM_FILE/$RESULT_FILE}
	verbose 'Creating pileup file: '$COMMAND
	eval $COMMAND

#    verbose "Delete sorted BAM file $RESULT_FILE$SNP_SORTED_BAM_EXT.bam"
#    rm -f $RESULT_FILE$SNP_SORTED_BAM_EXT.bam
#	verbose "Delete BCF file $RESULT_FILE.bcf"
#	rm -f $RESULT_FILE.bcf

    # If we converted the output, we delete it
	if [ $CONVERTED -eq 1 ]
	    then
	    verbose "Delete converted SAM output $RESULT_FILE"
	    rm -f $"RESULT_FILE"
	fi
    fi

    # If a RefSeq is given, we try to cross the SNP found
    if [ ! -z "$REFSEQ_FILE" ]; then
	check_genome_refseq SNP REFSEQ_DIR
	check_refseq SNP "$REFSEQ_FILE" REFSEQ_FILE

	if [  -z "$THRESHOLD_ACCEPT" ]; then
	    THRESHOLD_ACCEPT=$DEFAULT_THRESHOLD_SNP
	fi

	# If the tool is CRAC we first need to preprocess the SNP input file
	if [ "$SNP" == "CRAC" ]; then
	    CRAC_SNP_FILE=''
	    if [ "$INDEL" == "1" ]; then
		for mutation in $FLUX_CAUSES_CRAC_indel; do
		    echo $mutation
		    eval TMP=$RESULT_DIR/$TOOLS_CRAC_DIR/$LIBRARY/$TOOLS_CRAC_BASENAME$mutation
		    CRAC_SNP_FILE=$TMP' '$(echo $CRAC_SNP_FILE)
		done
	    else
		for mutation in $FLUX_CAUSES_CRAC_snp; do
		    eval TMP=$RESULT_DIR/$TOOLS_CRAC_DIR/$LIBRARY/$TOOLS_CRAC_BASENAME$mutation
		    CRAC_SNP_FILE=$TMP' '$(echo $CRAC_SNP_FILE)
		done
	    fi
	    FINAL_RESULT_FILES=''
	    for file in $CRAC_SNP_FILE; do
		if [ ! -f "$file" ]; then
		    verbose "SNP file for CRAC cannot be found: $file"
		    exit 29
		fi
		FINAL_RESULT_FILE=$file.raw
		FINAL_RESULT_FILES=$FINAL_RESULT_FILES" "$FINAL_RESULT_FILE
		if [ ! -f "$FINAL_RESULT_FILE" ]; then
		    verbose "Extracting SNPs from CRAC file $file"
		    $CRAC_SNP_EXTRACTION_SCRIPT $file | sort -k1,1 -k2n,2 | uniq > $FINAL_RESULT_FILE
		fi
	    done
	else
	    FINAL_RESULT_FILES=$FINAL_RESULT_FILE
	fi

	eval OUT_DIR=$RESULT_DIR/'$SNP_'$SNP'_OUT_DIR'/$LIBRARY

	if [ ! -d $"OUT_DIR" ] && ! mkdir -p "$OUT_DIR" ; then
	    echo "Unable to create directory $OUT_DIR. Abort." >&2
	    exit 30
	fi
	for FINAL_RESULT_FILE in $FINAL_RESULT_FILES; do
	    BASENAME=$(basename $FINAL_RESULT_FILE)
	    COMMAND="$CROSS_SNP_REFSEQ_SCRIPT -t $THRESHOLD_ACCEPT --false-positives $OUT_DIR/$BASENAME.t$THRESHOLD_ACCEPT.fp --true-positives $OUT_DIR/$BASENAME.t$THRESHOLD_ACCEPT.tp $FINAL_RESULT_FILE $REFSEQ_FILE"
	    verbose "Launching $COMMAND > $OUT_DIR/$BASENAME.t$THRESHOLD_ACCEPT.stats"
	    $COMMAND > $OUT_DIR/$BASENAME.t$THRESHOLD_ACCEPT.stats
	done
    fi
fi


if [ ! -z "$TOOLS" ]
then
    check_tool "$TOOLS"
    check_library '$'TOOLS_${TOOLS}_DIR CURR_LIBRARY
    check_CRAC_library

    eval OUTPUT_DIR=$RESULT_DIR'$'TOOLS_${TOOLS}_OUTPUT_DIR/$LIBRARY
    check_outdir $OUTPUT_DIR

    eval eval INPUT=$CURR_LIBRARY/'$'TOOLS_${TOOLS}_INPUT_FILES
    eval CRAC_INPUT=$RESULT_DIR$TOOLS_CRAC_DIR/$LIBRARY/$TOOLS_CRAC_INPUT_FILES
    eval OUT_FILE='$'TOOLS_${TOOLS}_OUT_FILE
    verbose "Launching the script $TOOLS_SCRIPT"
    verbose "\t-> $TOOLS_SCRIPT $INPUT $CRAC_INPUT > $OUT_FILE"
    cd $OUTPUT_DIR
    eval $TOOLS_SCRIPT $INPUT $CRAC_INPUT > $OUT_FILE
    cd -
fi


if [ ! -z "$SPLICES_CROSSED" ]
then
    check_library '$'TOOLS_${SPLICES_CROSSED}_OUTPUT_DIR 

    if [ ! -f $RESULT_DIR$SPLICE_CRAC_OUT_DIR/$LIBRARY/$SPLICE_SCRIPT_OUTPUT$SPLICES_CROSSED_JUNCTION ]
	then
	echo "No cross result ($SPLICE_SCRIPT_OUTPUT$SPLICES_CROSSED_JUNCTION) "
	" RefSeq/CRAC for $LIBRARY" >&2
	exit 12
    fi
    eval eval INPUT_FILE=$RESULT_DIR'$'TOOLS_${SPLICES_CROSSED}_OUTPUT_DIR/$LIBRARY/'$'SPLICES_CROSSED_${SPLICES_CROSSED}_FILES
    FIRST_FILE=`eval echo $INPUT_FILE`
    eval FIRST_FILE=`echo $FIRST_FILE | awk '{print $1}'`
    FIRST_FILE=`basename $FIRST_FILE`

    eval EXTRACT_PREFIX='$'SPLICES_CROSSED_${SPLICES_CROSSED}_EXTRACT_PREFIX
    eval PREFIX='$'{FIRST_FILE$EXTRACT_PREFIX}

    eval eval OUTPUT_FILE=$RESULT_DIR'$'TOOLS_${SPLICES_CROSSED}_OUTPUT_DIR/$LIBRARY/$'$'SPLICES_CROSSED_${SPLICES_CROSSED}_OUT_FILE
    
    verbose "Launching script $SPLICES_CROSSED_SCRIPT"

    verbose "\t-> $SPLICES_CROSSED_SCRIPT $RESULT_DIR$SPLICE_CRAC_OUT_DIR/$LIBRARY/$SPLICE_SCRIPT_OUTPUT$SPLICES_CROSSED_JUNCTION $INPUT_FILE > $OUTPUT_FILE"
    eval eval $SPLICES_CROSSED_SCRIPT $RESULT_DIR$SPLICE_CRAC_OUT_DIR/$LIBRARY/$SPLICE_SCRIPT_OUTPUT$SPLICES_CROSSED_JUNCTION $INPUT_FILE  > $OUTPUT_FILE
fi

if [ ! -z "$CHIMERA" ]
then
    check_CRAC_library CRAC_DIR
    check_threshold

    fileFoundSet $CRAC_DIR/$CHIMERA_INPUT_FILES INPUT_CHIMERA_FILE
    check_outdir $RESULT_DIR$CHIMERA_OUTPUT_DIR/$LIBRARY/
    
    PREFIX=${INPUT_CHIMERA_FILE##*/}
    PREFIX=${PREFIX%.*}
    eval OUTPUT_FILE=$RESULT_DIR$CHIMERA_OUTPUT_DIR/$LIBRARY/$CHIMERA_OUTPUT_FILE
    eval OUTPUT_FILE_REFSEQ=$RESULT_DIR$CHIMERA_OUTPUT_DIR/$LIBRARY/$CHIMERA_OUTPUT_FILE_REFSEQ

    COMMAND="$SCRIPT_DIR/$CHIMERA_SCRIPT '$INPUT_CHIMERA_FILE' $THRESHOLD --output-reads $OUTPUT_FILE $FLAG > $OUTPUT_FILE.tab"
    verbose "Launching $COMMAND"
    eval $COMMAND

    if [ ! -z "$REFSEQ_FILE" ]; then
	check_genome_refseq SPLICE CURR_REFSEQ 
    
	CROSS_REFSEQ_FILE=$RESULT_DIR$SPLICE_CRAC_OUT_DIR/$LIBRARY/$CHIMERA_REFSEQ_OUTPUT_FILE$CHIMERA_NB_BASES_TOLERANCE
	OUT_DIR=`dirname $CROSS_REFSEQ_FILE`
	if [ ! -d "$OUT_DIR" ]; then
	    mkdir -p "$OUT_DIR"
	fi
	check_refseq SPLICE '' REFSEQ_FILE
	verbose "$SPLICE_CRAC_OUT_DIR/$CHIMERA_REFSEQ_OUTPUT_FILE$CHIMERA_NB_BASES_TOLERANCE does not exist. Computing it"
	verbose "Launching $SPLICE_SCRIPT $REFSEQ_FILE $THRESHOLD $CHIMERA_NB_BASES_TOLERANCE $OUTPUT_FILE > $CROSS_REFSEQ_FILE"
	$SPLICE_SCRIPT $REFSEQ_FILE $THRESHOLD $CHIMERA_NB_BASES_TOLERANCE $OUTPUT_FILE > $CROSS_REFSEQ_FILE
	TMPFILE=$RESULT_DIR"/tmp$RANDOM$RANDOM"
	verbose "\t-> sed -n '1,/transcripts/p' $CROSS_REFSEQ_FILE | grep -E \"$CHIMERA_REFSEQ_GOOD_CATEGORIES\" | sed -r 's/^(\S+).*?\s+([0-9]+)$/\\1 \\2/' > $TMPFILE"
	sed -n '1,/transcripts/p' $CROSS_REFSEQ_FILE | grep -E "$CHIMERA_REFSEQ_GOOD_CATEGORIES" | sed -r 's/^(\S+).*?\s+([0-9]+)$/\1 \2/' > $TMPFILE
	verbose " Retrieving CRAC results with good IDs selected from the RefSeq"
	
	awk 'NR==FNR {_[$2]=$1;next} $1 in _{print}' $TMPFILE $OUTPUT_FILE > $OUTPUT_FILE_REFSEQ
	rm -f $TMPFILE $TMPFILE2
    
	verbose "Final results in $OUTPUT_FILE_REFSEQ"
    else
	verbose "Final results in $OUTPUT_FILE"
    fi
fi

if [ ! -z "$CAUSES" -o ! -z "$FLUX_TAGS" -o ! -z "$FLUX_TAGS_FP" ]
then
    check_genome
    check_libdir DIR

    LIB=''
    if [ -z $LIBRARY ]
	then
	echo "A library name (-l) should be given."
	exit 22
    fi
    if [ ! -z $LIB_NAME ]
	then
	LIB=$LIB_NAME
    else
	LIB=$LIBRARY
    fi

    for ext in $FLUX_ERRORS_EXTENSION $FLUX_MUTATION_EXTENSION $FLUX_NO_MUTATION_BED_EXTENSION
    do
      fileFound "$DIR/$LIB*$ext"
      if [ $? -eq 0 ]
	  then
	  echo "No $ext file for the $LIB library. Aborting." >&2
	  exit 23
      fi
    done

    TOOL_NAME=${CAUSES:-$FLUX_TAGS}
    TOOL_NAME=${TOOL_NAME:-$FLUX_TAGS_FP}
    OUTPUT_DIR=`eval echo $RESULT_DIR/$FLUX_CROSS_OUTPUT_DIR/$LIBRARY`
    check_outdir $OUTPUT_DIR

    THRESHOLD_EXT=''
    if [ ! -z "$THRESHOLD" ]; then
	THRESHOLD_EXT='t'$THRESHOLD'.'
    fi

    if [ ! -z "$CAUSES" ]
	then
	OUTPUT_NOT_FOUND=$FLUX_CAUSES_NOT_FOUND_OUTPUT
	OUTPUT_STATS=$FLUX_CAUSES_STAT_OUTPUT
	SCRIPT=$FLUX_CAUSES_SCRIPT
    else
	OUTPUT_NOT_FOUND=$FLUX_TAGS_NOT_FOUND_OUTPUT
	OUTPUT_STATS=$FLUX_TAGS_STAT_OUTPUT	
	SCRIPT=$FLUX_TAGS_SCRIPT
	if [ ! -z "$FLUX_TAGS_FP" ]; then
	    OUTPUT_STATS=$FLUX_TAGS_FP_STAT_OUTPUT
	fi
    fi

    OPTIONS=''
    if [ ! -z "$CAUSES" ]
	then
#	OPTIONS=" --chr-lengths $GENOME_DIR/SSA/$GENOME.conf"
	OPTIONS=" --chr-lengths $INDEX_DIR/crac/$GENOME.conf"
	if [ "$CAUSES" == "CRAC" -o "$CAUSES"=="TopHat" -o "$CAUSES"=="TopHat2" -o "$CAUSES" == "MapSplice" -o "$CAUSES" == "TopHatFusion" -o "$CAUSES" == "TopHatFusionPost" ]
	then
	    SUFFIX_OUTPUT=''
	    if [ ! -z "$THRESHOLD" ]; then
		SUFFIX_OUTPUT=.t$THRESHOLD
	    fi
	    if [ ! -z $THRESHOLD_ACCEPT ]
		then
		OPTIONS=$OPTIONS" --threshold-splice $THRESHOLD_ACCEPT"
		SUFFIX_OUTPUT=$SUFFIX_OUTPUT".threshold-$THRESHOLD_ACCEPT"
	    fi
	fi    
    fi

    OPTIONS=$OPTIONS" --bed $DIR/$LIB*.$FLUX_NO_MUTATION_BED_EXTENSION \
--flux-err $DIR/$LIB*.$FLUX_ERRORS_EXTENSION \
--mutations $DIR/$LIB*.$FLUX_MUTATION_EXTENSION"
    if [ -z "$FLUX_TAGS_FP" ]; then
	OPTIONS=$OPTIONS" --output-missing $OUTPUT_DIR/Cross.$LIB.$OUTPUT_NOT_FOUND$NO_DUP$SUFFIX_OUTPUT"
    fi

    if [ ! -z "$CAUSES" ] 
	then
	OPTIONS=$OPTIONS" --output-false $OUTPUT_DIR/Cross.$LIB.$FLUX_CAUSES_FP_OUTPUT$NO_DUP$SUFFIX_OUTPUT"
	if [ "$CAUSES" == "CRAC" ] 
	    then
	    OPTIONS=$OPTIONS" --output-tag-missing $OUTPUT_DIR/Cross.$LIB.$FLUX_CAUSES_TAG_MISSING_OUTPUT$NO_DUP$SUFFIX_OUTPUT"
	    OPTIONS=$OPTIONS" --output-bad-classification $OUTPUT_DIR/Cross.$LIB.$FLUX_CAUSES_BAD_CLASSIFICATION_OUTPUT$NO_DUP$SUFFIX_OUTPUT"
	fi	    
    fi


    if eval [ \"'$FLUX_CAUSES_'$CAUSES\" == 1 ] 
	then
	TOOL_NAME=${CAUSES:-$FLUX_TAGS}
	in_dir=''
	if eval [ ! -d '$RESULT_DIR$TOOLS_'$TOOL_NAME'_DIR/$LIBRARY' ]
	    then
	    echo "$TOOL_NAME has not been launched on $LIBRARY" >&2
	    exit 24
	else
	    eval in_dir='$RESULT_DIR$TOOLS_'$TOOL_NAME'_DIR/$LIBRARY'
	fi

	# Search a sub-dir which have the full lib name
	if [ -d $in_dir/$LIB_NAME ]
	    then
	    in_dir=$in_dir/$LIB_NAME
	fi
	for causes in $FLUX_CAUSES 
	do
	  eval ELEM='$'FLUX_CAUSES_${TOOL_NAME}'_'$causes
	  if [ ! -z "$ELEM" ] 
	      then
	      for elem in $ELEM
		do
		FILENAME=`eval eval echo '$in_dir/$TOOLS_'$TOOL_NAME'_BASENAME$elem'`
		echo $FILENAME
		fileFound $FILENAME
		if [ $? -ne 0 ]
		    then
		    
		    OPTIONS=$OPTIONS" --$causes $FILENAME"
		fi
	      done
	  fi
	done	
        if [ "$CAUSES" == "CRAC" ] && [ ! -z "$NO_DUP" ]; then
	    OPTIONS=$OPTIONS" --single"
        fi
#     elif [ "$CAUSES" == "TopHat" -o "$CAUSES" == "MapSplice" ]
# 	then
# 	current_tool=$CAUSES
# 	tophat_dir=
# 	for file in eval echo  '$RESULT_DIR$SPLICE_'$CAUSES'_DIR/$LIBRARY*/'
# 	  do
# 	  if [ -d $file ]
# 	      then
# 	      tophat_dir=$file
# 	      break
# 	  fi
# 	done
# 	if [ -z $tophat_dir ] 
# 	    then
# 	    echo "No $CAUSES result for $LIBRARY library" >&2
# 	    exit 26
# 	fi
# 	INPUT=`eval echo '$tophat_dir/$LIB_NAME/$SPLICE_'$CAUSES'_INPUT_FILES'`
# 	fileFound  $INPUT
# 	if [ $? -eq 0 ]
# 	    then
# 	    echo "$INPUT not found"
# 	    INPUT=`eval echo '$RESULT_DIR$SPLICE_'$CAUSES'_DIR/$LIBRARY/$SPLICE_'$CAUSES'_INPUT_FILES'`
# 	fi
# 	fileFound $INPUT
# 	if [ $? -eq 0 ]
# 	    then
# 	    echo "I cannot find the input files. Please check the library name.
#   Didn't find $INPUT " >&2
# 	    exit 27
# 	fi
# 	OPTIONS=$OPTIONS" --splicing $INPUT"
    else
	if eval [ ! -d $RESULT_DIR'$'TOOLS_${TOOL_NAME}_DIR/$LIBRARY ]
	    then
	    eval echo "Invalid library $LIBRARY: $RESULT_DIR"'$'TOOLS_${TOOL_NAME}_DIR"/$LIBRARY does not exist"  >&2
	    exit 28
	fi
	if [ ! -z "$FLUX_TAGS_FP" ]; then
	    # Check that --accuracy has been launched on that dataset.
	    fileFound $OUTPUT_DIR/Cross.$LIB$ACCURACY_STATS_FALSE_POS
	    if [ $? -eq 0 ]; then
		echo "I cannot find the false positives generated with --accuracy option. Please first launch the script with this option" >&2
		exit 29
	    fi
	    OPTIONS=$OPTIONS" --tool-name SensitivityFlux --tools-all "$OUTPUT_DIR/Cross.$LIB$ACCURACY_STATS_FALSE_POS
	else
	    # FLUX_TAGS normal
	    OPTIONS=$OPTIONS" --tool-name $TOOL_NAME "
	    if [ $TOOL_NAME == "CRAC" ]
		then
		eval INPUT=$RESULT_DIR$TOOLS_CRAC_DIR/$LIBRARY/$TOOLS_CRAC_BASENAME
		OPTIONS=$OPTIONS" --tools-multiple ${INPUT}multiple --tools-multiple ${INPUT}duplication --tools-single ${INPUT}single"
	    else
		INPUT=`eval eval echo $RESULT_DIR'$'TOOLS_${TOOL_NAME}_DIR/$LIBRARY/'$'TOOLS_${TOOL_NAME}_INPUT_FILES`
		for file in $INPUT; do
		    OPTIONS=$OPTIONS" --tools-all $file"
		done
	    fi
	fi
    fi
    verbose "Launching $SCRIPT $OPTIONS > $OUTPUT_DIR/Cross.$LIB.$OUTPUT_STATS$NO_DUP$SUFFIX_OUTPUT"
    $SCRIPT $OPTIONS > $OUTPUT_DIR/Cross.$LIB.$OUTPUT_STATS$NO_DUP$SUFFIX_OUTPUT
fi

if [ ! -z "$ACCURACY" ]
then

    check_libdir DIR 
    check_library '$'TOOLS_${ACCURACY}_DIR TOOLS_RESULT_DIR
    

    if [ ! -z $LIB_NAME ]
	then
	LIB=$LIB_NAME
    else
	LIB=$LIBRARY
    fi

    get_results $ACCURACY INPUT_FILENAME
    # if eval [ ! -z '$'TOOLS_${ACCURACY}_SINGLE_LOCATED ]
    #     then
    #     fileFoundSet $TOOLS_RESULT_DIR/'$'TOOLS_${ACCURACY}_SINGLE_LOCATED INPUT_FILENAME
    # else
    #     fileFoundSet $TOOLS_RESULT_DIR/'$'TOOLS_${ACCURACY}_INPUT_FILES INPUT_FILENAME
    # fi

    TOOL_NAME=$ACCURACY
    eval OUTPUT_DIR=$RESULT_DIR/$FLUX_CROSS_OUTPUT_DIR/$LIBRARY
    check_outdir $OUTPUT_DIR

    OUTPUT_BASENAME=$OUTPUT_DIR/Cross.$LIB

    fileFoundSet $DIR'/'$LIB'*.'$FLUX_NO_MUTATION_BED_EXTENSION BED_FILE

    OPTIONS="--bed $BED_FILE  --tool-name $ACCURACY\
 --false-pos $OUTPUT_BASENAME$ACCURACY_STATS_FALSE_POS"

    if [ ! -z $ACCURACY_STATS_TRUE_POS ]
	then
	OPTIONS=$OPTIONS" --true-pos $OUTPUT_BASENAME$ACCURACY_STATS_TRUE_POS"
    fi
    OPTIONS=$OPTIONS" $INPUT_FILENAME"

    verbose "Launching script $ACCURACY_SCRIPT $OPTIONS"
    $ACCURACY_SCRIPT $OPTIONS > $OUTPUT_BASENAME$ACCURACY_STATS_EXT
fi

