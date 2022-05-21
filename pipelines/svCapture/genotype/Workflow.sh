#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
source $MODULES_DIR/utilities/shell/create_shm_dir.sh
export MAX_SORT_RAM_INT=$(($TOTAL_RAM_INT - 4000000000))

# determine if this is a single-sample or a multi-sample analysis
# establish our list of input bam/cram files
echo "setting sample mode"
if [[ "$USE_CRAM" = "" || "$USE_CRAM" = "0" || "$USE_CRAM" = "null" ]]; then
    export NAME_BAM_SUFFIX=*.$GENOME.name.realigned.bam
else
    export NAME_BAM_SUFFIX=*.$GENOME.name.realigned.cram
fi
export SAMPLE_MODE=single
export NAME_BAM_FILES=`ls -1 $NAME_BAM_SUFFIX 2>/dev/null`
if [ "$NAME_BAM_FILES" = "" ]; then
    export SAMPLE_MODE=multiple
    export NAME_BAM_FILES=`ls -1 */$NAME_BAM_SUFFIX 2>/dev/null` 
    if [ "$NAME_BAM_FILES" = "" ]; then
        echo "could not find input bam files in:"
        echo "    $TASK_DIR"
        exit 1
    fi
fi
echo "  $SAMPLE_MODE"
export COORDINATE_BAM_FILES=`echo $NAME_BAM_FILES | sed 's/name\.realigned/coordinate\.realigned/g'`

# create a filterable BED file of padded target regions
awk 'BEGIN{OFS="\t"}{
    $2 -= '$REGION_PADDING';
    $3 += '$REGION_PADDING';
    print $1, $2, $3;
}' $TARGETS_BED > $PADDED_TARGETS_BED

# skip first pipeline actions if using any externally provided haplotype file
if [[ "$HAPLOTYPE_FILE" == "" || "$HAPLOTYPE_FILE" == "NA" ||  "$HAPLOTYPE_FILE" == "null" ]]; then

    # sort the input bam/cram file(s), in series
    runWorkflowStep 1 sort $GENOMEX_MODULES_DIR/align/sort_bam_file.sh # samtools sort

    # call and normalize constitutive variants common to all sample in padded target regions
    runWorkflowStep 2 genotype genotype/genotype_source.sh # bcftools mpileup/call/norm
    runWorkflowStep 3 build    genotype/build_haplotypes.sh 

    export HAPLOTYPE_FILE=$GENOTYPE_PREFIX.unphased_haplotypes.rds
fi

# find SV molecules with SNVs / indels that cannot be accounted for by the source alleles
runWorkflowStep 4 find genotype/find_SNVs.sh

# clean up
rm -fr $TMP_DIR_WRK
