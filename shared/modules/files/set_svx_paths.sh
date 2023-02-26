# action:
#     set file paths and prefixes for SV finding algorithsm
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/files/set_svx_paths.sh

# collate
export COLLATE_PREFIX=$DATA_GENOME_PREFIX.collate
export CONSENSUS_PREFIX=$DATA_GENOME_PREFIX.consensus

# extract
export EXTRACT_PREFIX=$DATA_GENOME_PREFIX.extract
export COMPILE_PREFIX=$DATA_GENOME_PREFIX.compile
export EXTRACT_GLOB_PREFIX=$TASK_DIR/*/*.$GENOME.extract # for multi-sample find
export COMPILE_GLOB_PREFIX=$TASK_DIR/*/*.$GENOME.compile

# coverage 
export COVERAGE_PREFIX=$DATA_GENOME_PREFIX.coverage

# find
export FIND_PREFIX=$DATA_GENOME_PREFIX.find
export FIND_GLOB_PREFIX=$TASK_DIR/*/*.$GENOME.find

# manifest
export MANIFEST_PREFIX=$DATA_GENOME_PREFIX.manifest

# genotype
export GENOTYPE_PREFIX=$DATA_GENOME_PREFIX.genotype
export CONSTITUTIVE_VCF=$GENOTYPE_PREFIX.constitutive.vcf.gz
export PADDED_TARGETS_BED=$DATA_FILE_PREFIX.padded_targets.bed
