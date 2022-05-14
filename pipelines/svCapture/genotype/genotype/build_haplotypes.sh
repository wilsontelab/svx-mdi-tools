# action:
#     make a bi-allelic, unphased haplotype map of all padded target region bases
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $CONSTITUTIVE_VCF
# outputs:
#     $GENOTYPE_PREFIX.unphased_haplotypes.gz

# log file feedback
echo "assembling genotype calls into constitutive base values per target position"
Rscript $ACTION_DIR/genotype/build_haplotypes.R
checkPipe

echo "done"
