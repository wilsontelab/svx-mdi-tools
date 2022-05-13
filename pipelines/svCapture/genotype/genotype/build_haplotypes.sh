# action:
#     find heterozygous and homozygous SNVs and other variants in the bulk
#     capture targets, i.e., in the source individual's constitutive genome,
#     using reads from all samples in the project, i.e., under $TASK_DIR
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $CONSTITUTIVE_VCF
# outputs:
#     

# log file feedback
echo "assembling genotype calls into constitutive base values per target position"
Rscript $ACTION_DIR/genotype/build_haplotypes.R
checkPipe

###############
exit 1

echo "done"
