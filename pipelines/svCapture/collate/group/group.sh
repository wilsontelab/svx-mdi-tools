# action:
#     group read pairs and make consensus per source molecule
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $NAME_BAM_FILE (possibly overridden by $BAM_FILE)
# outputs:
#     fastq-format file ready for realignment to genome

# log file feedback
echo "grouping read pairs to unique, consensus source molecules"
source $GENOMEX_MODULES_DIR/source/check_name_bam_file.sh
echo "input bam: $NAME_BAM_FILE"

# do the work, supporting one or more bam files in NAME_BAM_FILE
samtools cat $NAME_BAM_FILE |
samtools view -F 4 - | # strip SAM header; suppress unmapped reads (orphan reads persist)
perl $ACTION_DIR/group/parse_bam.pl | # extract signatures of all unique input source molecules
sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $MAX_SORT_RAM_INT"b" --compress-program=pigz -k1,1 -k2,2n | # sort by mol_key/pos1
perl $ACTION_DIR/group/group_reads.pl | # apply fuzzy position logic to establish molecule groups
perl $ACTION_DIR/consensus/make_consensus.pl # merge reads and build strand/molecule consensuses
checkPipe

echo "done"
