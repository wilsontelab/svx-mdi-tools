# action:
#     summarize information on final amplicon alignments
# optional:
#     $USE_CRAM         [default: creates .bam file]
# input:
#     realigment bam/cram file
# output:
#     $COLLATE_PREFIX.moleculeTypes.rds # for use in the app
#     $COLLATE_PREFIX.junctions.rds
#     data package

# log file feedback
echo "assembling molecule summary files"
echo "input bam: $NAME_REALIGNED_BAM_FILE"

# get ready
export INTERIM_FILE=$DATA_FILE_PREFIX.extract.interim.txt.gz
export MOLECULE_TYPES_FILE=$COLLATE_PREFIX.moleculeTypes.rds
export JUNCTIONS_FILE=$COLLATE_PREFIX.junctions.rds
export MANIFEST_FILE=$COLLATE_PREFIX.manifest.csv
MOLECULE_FILE=$COLLATE_PREFIX.molecule_summary.txt.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# parse alignments into productive pairs
echo "parsing alignments into paired outer endpoint nodes"
slurp -s 10M $NAME_REALIGNED_BAM_FILE |
samtools view - |
perl $ACTION_DIR/extract/extract_nodes.pl | 
pigz -p $N_CPU -c | 
slurp -s 10M -o $INTERIM_FILE
checkPipe

# apply sorting and grouping based on SVs within molecules
echo "sorting and grouping molecules based on the SVs they carry"
Rscript $ACTION_DIR/extract/sort_and_group.R
checkPipe

echo "done"
