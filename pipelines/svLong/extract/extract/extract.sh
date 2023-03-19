# action:
#     pull a list of inferred SV junctions
#     create a coverage map
# input:
#     minimap2 PAF file with CIGAR string
# output:
#     $EXTRACT_PREFIX.xxx
#     data package

# log file feedback
echo "assembling molecule summary files"
echo "input PAF: $NAME_PAF_FILE"

# get ready
# export INTERIM_FILE=$DATA_FILE_PREFIX.extract.interim.txt.gz
# export MOLECULE_TYPES_FILE=$COLLATE_PREFIX.moleculeTypes.rds
# export JUNCTIONS_FILE=$COLLATE_PREFIX.junctions.rds
# export MANIFEST_FILE=$COLLATE_PREFIX.manifest.csv
# MOLECULE_FILE=$COLLATE_PREFIX.molecule_summary.txt.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

export N_CPU=1

# NAME_PAF_FILE=/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svLong/test.paf.gz

# parse alignments into productive pairs
echo "parsing alignments into molecule node paths"
slurp -s 10M $NAME_PAF_FILE |
zcat |
perl $ACTION_DIR/extract/extract_nodes.pl | 

# sort -k1,1nr | 
head -n 100
exit 1

pigz -p $N_CPU -c | 
slurp -s 10M -o $INTERIM_FILE
checkPipe

# apply sorting and grouping based on SVs within molecules
echo "sorting and grouping molecules based on the SVs they carry"
Rscript $ACTION_DIR/extract/sort_and_group.R
checkPipe

echo "done"
