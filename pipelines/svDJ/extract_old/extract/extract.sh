# action:
#     summarize information on final amplicon alignments
# input:
#     $NAME_BAM_FILE (possibly overridden by $BAM_FILE)
# output:

# log file feedback
echo "assembling molecule summary files"
source $GENOMEX_MODULES_DIR/source/check_name_bam_file.sh
echo "input bam: $NAME_BAM_FILE"

# get ready
export INTERIM_FILE=$DATA_FILE_PREFIX.extract.interim.txt.gz
export MOLECULE_TYPES_FILE=$EXTRACT_PREFIX.moleculeTypes.rds
export JUNCTIONS_FILE=$EXTRACT_PREFIX.junctions.rds
export MANIFEST_FILE=$EXTRACT_PREFIX.manifest.csv
MOLECULE_FILE=$EXTRACT_PREFIX.molecule_summary.txt.gz
# SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
# SORT_RAM=$SORT_RAM"b"
# SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# parse alignments into junction graph
echo "parsing alignments into junction nodes and edges"
slurp -s 10M $NAME_BAM_FILE |
samtools view - |
perl $ACTION_DIR/extract/extract_nodes.pl | 

# # adjust the segment type when an insertion is present
# #   large del + any ins = D 
# #   short/no del + large ins = I
# awk 'BEGIN{OFS="\t"}{ 
#     if($7 > 0 && ($4 == "D" || $4 == "I")){ $4 = $6 >= '$MIN_SV_SIZE' ? "D" : "I"}; 
#     print $0;
# }' | 

# break the node signature apart
sed 's/ZZ/\t/g' | 

# parse and add some additional information
perl $ACTION_DIR/extract/append_base_qual.pl | 
pigz -p $N_CPU -c | 
slurp -s 10M -o $INTERIM_FILE
checkPipe

# apply sorting and grouping based on SVs within molecules
echo "sorting and grouping molecules based on the junctions they carry"
Rscript $ACTION_DIR/extract/sort_and_group.R
checkPipe

echo "done"
