# action:
#     summarize information on final amplicon alignments
# optional:
#     $USE_CRAM         [default: creates .bam file]
# input:
#     realigment bam/cram file
# output:
#     $EXTRACT_PREFIX.moleculeTypes.rds # for use in the app
#     $EXTRACT_PREFIX.junctions.rds
#     data package

# log file feedback
echo "assembling molecule summary files"
echo "input bam: $NAME_REALIGNED_BAM_FILE"

# get ready
export AMPLICONS_FILE=$EXTRACT_PREFIX.amplicons.txt
export INTERIM_FILE=$DATA_FILE_PREFIX.extract.interim.txt.gz
export MOLECULE_TYPES_FILE=$EXTRACT_PREFIX.moleculeTypes.rds
export JUNCTIONS_FILE=$EXTRACT_PREFIX.junctions.rds
export MANIFEST_FILE=$EXTRACT_PREFIX.manifest.csv
MOLECULE_FILE=$EXTRACT_PREFIX.molecule_summary.txt.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# parse alignments into amplicon alignments
echo "parsing alignments into nodes and edges"
slurp -s 10M $NAME_REALIGNED_BAM_FILE |
samtools view - |
perl $ACTION_DIR/extract/extract_nodes.pl | 

# break the node signature apart
sed 's/ZZ/\t/g' | 

# adjust the edge type when an insertion is present
#   any del + any ins = D 
#   no del  + any ins = I
awk 'BEGIN{OFS="\t"}{ 
    if($19 != "NA" && $19 > 0 && ($17 == "D" || $17 == "I")){ 
        $17 = $18 > 0 ? "D" : "I";
    }
    print $0;
}' | 

# assess microhomologies at inline (CIGAR) junctions
perl $ACTION_DIR/extract/adjust_inline_junctions.pl | 

# add base qualities to aligments for junction filtering
perl $ACTION_DIR/extract/add_base_quals.pl | 
pigz -p $N_CPU -c | 
slurp -s 10M -o $INTERIM_FILE
checkPipe

# apply sorting and grouping based on SVs within molecules
echo "sorting and grouping molecules based on the SVs they carry"
Rscript $ACTION_DIR/extract/sort_and_group.R
checkPipe

echo "cleaning up"
rm -f ${INTERIM_FILE}

echo "done"
