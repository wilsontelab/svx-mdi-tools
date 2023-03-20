# action:
#     pull a list of inferred SV junctions
#     create a coverage map
# input:
#     minimap2 PAF file with CIGAR string
#     fastq files for sequence retrieval during junction analysis
# output:
#     $EXTRACT_PREFIX.nodes.txt.gz
#     $EXTRACT_PREFIX.windowCoverage.txt.bgz[.tbi]

# log file feedback
echo "assembling molecule summary files"
echo "input PAF: $NAME_PAF_FILE"

##########################
export N_CPU=1

# get ready
NODES_FILE=$EXTRACT_PREFIX.nodes.txt.gz
COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# assemble the node path of every molecule, one line per segment
echo "parsing alignments into molecule node paths"
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $ACTION_DIR/extract/extract_nodes.pl | 

head -n 200
exit 1

# adjust the segment type when an insertion is present
#   large del + any ins = D 
#   short/no del + large ins = I
awk 'BEGIN{OFS="\t"}{ 
    if($7 > 0 && ($4 == "D" || $4 == "I")){ $4 = $6 >= '$MIN_SV_SIZE' ? "D" : "I"}; 
    print $0;
}' | 

# TODO: print windowCoverage using the merged stream

pigz -p $N_CPU -c | 
slurp -s 10M -o $NODES_FILE
checkPipe

# # aggregate and index coverage map over all molecules
# echo "assembling window coverage map"
# zcat $EXTRACT_PREFIX.windowCoverage.*.txt.gz |
# perl $ACTION_DIR/extract/window_coverage.pl | 
# bgzip -c | 
# slurp -s 10M -o $COVERAGE_FILE
# checkPipe
# tabix -s 1 -b 2 -e 2 $COVERAGE_FILE
# checkPipe

echo "done"
