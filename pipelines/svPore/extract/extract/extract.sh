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

# get ready
export NODES_FILE=$EXTRACT_PREFIX.nodes.txt.gz
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
export SEQUENCES_FILE=$EXTRACT_PREFIX.sequences.txt.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

N_CPU_HOLD=$N_CPU
export N_CPU=1

# assemble the node path of every molecule, one line per segment
echo "parsing alignments into molecule node paths"
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $ACTION_DIR/extract/extract_nodes.pl | 

# adjust the segment type when an insertion is present
#   large del + any ins = D 
#   short/no del + large ins = I
# pad out any missing QSTART and QEND with NA
awk 'BEGIN{OFS="\t"}{ 
    if($7 > 0 && ($4 == "D" || $4 == "I")){ 
        $4 = $6 >= '$MIN_SV_SIZE' ? "D" : "I";
    }
    if(NF == 8){
        $10 = $8;
        $8 = "NA";
        $9 = "NA";
    }
    print $0;
}' | 

# finalize a fragment coverage map over all aggregated molecules
perl $ACTION_DIR/extract/window_coverage.pl | # repeats nodes to stream
pigz -p $N_CPU -c | 
slurp -s 10M -o $NODES_FILE
checkPipe

# index coverage map
echo "indexing window coverage map"
tabix -s 1 -b 2 -e 2 $COVERAGE_FILE
checkPipe

# extract read sequences for adapter discovery
export N_CPU=$N_CPU_HOLD
echo "extracting read sequences for adapter discovery"
zcat $NODES_FILE | 
perl $ACTION_DIR/extract/extract_reads.pl | 
pigz -p $N_CPU -c | 
slurp -s 10M -o $SEQUENCES_FILE
checkPipe

echo "done"
