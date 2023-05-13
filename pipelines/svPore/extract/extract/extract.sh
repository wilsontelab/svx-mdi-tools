# action:
#     assemble the node path of every molecule
#     create a coverage map
# input:
#     minimap2 PAF file with CIGAR string
#     fastq files for sequence retrieval during junction analysis
# output:
#     $EXTRACT_PREFIX.windowCoverage.txt.bgz[.tbi]
#     $EXTRACT_PREFIX.edges.no_sv.txt.gz
#     $EXTRACT_PREFIX.edges.sv.txt.gz
#     $EXTRACT_PREFIX.sequences.txt.gz

# log file feedback
echo "assembling molecule summary files"
echo "input PAF: $NAME_PAF_FILE"

# assemble the node path of every molecule to yield 
# one line per 2-node edge, one or more edges per molecule, in molecule order
echo "parsing alignments into molecule node paths"
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $ACTION_DIR/extract/extract_nodes.pl | 

# adjust the edge type when an insertion is present
#   any del + any ins = D 
#   no del + large ins = I
awk 'BEGIN{OFS="\t"}{ 
    if($7 > 0 && ($4 == "D" || $4 == "I")){ 
        $4 = $6 > 0 ? "D" : "I";
    }
    print $0;
}' | 

# finalize a fragment coverage map over all aggregated molecules
sed 's/ZZ/\t/g' | 
perl $ACTION_DIR/extract/window_coverage.pl # also prints edges to various files
checkPipe

# index coverage map
echo "indexing window coverage map"
tabix -s 1 -b 2 -e 2 $COVERAGE_FILE
checkPipe

# extract read sequences for adapter discovery
echo "extracting and indexing read sequences for adapter discovery"
zcat $EDGES_TMP_FILE $EDGES_SV_FILE | 
perl $ACTION_DIR/extract/extract_reads.pl |  
slurp -s 10M -o $SEQUENCES_FILE
checkPipe

echo "done"
