# action:
#     assemble the node path of every target molecule from extract/rebasecall/realign
# input:
#     minimap2 PAF file with CIGAR string
#     fastq files for sequence retrieval during junction analysis
# output:
#     $REEXTRACT_PREFIX.edges.tmp.txt.gz # ~10K reads for adapter training only
#     $REEXTRACT_PREFIX.edges.sv.txt.gz

# log file feedback
echo "executing second-pass high-resolution analysis of target reads"
echo "input PAF: $NAME_PAF_FILE"
export EXTRACT_PASS=2
EDGE_GLOB=$EXTRACT_PREFIX.extend_edges.*.txt.gz

# assemble the node path of every molecule to yield 
# one line per 2-node edge, one or more edges per molecule, in molecule order
echo "parsing alignments into edges"
rm -f $EDGE_GLOB
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $PIPELINE_DIR/extract/extract_nodes.pl | 
sed 's/ZZ/\t/g' | 

# adjust the edge type when an insertion is present
#   any del + any ins = D 
#   no del + large ins = I
awk 'BEGIN{OFS="\t"}{ 
    if($11 != "NA" && $11 > 0 && ($9 == "D" || $9 == "I")){ 
        $9 = $10 > 0 ? "D" : "I";
    }
    print $0;
}' | 

# add information used for adapter finding, SV analysis, etc.
perl $PIPELINE_DIR/extract/extend_edges.pl 
checkPipe

# print edges to various files in preparation for analyze_edges.R
echo "splitting edge groups for downstream analysis"
zcat $EDGE_GLOB | 
perl $PIPELINE_DIR/extract/split_edge_groups.pl 
checkPipe

# rm -f $EDGE_GLOB

# exit 1

# # reformat and index read sequences for adapter discovery
# echo "indexing primary read data for adapter discovery and squiggle examination"
# samtools view $DATA_FILE_PREFIX.target_reads.unaligned.bam | 
# perl $PIPELINE_DIR/extract/index_reads.pl | 
# slurp -s 10M -o $SEQUENCES_FILE
# checkPipe

echo "done"
