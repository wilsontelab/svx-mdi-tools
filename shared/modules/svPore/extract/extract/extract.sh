# action:
#     collapse duplex reads into a single output read
#     create a coverage map
#     assemble the node path of every SV molecule
# input:
#     minimap2 PAF file
# output:
#     $EXTRACT_PREFIX.windowCoverage.txt.bgz[.tbi]
#     $EXTRACT_PREFIX.edges.no_sv.txt.gz
#     $EXTRACT_PREFIX.edges.sv.txt.gz
#     $EXTRACT_PREFIX.edges.tmp.txt.gz # ~10K reads for adapter training only

# log file feedback
echo "extracting SV information from aligned reads"
echo "input PAF: $NAME_PAF_FILE"
EDGE_GLOB=$EXTRACT_PREFIX.extend_edges.*.txt.gz # TODO: use /dev/shm for this buffer?
rm -f $EDGE_GLOB

# assemble the node path of every molecule to yield 
# one line per 2-node edge, one or more edges per molecule, in molecule order
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $EXTRACT_STEP_DIR/extract_nodes.pl | 
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

# finalize a low-resolution fragment coverage map over all aggregated molecules
# prints non-sv edges to file and sv and training edges to stream
perl $EXTRACT_STEP_DIR/window_coverage.pl |

# add information used for adapter finding, SV analysis, etc.
perl $EXTRACT_STEP_DIR/extend_edges.pl
checkPipe

# print edges to various files in preparation for analyze_edges.R
echo "splitting edge groups for downstream analysis"
zcat $EDGE_GLOB | 
perl $EXTRACT_STEP_DIR/split_edge_groups.pl 
checkPipe
rm -f $EDGE_GLOB

# index coverage map
echo "indexing window coverage map"
tabix -s 1 -b 2 -e 2 $COVERAGE_FILE
checkPipe

echo "done"
