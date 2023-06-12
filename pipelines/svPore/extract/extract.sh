# action:
#     assemble the low resolution node path of every molecule
#     create a coverage map
#     extract all SV and a representative subset of non-SV qNames for high accuracy analysis
# input:
#     minimap2 PAF file (CIGAR string not required)
# output:
#     $EXTRACT_PREFIX.windowCoverage.txt.bgz[.tbi]
#     $EXTRACT_PREFIX.edges.no_sv.txt.gz
#     $EXTRACT_PREFIX.target.qNames.txt (cannot be compressed)

# log file feedback
echo "executing first-pass low-resolution analysis of all reads"
echo "input PAF: $NAME_PAF_FILE"
export EXTRACT_PASS=1

# assemble the node path of every molecule to yield 
# one line per 2-node edge, one or more edges per molecule, in molecule order
slurp -s 10M gunzip -c $NAME_PAF_FILE |
perl $PIPELINE_DIR/extract/extract_nodes.pl | 
sed 's/ZZ/\t/g' | 

# finalize a low-resolution fragment coverage map over all aggregated molecules
perl $PIPELINE_DIR/extract/window_coverage.pl # also prints edges and/or qNames to various files
checkPipe

# index coverage map
echo "indexing window coverage map"
tabix -s 1 -b 2 -e 2 $COVERAGE_FILE
checkPipe

echo "done"
