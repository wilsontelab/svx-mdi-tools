# action:
#     summarize information on final amplicon alignments
# optional:
#     $USE_CRAM         [default: creates .bam file]
# input:
#     realigment bam/cram file
# output:
#     $COLLATE_PREFIX.molecule_summary.txt.gz, for use in the app
#     data package

# log file feedback
echo "assembling molecule summary files"
echo "input bam: $NAME_REALIGNED_BAM_FILE"

# get ready
MOLECULE_FILE=$COLLATE_PREFIX.molecule_summary.txt.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

#########################
# export N_CPU=1

# parse alignments into productive pairs
echo "parsing alignments into paired outer endpoint nodes"
slurp -s 10M $NAME_REALIGNED_BAM_FILE |
samtools view - |
perl $ACTION_DIR/extract/extract_nodes.pl | 

$SORT -k8,8nr | 
# awk '$8 ~ /chr/' | 
# awk '$7 > 0 && $9 > 0' | 
cat
# head -n 100
exit 1

$SORT -k4,4nr | 
cat > /nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/lukaszewicz/svAmplicon/TMP.txt
exit 1

# awk '$5 == 1 && $7 == "P"' | #  && $8 > 20
# sed 's/~~~/\n/g' | 
# grep "N:" | 



# head -n 500
cat # > /dev/null

exit 1

# pigz -p $N_CPU -c | 
# slurp -s 10M -o $MOLECULE_FILE
# checkPipe

echo "done"
