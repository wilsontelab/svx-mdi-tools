# action:
#     group read pairs to unique molecule sequences, merge additional read pairs
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $NAME_BAM_FILE (possibly overridden by $BAM_FILE)
# outputs:
#     $COLLATE_PREFIX.amplicons.txt
#     $COLLATE_PREFIX.fastq.gz, for subsequent realignment to genome

# log file feedback
echo "grouping read alignments into amplicons"
source $GENOMEX_MODULES_DIR/source/check_name_bam_file.sh
echo "input bam: $NAME_BAM_FILE"

# get ready
INTERIM_FILE=$DATA_FILE_PREFIX.interim.txt.gz
export DISCOVERY_FILE=$DATA_FILE_PREFIX.discovery.txt
export ALLOWED_FILE=$DATA_FILE_PREFIX.allowed.txt
export AMPLICONS_FILE=$COLLATE_PREFIX.amplicons.txt
FASTQ_FILE=$COLLATE_PREFIX.fastq.gz
SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
SORT_RAM=$SORT_RAM"b"
SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# parse alignments into productive pairs
echo "parsing alignments into paired outer endpoint nodes"
slurp -s 10M $NAME_BAM_FILE |
samtools view -F 4 - | # strip SAM header; suppress unmapped reads (orphan reads persist for now)
perl $ACTION_DIR/group/parse_bam.pl | # extract signatures of all productively aligned molecules
pigz -p $N_CPU -c | 
slurp -s 10M -o $INTERIM_FILE
checkPipe

# determine if the input used UMIs to set the proper molecule counting method
echo "discovering UMI utilization status"
export N_UMI_PAIRS=`zcat $INTERIM_FILE | head -n 1000 | cut -f3 | sort | uniq | wc -l`
export GROUP_COUNT="count_distinct"
if [ "$N_UMI_PAIRS" == "1" ]; then export GROUP_COUNT="count"; fi

# count molecules for all pairs of outer endpoint nodes in the amplicon pool
echo "counting molecules for each pair of outer endpoint nodes"
zcat $INTERIM_FILE | 
cut -f1-3 | 
$SORT -k1,1 -k2,2 | 
bedtools groupby -g 1,2 -c 3 -o $GROUP_COUNT | 
$SORT -k3,3nr | 
sed -e 's/,/\t/g' | 
slurp -s 10M -o $DISCOVERY_FILE
checkPipe

# iteratively discover amplicon endpoint peaks
# apply fuzzy position logic to associate molecules with the most abundant index node pair(s)
echo "discovering amplicon outer endpoint index nodes"
Rscript $ACTION_DIR/group/discover-amplicons.R
checkPipe
cat $AMPLICONS_FILE

# reject molecules with too-large a position difference from and index and too low a count to become a new index
echo "counting and saving unique DNA sequences from kept amplicon molecules"
zcat $INTERIM_FILE | 
sed -e 's/,/\t/g' | 
perl $ACTION_DIR/group/filter-read-pairs.pl | 

# identify and count the unique sequences within the pool of kept ampliconic read pairs
# parse the corresponding QUAL for passing to realignment
# at this stage, N bases and sequencing errors still break molecule groups (will group later based on junctions)
$SORT -k1,1n -k2,2 -k3,3 | 
bedtools groupby -g 1,2,3 -c 4,5,6,7,7 -o first,first,first,first,count | # TODO: use better aggregate QUAL, e.g., via $val = ord ($char) - 33 ??
perl $ACTION_DIR/group/assemble-fastq.pl | 
pigz -p $N_CPU -c | 
slurp -o $FASTQ_FILE
checkPipe

echo "done"
