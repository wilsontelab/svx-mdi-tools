# action:
#     set environment variables for various paired-end library properties
# expects:
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/files/set_svx_paths.sh
# usage:
#     source $MODULES_DIR/library/set_library_vars.sh

# to support bam merging in the collate step
#     strips quotes from the list of name-sorted bam files
#     extract just the first file for metadata retrieval (subsequent files expected to be sufficiently similar)
if [[ "$BAM_FILE" != "" && "$BAM_FILE" != "null" ]]; then
    export NAME_BAM_FILE=$(echo "$BAM_FILE" | sed 's/"//g')
fi
FIRST_NAME_BAM_FILE="${NAME_BAM_FILE%% *}"

for BAM_FILE_CHECK in $NAME_BAM_FILE; do
    if [ ! -e $BAM_FILE_CHECK ]; then
        echo -e "bam file not found:\n    $BAM_FILE_CHECK"
        exit 1
    fi
done
echo "alignment file(s): $NAME_BAM_FILE"

echo "SAMPLE: $DATA_NAME" > $EXTRACT_PREFIX.library_stats.yml

echo "extracting READ_LEN"
# TODO: for unknown reasons, samtools view of cram files works fine for hg38 but hangs indefinitively for mm10
# it most likely relates to name sorting limitations describe here: https://www.htslib.org/workflow/cram.html
# it is unclear why this would differ between genomes
export READ_LEN=`
    samtools view $SAMTOOLS_VIEW_REFERENCE -f 1 $FIRST_NAME_BAM_FILE | 
    awk '$10 != "*"' | 
    head -n 10000 | 
    awk 'BEGIN{rl=0}{if(length($10) > rl) rl = length($10)}END{print rl}'
` # assumes smart-paired alignments
echo "READ_LEN = $READ_LEN"
echo "READ_LEN: $READ_LEN" >> $EXTRACT_PREFIX.library_stats.yml

echo "extracting MAX_TLEN"
export MAX_TLEN=`samtools view $SAMTOOLS_VIEW_REFERENCE $FIRST_NAME_BAM_FILE | head -n 100000 | perl $MODULES_DIR/library/get_max_TLEN.pl`
echo "MAX_TLEN = $MAX_TLEN"
echo "MAX_TLEN: $MAX_TLEN" >> $EXTRACT_PREFIX.library_stats.yml
