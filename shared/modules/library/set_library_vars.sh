# action:
#     set environment variables for various paired-end library properties
# expects:
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
# usage:
#     source $MODULES_DIR/library/set_library_vars.sh

if [ ! -e $NAME_BAM_FILE ]; then
    echo "missing alignment file: $NAME_BAM_FILE"
    exit 1
fi
echo "alignment file: $NAME_BAM_FILE"

echo "extracting READ_LEN"
export READ_LEN=`samtools view -f 1 $NAME_BAM_FILE | head -n 1 | awk '{print length($10)}'`
echo "READ_LEN = $READ_LEN"

echo "extracting MAX_TLEN"
export MAX_TLEN=`samtools view $NAME_BAM_FILE | head -n 100000 | perl $MODULES_DIR/library/get_max_TLEN.pl`
echo "MAX_TLEN = $MAX_TLEN"
