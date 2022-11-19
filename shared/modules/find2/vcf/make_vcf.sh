# action:
#     export svCalls to simplified VCF format for upload and other programs
# expects:
#     find/find_svs.sh
# outputs:
#     $FIND_PREFIX.structural_variants.vcf.bgz
#     $FIND_PREFIX.structural_variants.vcf.bgz.tbi

# get ready
export SAMPLES=`cat $LIBRARY_STAT_FILES | grep SAMPLE | sed 's/SAMPLE:\s//'`
export SAMPLES=`echo "$SAMPLES" | tr '\n' ' '`
OUT_FILE=$FIND_PREFIX.structural_variants.vcf.bgz

# write and compress the VCF file
echo "writing VCF file using breakends"
Rscript $ACTION_DIR/vcf/make_vcf.R | 
bgzip -c > $OUT_FILE
checkPipe

# index using tabix
echo "indexing VCF file using tabix"
rm -f $OUT_FILE.tbi
tabix -p vcf $OUT_FILE
checkPipe

echo "done"
