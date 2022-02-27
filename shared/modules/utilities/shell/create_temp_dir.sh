# action:
#     create a temporary directory for large data files
# sets:
#     $TMP_DIR_WRK
#     $TMP_FILE_PREFIX
# creates:
#     $TMP_DIR_WRK

export TMP_DIR_WRK=$TMP_DIR_LARGE/$SUITE_NAME.$PIPELINE_NAME.$PIPELINE_ACTION/$DATA_NAME
export TMP_FILE_PREFIX=$TMP_DIR_WRK/$DATA_NAME
mkdir -p $TMP_DIR_WRK
trap "rm -rf $TMP_DIR_WRK" EXIT # makes sure we clean up the tmp dir however script exits
