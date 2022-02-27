# action:
#     create a temporary directory in shared memory, i.e., SHM_DIR_WRK
# sets:
#     $SHM_DIR_WRK
#     $SHM_FILE_PREFIX
# creates:
#     $SHM_DIR_WRK

export SHM_DIR_WRK=/dev/shm/$SUITE_NAME.$PIPELINE_NAME.$PIPELINE_ACTION/$DATA_NAME
export SHM_FILE_PREFIX=$SHM_DIR_WRK/$DATA_NAME
mkdir -p $SHM_DIR_WRK
trap "rm -rf $SHM_DIR_WRK" EXIT # makes sure we clean up the tmp dir however script exits
