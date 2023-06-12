#!/bin/bash

# process and set working paths
EXPANDED_INPUT_DIR=`echo ${INPUT_DIR}`
FASTQ_FILE=${DATA_FILE_PREFIX}.fastq.gz

# log report
echo "calling bases"
echo "  Dorado version: "${DORADO_VERSION}
echo "  model: "${ONT_MODEL}
echo "  input: "${EXPANDED_INPUT_DIR}
echo "  output: "${FASTQ_FILE}

# call bases from pod5 to fastq.gz
${DORADO_EXECUTABLE} basecaller ${ONT_MODEL_DIR} ${EXPANDED_INPUT_DIR} | # outputs unaligned SAM
samtools fastq - |
pigz -p $N_CPU -c | 
slurp -s 10M -o ${FASTQ_FILE}
checkPipe


# $ dorado --help
# Usage: dorado [options] subcommand

# Positional arguments:
# basecaller
# download
# duplex

# Optional arguments:
# -h --help               shows help message and exits
# -v --version            prints version information and exits
# -vv                     prints verbose version information and exits


# $ dorado basecaller --help
# Usage: dorado [-h] [--device VAR] [--read-ids VAR] [--max-reads VAR] [--min-qscore VAR] [--batchsize VAR] [--chunksize VAR] [--overlap VAR] [--num_runners VAR] [--modified-bases VAR...] [--modified-bases-models VAR] [--emit-fastq] [--emit-moves] model data

# Positional arguments:
#   model                         the basecaller model to run. 
#   data                          the data directory. 

# Optional arguments:
#   -h, --help                    shows help message and exits 
#   -v, --version                 prints version information and exits 
#   -v, --verbose          
#   -x, --device                  device string in format "cuda:0,...,N", "cuda:all", "metal" etc.. [default: "cuda:all"]
#   -l, --read-ids                A file with a newline-delimited list of reads to basecall. If not provided, all reads will be basecalled [default: ""]
#   -n, --max-reads               [default: 0]
#   --min-qscore                  [default: 0]
#   -b, --batchsize               if 0 an optimal batchsize will be selected [default: 0]
#   -c, --chunksize               [default: 10000]
#   -o, --overlap                 [default: 500]
#   -r, --num_runners             [default: 2]
#   --modified-bases              [nargs: 1 or more] 
#   --modified-bases-models       a comma separated list of modified base models [default: ""]
#   --emit-fastq           
#   --emit-moves

