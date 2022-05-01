#!/usr/bin/perl
use strict;
use warnings;

# use pipeline.pl to set environment variables that are:
#   - derived from user options
#   - shared across multiple actions
# and/or to check the validity of existing environment variables as needed

# this script can be empty or deleted entirely if not needed

# command line option '--abc-def' becomes environment variable 'ABC_DEF'
# perl environment variable '$ENV{ABC_DEF}' becomes $ABC_DEF or ${ABC_DEF} in shell

# all upstream variables may be accessed here, e.g, PIPELINE_ACTION, PIPELINE_DIR, etc.

require "$ENV{MODULES_DIR}/genome/set_env_vars.pl";

$ENV{DATA_GENOME_PREFIX} = "$ENV{DATA_FILE_PREFIX}.$ENV{GENOME}";

$ENV{RATES_FILE} = "$ENV{DATA_GENOME_PREFIX}.rates.txt";
$ENV{COUNT_MATRIX_FILE} = "$ENV{DATA_GENOME_PREFIX}.bins.size_$ENV{BIN_SIZE}.bed.gz";   

1;
