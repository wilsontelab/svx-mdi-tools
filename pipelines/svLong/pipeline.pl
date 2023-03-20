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

# all task environment variables may be accessed here

# set options for paired read alignment
defined $ENV{WINDOW_POWER} and $ENV{WINDOW_SIZE} = 10 ** $ENV{WINDOW_POWER};
defined $ENV{SV_SIZE_POWER} and $ENV{MIN_SV_SIZE} = 10 ** $ENV{SV_SIZE_POWER};
$ENV{USE_CHR_M} = 1;

1;
