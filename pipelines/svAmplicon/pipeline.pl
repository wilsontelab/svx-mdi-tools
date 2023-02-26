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

# force paired-end reads to be aligned as two interleaved single reads
$ENV{SUPPRESS_SMART_PAIRING} = 1; 

1;
