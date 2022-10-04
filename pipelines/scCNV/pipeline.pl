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

# example:
# $ENV{XXX_FILE} = "$ENV{DATA_FILE_PREFIX}.xxx";

1;
