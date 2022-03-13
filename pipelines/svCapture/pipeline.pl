#!/usr/bin/perl
use strict;
use warnings;

# force the behavior of svx find algorithms
$ENV{IS_COLLATED} = 1;
$ENV{TARGET_SCALAR} = 10;

1;
