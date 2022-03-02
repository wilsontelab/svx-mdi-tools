#!/usr/bin/perl
use strict;
use warnings;

# force name-sorted source file to the re-alignment bam
$ENV{BAM_FILE} = "$ENV{DATA_FILE_PREFIX}.$ENV{GENOME}.name.realigned.bam";

# force the behavior of svx find algorithms
$ENV{IS_SV_CAPTURE} = 1;
$ENV{LIBRARY_ORIENTATION} = "FF";

1;
