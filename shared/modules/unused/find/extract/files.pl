use strict;
use warnings;

# extract output file handling
# defines the files that collect the extraction output streams, per thread

# operating parameters
use vars qw($error $isTargeted $EXTRACT_PREFIX);

# file handle variables
our (
    $nodesH,    
    $spansH,
    $endpointsH
);

# open only the file handles needed by the extraction run
# operates per thread; files are concatenated later
sub openFileHandles {
    my ($childN) = @_;

    # candidate SV data nodes
    # one line per source molecule per node in the graph
    # used later by all extraction runs to find and call SVs
    my $file = getFileStream('nodes', $childN);
    open $nodesH, "|-", $$file{stream} or die $$file{error};

    # end-to-end spans of all molecules, i.e., both proper and SV
    # used to create a genome-wide read coverage map
    # only used if NOT targeted, i.e., if whole-genomic, to help find ~clonal SVs
    $file = getFileStream('spans', $childN);
    $isTargeted or open $spansH, "|-", $$file{stream} or die $$file{error};

    # outer endpoint nodes of all molecules, i.e., both proper and SV
    # used to set the SHARED_PROPER count in later steps
    # only used if targeted, i.e., if svCapture, where we expect large strand family sizes
    $file = getFileStream('endpoints', $childN);
    $isTargeted and open $endpointsH, "|-", $$file{stream} or die $$file{error};
}

# close all open file handles
sub closeFileHandles {
    $nodesH      and close $nodesH;    
    $spansH      and close $spansH;
    $endpointsH  and close $endpointsH;
}

# helper function 
sub getFileStream {
    my ($type, $childN) = @_;
    my $file = "$EXTRACT_PREFIX.$type.$childN.gz";
    {
        stream => "gzip -c | slurp -s 50M -o $file",
        error  => "$error: could not open $file: $!\n"
    }
}

1;
