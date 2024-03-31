use strict;
use warnings;

# write alignments in the bam stream, one file per cell index
my ($fileCellIndex, $fileH) = (0);
while (my $line = <STDIN>){
    my ($cellIndex, $aln) = split("\t", $line, 2);
    if($cellIndex != $fileCellIndex){
        $fileH and close $fileH;
        $fileCellIndex = $cellIndex;
        my $file = "$ENV{ALIGNMENT_DIR}/$fileCellIndex.name.bam";
        open $fileH, "|-", "samtools view -b - > $file" or die "could not open: $file\n";
        print $fileH $ENV{BAM_HEADER}, "\n";
    }
    print $fileH $aln;
}
close $fileH;
