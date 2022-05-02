---
published: false
---

## wgaSeq: Whole Genome Amplification Sequencing

**wgaSeq** takes the output of low-depth, next-generation
sequencing as applied to single-cell, whole-genome
amplification (WGA). Reads are trimmed, merged, and aligned to
the reference genome and counts are made in fixed-width bins.

Given the expectations of relatively low-depth sequencing
combined with lower quality WGA input data (as compared to
standard bulk libraries), wgaSeq use a fast, low resolution 
model designed to support downstream detection of chromosome
aneuploidy and large segmental CNVs.

It is expected that a set of individual cells will have been
sequenced together in a single sequencing run and
thus subjected to common conditions and batch effects.
Different batches of cells should be analyzed in separate
calls to wgaSeq and compared later.

### Required Data Inputs

You must provide a set of fastq.gz format files, with one
pair of files for each cell corresponding to the two reads
of a paired-end sequencing run. These must be provided in one 
of the following two file structures:

```
--input-dir /path/to/project
```

```sh
# one cell per named sub-folder (the default, checked first)
/path/to/project/cell-1/*_R1_*.fastq.gz
/path/to/project/cell-1/*_R2_*.fastq.gz
/path/to/project/cell-2/*_R1_*.fastq.gz
/path/to/project/cell-2/*_R2_*.fastq.gz
```

```sh
# all cells in the same folder
/path/to/project/cell-1_*_R1_*.fastq.gz
/path/to/project/cell-1_*_R2_*.fastq.gz
/path/to/project/cell-2_*_R1_*.fastq.gz
/path/to/project/cell-2_*_R2_*.fastq.gz
```

either of which would find and analyze two cells, named cell-1 and cell-2.

Optionally, you may also provide a path to a Michigan Advanced
Genomics Core (AGC) manifest file that provides additional descriptions 
of the samples (ignore this option if your data did not come from the AGC).

```
--manifest-file /path/to/<leader>_DemuxStats.csv
```

### Required Genome Inputs

Additional required files provide information about the
genome in use. Option '--genome' might typically be 'hg38'.
Option '--genomes-dir' tells the pipeline where to find the
required genome reference files, which must include
sub-directory 'iGenomes' properly populated with the named
genome. See:

<https://support.illumina.com/sequencing/sequencing_software/igenome.html>

Additionally, pipeline 'prepareBins' must have previously
been run for the combination of values of inputs '--bin-size',
'--kmer-length' and '--n-errors', with relevant files found in
'\<genomes-dir\>/bins/\<genome\>/fixed_width_bins'.
