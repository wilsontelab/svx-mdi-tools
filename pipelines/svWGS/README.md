---
published: false
---

## svWGS: Structural Variant Analysis by Whole Genome Sequencing

**svWGS** takes the output of typical, paired-end whole-genome
sequencing data and uses it to call structural variants, with 
downstream visualization in a Stage 2 Shiny app.

The pipeline can be applied to one or multiple related samples in 
a run. When multiple samples are analyzed together, the 'find' 
algorithm will provide an implicit comparison of all samples to 
support the identification of SVs unique to one sample.

### Required Data Inputs

The svWGS pipeline can align your reads to the genome for you.
In this case, you must provide a set of fastq.gz format files, 
with one pair of files for each sample corresponding to the two 
reads of a paired-end sequencing run. These are specified using 
the '--input-dir' and '--input-name' options of the 'align' action.
Your job configuration should execute actions 'align', 'extract''
and 'find'.

Alternatively, if you pre-aligned your reads in some other pipeline,
you can skip the 'align' action and provide the path to a 
name-sorted (not coordinate-sorted) bam or cram file using the
'--bam-file' and '--use-cram' options of the 'extract' action.
Your job configuration should execute actions 'extract' and 'find'.

### Required Genome Inputs

Additional required files provide information about the
genome in use. Option '--genome' might typically be 'hg38'.
Option '--genomes-dir' tells the pipeline where to find the
required genome reference files, which must include
sub-directory 'iGenomes' properly populated with the named
genome. See:

<https://support.illumina.com/sequencing/sequencing_software/igenome.html>

Additionally, pipeline 'prepareBins' must have previously
been run exactly as follows, otherwise the 'find' action will 
throw an error:

```yml
---
pipeline: genomex-mdi-tools/prepareBins
variables:
    GENOMES_DIR:    ${MDI_DIR}/resources/genomes
    GENOME:         hg38 # or another genome, as required
do:
  fixed-width-bins:
    bin-size:       65536
    kmer-length:    100
    n-errors:       1
  output:
    output-dir:     ${GENOMES_DIR}/bins
    data-name:      ${GENOME}

execute:
  - do
```

The genomex-mdi-tools suite was co-installed with svx-mdi-tools.
