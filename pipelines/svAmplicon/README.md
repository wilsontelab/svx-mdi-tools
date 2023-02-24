---
published: false
---

## svAmplicon: Structural variant analysis in pooled amplicon libraries

**svAmplicon** takes the output of paired-end whole-genome
sequencing data derived from one or more pooled PCR amplicons
and uses it to call structural variants, with 
downstream visualization in a Stage 2 Shiny app.
Typical applications are in the study of targeted mutagenesis
by induced or natural DNA double-strand breaks.

The unique features of the svAmplicon pipeline are first that it expects
many reads to share the same outer molecule endpoints, i.e., to have
the expected endpoints derived from the genomic primer binding sites.
Thus, source molecules cannot be distinguished as unique based
on their outer endpoints.

Moreover, svAmplicon expects junctions that have been sequenced
completely, i.e., SV junctions in the gaps between read pairs
are less valuable. Therefore, it is imperative with paired-end 
sequencing that the read length be at least half as long as the
the longest expected molecule to be sequenced. Molecules with longer
internal insertions will be identified by not characterized further.

### Required data inputs

The svAmplicon pipeline will align your reads to the genome for you.
In this case, you must provide a set of fastq.gz format files, 
with one pair of files for each sample corresponding to the two 
reads of a paired-end sequencing run. These are specified using 
the '--input-dir' and '--input-name' options of the 'align' action.
Your job configuration should execute actions 'align', 'extract'
and 'find'.

Alternatively, if you pre-aligned your reads in some other pipeline,
you can skip the 'align' action and provide the path to a 
name-sorted (not coordinate-sorted) bam or cram file using the
'--bam-file' and '--use-cram' options of the 'extract' action.
Your job configuration should execute actions 'extract' and 'find'.

### Required genome inputs

Additional required files provide information about the
genome in use. Option '--genome' might typically be 'hg38'.
Option '--genomes-dir' tells the pipeline where to find the
required genome reference files, which must include
sub-directory 'iGenomes' properly populated with the named
genome. See:

<https://support.illumina.com/sequencing/sequencing_software/igenome.html>
