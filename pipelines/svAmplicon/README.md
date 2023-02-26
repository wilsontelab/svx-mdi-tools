---
published: false
---

## svAmplicon: Structural variant analysis in pooled amplicon libraries

**svAmplicon** takes the output of paired-end whole-genome
sequencing data derived from one or more pooled PCR amplicons
and uses it to call structural variants (SVs), with 
downstream visualization in a Stage 2 Shiny app.
Typical applications are in the study of targeted mutagenesis
by induced or natural DNA double-strand breaks.

The distinct feature of svAmplicon relative to other svx pipelines
is that it expects
reads to share common outer molecule endpoints, i.e., to have
the expected endpoints derived from the genomic primer binding sites.
Thus, source molecules cannot be distinguished as unique based
on their outer endpoints.

### Operating modes and assocated experimental designs

svAmplicon is built on the same algorithms used for _de novo_ genomic
SV calling so is highly generalized. However, only two experimental
approaches will yield optimal results.

The more common approach is to use **overlapping paired reads** that
sequence an amplicon with relatively close primer positions.
Here, svAmplicon expects junctions to have been sequenced
completely. Therefore, it is imperative that the read length be more than 
half as long as the longest expected molecule. Molecules with longer
internal insertions will be identified by not characterized in detail.
This mode is called **expectOverlap** in amplicon reports.

Less commonly, you might design **distant primer pairs** that 
are so widely spaced that you don't expect unrearranged genomic
DNA to give you any products. Thus, any sequenced products are assumed
to represent SVs relative to the primer design.
This mode is called **notPossible** in amplicon reports.

The pipeline will also process, up to a point, a third configuration
where you do expect to amplify and sequence DNA molecules that are 
too big for the reads to overlap in the middle. This mode is 
discouraged as you will fail to characterize many junctions. The
pipeline won't report as much detail here. 
This mode is called **expectGaps** in amplicon reports.

### Required data inputs

The svAmplicon pipeline aligns your reads to the genome for you.
Accordingly, you must provide a set of fastq.gz format files, 
with one pair of files for each sample corresponding to the two 
reads of a paired-end sequencing run. These are specified using 
the '--input-dir' and '--input-name' options of the 'align' action.
Your job configuration should execute actions 'align' and 'collate'.

It is not possible to use externally provided bam files due to the special 
handling of the alignment process in svAmplicon. Amplicon sequencing
alignment is generally small and fast so this is not a significant limitation.

It is not necessary to provide up-front information on the amplicon(s)
your are sequencing. svAmplicon uses an agnostic approach in which
reads are aligned to the reference genome with heuristic algorithms 
used to perform _de novo_ amplicon discovery. There are a few options you
can set to control the number and quality of found amplicons. 

The amplicon discovery approach makes it routinely possible to sequence
pools of distinct amplicons without barcoding them. The ends of the sequenced
molecules identify the amplicon they correspond to. 
Of course, if you wish to pool the same amplicon from different samples,
they must be barcoded and demultiplexed before running svAmplicon. 

### Required genome inputs

Additional required files provide information about the
genome in use. Option '--genome' might typically be 'hg38'.
Option '--genomes-dir' tells the pipeline where to find the
required genome reference files, which must include
sub-directory 'iGenomes' properly populated with the named
genome. See:

<https://support.illumina.com/sequencing/sequencing_software/igenome.html>
