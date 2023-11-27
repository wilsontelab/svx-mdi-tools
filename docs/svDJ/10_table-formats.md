---
title: Table Formats
parent: svDJ
has_children: false
nav_order: 10
published: true
---

## {{page.title}}

The following are detailed descriptions of the table formats
used to encode amplicon alignments and SVs. Columns are ordered for ease of understanding the information they encode;
the actual order of the columns in R data objects or text files may differ.
Tables are produced as data.tables in R and saved as RDS files.

### Primers

The `primers` table has exactly two rows with metadata about the 
two primers inferred to have generated the read data. If it is 
a F/R primer pair relative to the genome, the forward (F) primer
is always primer 1 and the reverse (R) primer is primer 2. Otherwise primers are sorted according to
position along the concatenated genome. 

| column | data_type | description |
| ------ | --------- | ----------- |
| primer | integer | the primer number, either 1 or 2; primer 1 defines side 1 (e.g., chrom1) in other tables |
| node | integer64 | the 64-bit signed integer identifying the genome position and strand of the 5'-most primer base |
| chrom | character | the chromosome name |
| chromIndex | integer | the chromosome index |
| refPos | integer | the unsigned reference position on `chrom` |
| strand | character | the strand of the primer, either "+" or "-" |
| sequence | character | the genome reference sequence starting from (if `strand` == "+") or ending at (if `strand` == "-") `refPos` |
| rcSequence | character | reverse complement of `sequence` |
| seqLength | integer | length of `sequence`, equivalent to the longest amplicon sequence from refPos |

### Edges

The `edges` table has one row per edge, i.e., a connection between two nodes, with an odd number of edge rows per
sequenced read segment, in format A-[J-A...]-A. The columns are 
listed here as found in the final edges table, i.e., after splitting of 
reads into segments and other filtering actions. The edges table carries the most comprehensive
information about each read and is accordingly the most complex table.

| column | data_type | description |
| ------ | --------- | ----------- |
| **source identifiers** | | |
| sample | character | the name of the source sample, obtained from option `--data-name` |
| channel | integer | the nanopore channel on the ONT flow cell where the read was obtained |
| pod5File | character | path to the pod5/fast5 file that carries the nanopore sampling data |
| qName | character | name of the read as assigned by nanopore code |
| segmentName | character | extension of `qName` with `segmentN` to provide a unique identifier after chimeras are split |
| **edge numbering** | | |
| readI | integer | unique number of the read over all sample data | 
| segmentN | integer | sequential unique number of the segment within a given read |
| blockN | integer | sequential unique number of a contiguity block within a given segment; blocks may disrupted by low quality, unused, and untrusted junctions |
| edgeN | integer | sequential unique number of this edge within a given segment |
| **reference alignment** | | |
| qStart | integer | start position of the edge in the original read, i.e., query sequence |
| qEnd | integer | end position of the edge in the original read |
| node1 | integer64 | stranded start position of the edge in the reference genome; negative values are on the bottom strand |
| node2 | integer64 |stranded end position of the edge in the reference genome  |
| cigar | character | for alignments, the CIGAR string for the alignment of the query read to the reference genome |
| **edge metadata** | | |
| isCanonical | logical | TRUE if the read corresponds to the canonical orientation of a junction |
| edgeType | character | the kind of edge this is, one of A:alignment, D:deletion, U:duplication, V:inverstion, T:translocation (all values other than "A" are junction edges) |
| eventSize | integer | the length of an alignment on on the reference genome, or the size of the SV associated with a junction |
| insertSize | integer | if positive, the number of non-reference bases found at the junction; if negative, the number of bases of microhomology; if 0, it is a blunt joint |
| jxnSeq | character | the base sequence of the inserted or microhomologous bases |
| cJxnSeq | character | the same as `jxnSeq`, but reverse-complemented if `isCanonical` is FALSE |
| **quality metrics** | | |
| mapQ | integer | alignment MAPQ, i.e., mapping quality, as determined by the aligner; for junctions, the minimum value of the two flanking alignments  |
| gapCompressedIdentity | double | the fraction of reference bases that matched the query as an indicator of sequence quality; for junctions, the minimum value of the two flanking alignments |
| baseQual | double | the average base quality score over all bases in the edge |
| alnBaseQual | double | for junctions, the minimum `baseQual` of the two flanking alignments |
| alnSize | integer | for junctions, the minimum `eventSize` of the two flanking alignments |
| passedBandwidth | logical | whether a series of junctions including this junction displaced the continuity of the read on the reference genome by more than `--min-sv-size` |
| **read summaries** | | |
| nTotalJunctions | integer | the number of junctions found in the read containing this edge |
| nKeptJunctions |integer  | the number of high quality junctions in the read containing this edge that were used for junction analysis |
| nKeptSegments | integer | the number of segments contributed by this read after adapter splitting and quality filtering |
| **junction comparison** | | |
| jxnKey | character | a unique identifier for this distinct junction sequence, shared with all other junction edges with the same sequence |
| nSegments | integer | the number of edges with the same `jxnKey` |
 
Edge columns `clip5/3`, `score5/3`, `nBases5/3`, `start5/3`, `end5/3`, and `hasAdapter5/3` are used in the svPore pipeline for finding adapter
sequences in reads. They should be ignored in svDJ tables.
Edge columns `node1/2IsPrimer` and `isChimeric` were used for filtering the edges table prior to the final output. They have constant values and can also be ignored. Finally, column `nStrands`
is not consistently used by svDJ and should be ignored.

### Junctions

The `junctions` table has one row for every distinct junction sequence
found in the edges table, after applying junction quality filters.
Junction edges aggregated into the same junction row did _not_ 
necessarily have the same sequence throughout the entire read (where base differences are most likely errors),
but _did_ have the same sequence at the junction itself, as defined by the
combination of `node1`, `node2`, `insertSize` and, if insertSize > 0, `jxnBases`. 

| column | data_type | description |
| ------ | --------- | ----------- |
| **junction identifiers** | | |
| jxnKey | character | a unique identifier for this distinct junction sequence, matching the same column in `edges` |
| **junction metadata** | | |
| nMatchingSegments | integer | the number of read segments that contained this junction |
| nCanonical | integer | `nMatchingSegments` where `isCanonical` is TRUE |
| nNonCanonical | integer | `nMatchingSegments` where `isCanonical` is FALSE |
| **junction structure** | | |
| node1 | integer64 | stranded start position of the junction in the reference genome; negative values are on the bottom strand |
| node2 | integer64 |stranded end position of the edge in the reference genome  |
| insertSize | integer | if positive, the number of non-reference bases found at the junction; if negative, the number of bases of microhomology; if 0, it is a blunt joint |
| jxnSeq | character | the base sequence of the inserted or microhomologous bases, on the canonical strand |
| chrom1 | character | the chromosome name on side 1 of the junction |
| chromIndex1 | integer | the chromosome index on side 1 of the junction |
| refPos1 | integer | the unsigned reference position on `chrom1` |
| strand1 | character | the strand on side 1 of the junction, either "+" or "-" |
| chrom2 | character | the chromosome name on side 2 of the junction |
| chromIndex2 | integer | the chromosome index on side 2 of the junction |
| refPos2 | integer | the unsigned reference position on `chrom2` |
| strand2 | character | the strand on side 2 of the junction, either "+" or "-" |
| fakeSeq | character | an idealized sequence for this read where flanking alignments are replaced with reference genome bases but the junction itself is as it was sequenced |
| fakeLength | integer | number of bases in fakeSeq, best estimate of the original source molecule length |
| **network membership** | | |
| junctionI | integer | sequential number of this junction, sorted by decreasing value of `nMatchingSegments` |
| parentJunctionI | integer | `junctionI` of the parent junction of the network to which this junction belongs; NA for low frequency junctions not added to any network |
| parentDistance | integer | the number of edits from this junction to the parent junction of the network; NA for low frequency junctions not added to any network |
| junctionOnParent | character | map of this junction sequence onto its network parent, a string of same length as parent `fakeSeq` | 

### Networks

The `networks` table further aggregates the junctions table after the 
construction of junction networks. The junction rows aggregated into
a given network are inferred to have arisen from the same
source molecule but diverged due to PCR or sequencing errors. 
There may be some "collisions" where two independent source molecules
fortuitously had a small edit distance and were inappropriately placed into the same network. 
This outcome is infrequent due to the consideration of junction coverage during network analysis,
but when it happens will reduce the true molecular complexity of the source sample in the networks table. 

However, especially in the case of nanopore sequencing with its higher error rate, it is also likely that some reads diverged from their parent source molecules because they acquired multiple errors. Such reads fail
to be merged into their true network, falsely increasing the apparent molecular complexity of the source sample. Importantly, these often result in junction networks with low coverage that can be filtered against. They also likely arise from lower quality reads, and can be filtered against by option `--min-alignment-identity`, which demands that the alignments flanking a junction
had a high percentage of bases that matched the reference genome (implying that the junction sequence itself was also high quality).

It is up to the user to decide whether junctions or networks give the most
reliable quantification of the data, whether to trust low-coverage networks, how aggressively to set quality filters, etc. In general, we think the network analysis is robust such that the
networks table provides the fairest assessment of detected junctions and their relative abundance. 

| column | data_type | description |
| ------ | --------- | ----------- |
| networkKey | character | the `jxnKey` of the parent junction of the network |
| parentJunctionI | integer | the `junctionI` of the parent junction of the network |
| maxDistance | integer | max(`parentDistance`) over all junctions in the network |
| nMatchingJunctions | integer | number of `junctions` rows aggregated into this network | 
| parentNMatchingSegments | integer | `nMatchingSegments` of the parent junction | 
| nextNMatchingSegments | integer | `nMatchingSegments` of the first non-parent junction added to the network |
| nMatchingSegments | integer | sum(`nMatchingSegments`) over all contributing junctions |
| nCanonical to fakeLength | various | as defined for the `junctions` table |
