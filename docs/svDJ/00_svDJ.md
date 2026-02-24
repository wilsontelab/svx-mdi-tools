---
title: svDJ
has_children: true
nav_order: 60
published: true
---

## {{page.title}}

The svDJ pipeline and app analyzes long-read amplicon sequencing
libraries for structural variant (SV) junctions
relative to a reference genome.

Early pipeline actions support file conversion and offline basecalling of
Oxford Nanopore long-read data.

Later actions support alignment to a reference genome and SV junction
discovery, collation, and inter-molecule comparison. These actions
might reasonably be applied to any long read data with minor modifications,
but at present have only been applied to nanopore data.

### Expectations
svDJ assumes that each input library of molecules represents amplification
of a **single amplicon family** subjected to adapter ligation and sequencing.
A set of such libraries may be barcoded and sequenced together on a
single flow cell, but upstream processing (e.g., by MinKNOW) must have
already demultiplexed the libraries into separate file sets.
Thus, there should be separate folders that each hold one or more sequence
files with reads from one or more PCR reactions performed with one forward
and one reverse primer. It is acceptable, in fact expected, that those primers
will have amplified more than more one type of source molecule.
In the future, it could be possible to have 
mixed amplicons, i.e., primer sets, but this is not presently supported.

As its name implies, svDJ was developed to analyze V(D)J recombination
junctions. However, it could be applied equally well to any amplicon
configuration that matches the expectations described above and below.

### Sequence read configurations

**Single-junction reads**

The following diagram illustrates the expected amplicon structure:

<pre style="width: 1000px;">
outerNode1  junctionNode1  junctionNode2  outerNode2 
|                       |  |                       | 
*=========>-------------*//*-------------<=========* 
primer1               junction               primer2 
</pre>

Importantly, svDJ will only yield proper outcomes when there is one
junction in between two flanking alignments. V(D)J amplicons
satisfy this expectation when D segments are too small to be aligned
to the referene on their own. The D segments are found within the 
insertion bases at a junction called from the V to J segments.

Aberrant and more complex molecules are also invariably present, such as:

**Multi-junction reads**, such as might arise by novel biological insertions of other genomic segments, or by PCR artifact.

<pre style="width: 1000px;">
outerNode1  junctionNode1  junctionNode2 junctionNode3  junctionNode4  outerNode2 
|                       |  |                         |  |                       | 
*=========>-------------*//*-------------------------*//*-------------<=========* 
</pre>
**Chimeric amplicon fusions**, resulting from either intermolecular ligation or miscalling of two DNA molecules in a single read:

<pre style="width: 1000px;">
outerNode1  junctionNode1  junctionNode2  falseNode1  falseNode2  junctionNode3  junctionNode4  outerNode2 
|                       |  |                       |  |                       |  |                       | 
*=========>-------------*//*-------------<=========*//*=========>-------------*//*-------------<=========* 
</pre>
**Duplex reads**, a special case of chimeric amplicon fusion in which the two strands of a single source molecule were both sequenced and fused (note the inverted repetition of the same junction nodes):

<pre style="width: 1000px;">
outerNode1  junctionNode1  junctionNode2  falseNode1  falseNode2  junctionNode2  junctionNode1  outerNode2 
|                       |  |                       |  |                       |  |                       | 
*=========>-------------*//*-------------<=========*//*=========>-------------*//*-------------<=========* 
</pre>
**Truncated reads**, resulting from strand breaks, etc.:

<pre style="width: 1000px;">
outerNode1  junctionNode1  junctionNode2   outerNode2 
|                       |  |               | 
*=========>-------------*//*-------------<=* 
</pre>

svDJ seeks to:
- initially keep and characterize single-junction and multi-junction reads 
- split chimeric amplicon fusions into their respective parts, which are each kept and analyzed as independent molecules
- keep just one half of duplex reads, to prevent double-counting of non-independent source molecules
- discard truncated and off-target reads as untrustworthy

### Read processing
To facilitate analysis and characterization of all molecule types, svDJ uses
a node-edge graph structure. Each aligned segment
of a read has two node positions, one at each end (see above). These nodes define two kinds
of edges between them:

- **alignments (A)**, the expected reference genome contiguity found by the aligner
- **junctions (J)**, the unexpected fusion of two distant genomic positions, called as a split read by the aligner

Alignments have properties such as the CIGAR strings that describe their
exact base matches to the reference genome.
Junctions have properties such as the nature of any
microhomologies and inserted bases.

During processing by action <code>extract</code>, each read is converted to a set of rows in a table of <code>edges</code>, as follows:

<table>
<thead>
<tr>
<th>readI</th>
<th>edgeN</th>
<th>edgeType</th>
<th>node1</th>
<th>node2</th>
<th>additional columns</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>1</td>
<td>A</td>
<td><strong>12</strong></td>
<td>34</td>
<td>alignment properties</td>
</tr>
<tr>
<td>1</td>
<td>2</td>
<td>J</td>
<td>34</td>
<td>-78</td>
<td>junction properties</td>
</tr>
<tr>
<td>1</td>
<td>3</td>
<td>A</td>
<td>-78</td>
<td><strong>-56</strong></td>
<td>etc.</td>
</tr>
<tr>
<td>2</td>
<td>1</td>
<td>A</td>
<td><strong>56</strong></td>
<td>78</td>
<td>etc.</td>
</tr>
<tr>
<td>2</td>
<td>2</td>
<td>J</td>
<td>78</td>
<td>-34</td>
<td>etc.</td>
</tr>
<tr>
<td>2</td>
<td>3</td>
<td>A</td>
<td>-34</td>
<td><strong>-12</strong></td>
<td>etc.</td>
</tr>
<tr>
<td>3</td>
<td>1</td>
<td>etc.</td>
<td></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

where:
- a read always has an odd number of edges of the form A-J[-A-J...]-A
- node positions are 64-bit base coordinates along the entire concatenated genome
- negative node positions denote reads corresponding to the bottom reference strand
- most reads are expected to have one of the two possible strandedness configurations of the same outer node positions, as fixed by the amplification primers (shown in bold)

After extracting the complete edges table, svDJ <code>parse</code>:
- performs <em>ab initio</em> discovery of the two primer outer node positions 
- splits chimeric reads at internal junctions that match abutted primer nodes, to yield one or more <code>segments</code>
- discards segments where the outer nodes do not match the expected primer1/primer2 amplicon
- discards one split segment when a duplex chimeric read repeats the same SV junction(s) in opposite orientation
- performs exact junction matching between segments to highlight biological or technical replicates
- assembles unique junction sequences into networks of similar junction sequences (see below)

### Matching nodes and junctions
Importantly, sequencing errors lead to imprecision in the assigned values of node positions.
svDJ therefore sometimes uses an adjancency tolerance when matching nodes. For example,
molecules with outer nodes 12/-56 and 12/-57 would both be considered a match to
primer nodes 12/-56.

In its first layer of junction matching, <code>parse</code> uses exact matching based on the
actual node positions and any novel inserted bases. This yields the following table of
<code>junctions</code>.

<table>
<thead>
<tr>
<th>jxnKey</th>
<th>nMatchingSegments</th>
<th>node1</th>
<th>node2</th>
<th>insertSize</th>
<th>jxnSeq</th>
<th>etc.</th>
</tr>
</thead>
<tbody>
<tr>
<td>xxx</td>
<td>23</td>
<td>34</td>
<td>-78</td>
<td>4</td>
<td>CGTA</td>
<td></td>
</tr>
<tr>
<td>yyy</td>
<td>1</td>
<td>33</td>
<td>-78</td>
<td>5</td>
<td>ACGTA</td>
<td></td>
</tr>
<tr>
<td>etc.</td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

However, this might undercall true replicates that diverged due to PCR or
sequencing errors and are thus still found in different junction rows.
The following examples illustrate some of the ways that even a single basecalling error
can lead to a substantial change in the apparent junction call:

<pre>
correct junction, 1-base microhomology:
    alignment 1:   -----*
    read sequence: ~~~~~A~~~~~
    alignment 2:        *-----
junction with A to T base error, called as 1-base insertion:
    alignment 1:   ----*
    read sequence: ~~~~~t~~~~~
    alignment 2:         *----
</pre>
<pre>
correct junction, 3-base insertion:
    alignment 1:   ----*
    read sequence: ~~~ATacg~~~~~
    alignment 2:           *----
junction with missing base "A", called as 4-base insertion:
    alignment 1:   --*
    read sequence: ~~~Tacg~~~~~
    alignment 2:          *----
</pre>

To help match junctions impacted by these factors, a second
round of matching is performed based on the edit distance between junctions,
irrespective of the node positions and insert size.
The logic is that highly abundant junctions (based on exact matching)
are likely to have spawned some reads that diverged from the parent molecule
due to a method error at the junction.

The <code>parse</code> action builds
adjacency networks by analyzing a graph where 
all unique junctions above a user-defined coverage threshold
are the nodes and edge weights are the edit distance
between them. The graph is directed, such that lower frequency
junctions point toward higher frequency junctions. 

The most abundant junction is declared as the parent node
of a first subgraph expected to have less frequent derivative 
junctions that arose by error processes.
Each less frequent junction node is then analyzed to assess 
whether it is best modeled as a new parent node representing
a distinct junction or as a child, i.e., derivative, of
an existing parent node. This decision is controlled
by user parameters of the maximum edit distance from a parent
to a child and how much to weight the junction coverage.
Larger edit distances and higher junction coverage favor
calling a query node as the parent of a new junction subgraph.

Once all candidate junction nodes have been evaluated, `parse` can
continue on, if instructed, to add the lowest coverage nodes to 
existing junction networks, but they will not create new junction calls.
Because this can entail processing large numbers of singleton sequences,
the default behavior is to skip this process as it is slow and often has little impact
on the final result.

Importantly, only the bases at the junction are used when comparing junctions to each other.
Base differences in the amplicon "stems" on each side of the junction are assumed to be 
sequencing errors and are disregarded by assembling a "fake" representative sequence for each read segment,
where the flanks are replaced with the corresponding amount of reference genome sequence 
(i.e., from outer primer node to the junction node).

Also important is the fact that a given read might correspond to either strand
of a source DNA amplicon molecule, which must be accounted for during junction matching.
svDJ achieves this by defining a **canonical orientation** for every junction, i.e.,
for every possible combination of node pairs. A read segment is on the canonical strand
if it begins at the primer1 node, as defined by the primers table. Noncanonical segments
are reverse complemented prior to junction comparison. Notably, the caononical
orientation of a read is not necessarily in V(D)J order, i.e., primer 1 might be on
the J side of the amplicon depending on the orientation of the V and J segments
in the genome.

The resulting final level of aggregation is thus the table of junction <code>networks</code>, as follows.
All tables (primers, edges, junctions, and networks) are returned as RDS files for further
exploration.

<table>
<thead>
<tr>
<th>networkKey</th>
<th>nMatchingJunctions</th>
<th>nMatchingSegments</th>
<th>node1</th>
<th>node2</th>
<th>insertSize</th>
<th>jxnSeq</th>
<th>etc.</th>
</tr>
</thead>
<tbody>
<tr>
<td>xxx</td>
<td>2</td>
<td>24</td>
<td>34</td>
<td>-78</td>
<td>4</td>
<td>CGTA</td>
<td></td>
</tr>
<tr>
<td>etc.</td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

</body>
</html>
