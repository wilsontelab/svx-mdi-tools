---
title: svWGS
has_children: true
nav_order: 40
published: true
---

## {{page.title}}

The svWGS pipeline and app analyze paired-end, short-read whole genome sequencing (WGS)
data to find individual high-confidence SV junctions. 

A premium is placed
on building a comprehensive evidence base to support individual junction calls
as being real above the ever-present baseline SV artifact levels.  This might include either
clonal or mosaic SVs with multi-molecule evidence, or ultra-rare mosaic molecules
only ever expected to be found once in the input library.

While the focus is on careful characterization of junctions, the
pipeline and app also monitor copy number status to support
clonal deletion and duplication junctions by parallel CNV detection. 

svWGS is not as useful for finding stretches of DNA containing multiple SVs
and characterizing their paths - the reads are too short for this (see svPore). 
It is also not generally intended to produce
a final comprehensive "true call set" of all SVs in a sample, although one
can certainly be produced after applying some desired set of quality filters
to the complete table of all SV junctions that is produced.
