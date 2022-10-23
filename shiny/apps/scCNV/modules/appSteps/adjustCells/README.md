---
title: adjustCells
parent: App Steps
has_children: false
nav_order: 10
---

## {{page.title}}

The **adjustCells** appStep allows users to:

- remove cells from further consideration, e.g., replicating or otherwise bad cells
- change the modal autosome copy number
- mark cells as doublets to allow proper copy number calls

Visualization plots help identify cells with patterns
that are difficult to discern algorthmically but that 
are easily identified by a trained eye, most importantly
replicating cells.

Only cells with a "keep" designation are moved forward 
into CNV calling.
