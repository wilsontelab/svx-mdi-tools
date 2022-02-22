---
title: Stage 1 Pipelines
has_children: true
nav_order: 1
published: true # set to false to remove this tab from your suite's doc site
---

## Stage 1 Pipelines

A **pipeline** (also called a workflow) is a single coordinated set
of data analysis instructions that can be called by name using the
mdi command line utility, e.g.,

```bash
mdi <pipeline> ...
```

A pipeline is distinct from many standalone programs in
that a pipeline might have many inputs and outputs and typically makes 
calls to many installed programs. Thus, rather than being a program itself, 
a pipeline is a means of coordinating a reproducible set of calls to 
other programs.
