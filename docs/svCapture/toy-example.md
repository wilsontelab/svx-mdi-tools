---
title: Example data and code
parent: svCapture
has_children: false
nav_order: 10
published: true
---

## {{page.title}}

We provide a complete working example
of a pipeline configuration file and associated data set
for testing your installation.

This section assumes you have successfully installed
the svx-mdi-tools script repository into your MDI installation
as described here:
- <>

### Obtain the data

Example data can be downloaded from Mendeley Data:
- <>

Reads in these fastq files were obtained from 
a tagmentation svCapture library in which the central
250 kb of the WWOX gene on human chr16 were subjected to probe capture.
Reads were filtered to include only chr16 and downsampled
to keep things small and fast for demonstration purposes.

### Create and examine the job configuration file

The svCapture pipeline can be executed directly from the command line
or by constructing a YAML-format job configuration file, which is recommended.

Copy and paste these lines to create your own demo job file.

```yml
# demo.yml

```

Notice that all pipeline options are specific in YAML format in sections
whose name matches the corresponding pipeline action. The last section
labeled 'execute' names the specific actions that will be executed 
when the pipeline is submitted.

You will need to edit the file to include appropriate paths to
the genome reference fasta file and data output paths.

The following command will print a template you can use to construct a job file.

```sh
mdi svCapture template
```

Complete details of constructing job files are described here:
- <>

As noted, you can also call the pipeline directly from the command line
in a command such as:

```sh
mdi svCapture --xxxx
```
