---
title: Demo data and code
parent: svCapture
has_children: false
nav_order: 30
published: true
---

## {{page.title}}

We provide a complete working example of a job configuration file and 
associated data set for testing and demonstrating your svCapture installation. 

### Obtain the demo data, scripts, and support files

Working from whatever folder you'd like, download and unpack the demo archive (file size = 135 MB):

```sh
wget https://PENDING/svCapture-demo.tar.gz
tar -xzvf svCapture-demo.tar.gz
rm svCapture-demo.tar.gz
cd svCapture-demo
```

The entire demo will take place in the 'svCapture-demo' directory so you can easily delete it later.

Reads in the FASTQ files were obtained from cell line HCT116 from
a tagmentation svCapture library in which the central
400 kb of the WWOX gene on human chr16 was subjected to probe capture.
Reads were filtered to include only chr16 and downsampled
to 1M read pairs to keep the demo small and fast.

### Install everything else

Follow the [installation instructions](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/code.html)
to create:
- a multi-suite MDI installation
- an alias to the MDI utility called _mdi_

If you choose a different type of installation or don't make an alias, 
please adjust all commands as needed.

Then, [build the required conda runtime environments](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/runtime.html)
and [download the hg38 reference genome](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/genome.html)
into the demo directory by following the installation instructions.
If you install the genome into a different directory, you will need 
to edit the job file.

### Examine the svCapture job configuration file

```sh
cat svCapture-demo.yml
mdi inspect svCapture-demo.yml # check syntax, directories, and report all options
```

Pipeline options are specified in an extended YAML format 
that supports variables and option declarations
common to multiple pipeline actions. See the file comments for details.
 
The demo job file is configured to work entirely
from your working demo directory - change paths when doing real work,
or if you installed the hg38 genome into a different location, above.

### Run the demo pipeline

To execute the demo in the command shell, use:

```sh
mdi svCapture svCapture-demo.yml --dry-run
mdi svCapture svCapture-demo.yml
```

To submit the demo to your cluster job scheduler, use:

```sh
mdi submit --dry-run svCapture-demo.yml
mdi submit svCapture-demo.yml
```

Depending on your HPC server, you might need to specify additional options,
which you can edit into the job file or add at the command line. For 
example, the following would specify a user account that takes precedence 
over any value found in the job file.

```sh
mdi submit svCapture-demo.yml --account <userAcount>
```
