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
URL=https://data.mendeley.com/public-files/datasets/38zrkbsbph/files/8b11a546-dac0-4610-9da4-1278936bb16f/file_downloaded
wget ${URL} -o svCapture-demo.tar.gz
tar -xzvf svCapture-demo.tar.gz
rm svCapture-demo.tar.gz
cd svCapture-demo
```

As an alternative to `wget`, you can also download the tarball here:
- <https://data.mendeley.com/datasets/38zrkbsbph>

The entire demo will take place in the 'svCapture-demo' directory so you can easily delete it later.

Reads in the FASTQ files were obtained from a commercially available human cell line from
a tagmentation svCapture library in which the central
400 kb of the WWOX gene on chr16 was subjected to probe capture.
Reads were filtered to include only chr16 and downsampled
to 1M read pairs to keep the demo small and fast.

### Install everything else

Follow the [installation instructions](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/code.html)
to create:
- a multi-suite MDI installation
- an alias to the MDI utility called _mdi_

If you choose a different type of installation or don't make an alias, 
please adjust all commands below as needed.

Then, [build any required conda runtime environments](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/runtime.html)
and [download the hg38 reference genome](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/genome.html)
into the demo directory by following the installation instructions.
Building conda environments is only required if you do not have Singularity available
on your server as the svCapture pipeline supports containers.
If you install the genome into a different directory, you will need 
to edit the job file.

### Explore the svCapture job configuration file

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
mdi svCapture svCapture-demo.yml --dry-run # check the configuration but don't do anything
mdi svCapture svCapture-demo.yml
```

To submit the demo to your cluster job scheduler, use:

```sh
mdi submit --dry-run svCapture-demo.yml
mdi submit svCapture-demo.yml # add options such as --account if needed on your server
```

Depending on your HPC server, you might need to specify additional options,
which you can edit into the job file or add at the command line. For 
example, the following would specify a user account that takes precedence 
over any value found in the job file.

```sh
mdi submit svCapture-demo.yml --account <userAcount>
```

### Examine the pipeline log

If you ran svCapture as a command above, the log stream will have printed to stdout.

If you submitted the demo to your job scheduler, you can view job status
and a log report using commands:

```sh
mdi status svCapture-demo.yml
mdi report -j all svCapture-demo.yml
```

[This repository file](https://github.com/wilsontelab/svx-mdi-tools/blob/main/docs/svCapture/svCapture-demo.log) 
shows a complete log of the demo job sequence as it executed on our server.

### Access and visualize the results

The output of the demo pipeline will be in folder `.../svCapture-demo/svCapture-demo`. 

```sh
ls -l svCapture-demo
```

You may examine the vcf file, etc.

File `svCapture-demo.svCapture.find.mdi.package.zip` is the data package ready
to be uploaded into the svCapture R Shiny app. The best way to run the app server
is by [installing the MDI Desktop](https://midataint.github.io/mdi-desktop-app/docs/installation).

To quickly see the basics, you can view the results of the svCapture demo in our
[publicly accessible demo server](https://mdi-demo.wilsonte-umich.io/). 
Login using the passphrase `mdi-demo` and 
use the `Load from Server` button to find the svCapture demo package and bookmark.

Please be aware that the app will be functional but not very informative when looking at
the very much reduced demo data set.
