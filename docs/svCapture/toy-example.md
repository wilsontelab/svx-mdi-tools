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

Working from whatever folder you'd like on a Linux computer, 
download and unpack the demo archive (file size = 270 MB):

```sh
wget https://mdi-demo.wilsonte-umich.io/files/svx-mdi-tools/svCapture-demo2.tar.gz
tar -xzvf svCapture-demo2.tar.gz
rm svCapture-demo2.tar.gz
cd svCapture-demo2
```

The entire demo will take place in the 'svCapture-demo2' directory so you can easily delete it later.

Reads in the FASTQ files were obtained from cell line GM12878 from
tagmentation svCapture libraries in which the central
400 kb of the WWOX gene on chr16 was subjected to probe capture.
Reads were filtered to include only chr16 and downsampled
to 1M read pairs per sample to keep the demo small and fast. One sample (sv_high)
was induced to have a higher SV burden than the other sample (sv_low).

### Install everything else

Follow the 
[installation instructions](https://wilsontelab.github.io/svx-mdi-tools/docs/20_installation/10_code.html)
to create:
- a multi-suite MDI installation
- an alias to the MDI utility called _mdi_

If you choose a different type of installation or don't make an alias, 
please adjust all commands below as needed.

Next, 
[build the required conda runtime environments](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/runtime.html)
by following the installation instructions. All required support software is installed in these environments
with appropriate versions.

Finally, 
[download the hg38 reference genome](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/genome.html)
into the demo directory by following the installation instructions. If desired, the demo 
includes a job file, `download-hg38.yml`, which will properly install the required genome into the demo folder
so you may subsequently easily delete it.

```sh
mdi download download-hg38.yml --dry-run # check the configuration but don't do anything
mdi download download-hg38.yml
```

It will take many minutes to complete the download and extraction.

If you install the genome into a different directory, you will need to edit the svCapture job file.

### Explore the svCapture job configuration file

```sh
cat svCapture-demo2.yml
mdi inspect svCapture-demo2.yml # check syntax, directories, and report all options
```

Pipeline options are specified in an extended YAML format 
that supports variables and option declarations
common to multiple pipeline actions. See the file comments for details.
 
The demo job file is configured to work entirely
from your working demo directory - change paths when doing real work,
or if you installed the hg38 genome into a different location.

### Run the demo pipeline

To execute the demo in the command shell, use:

```sh
mdi svCapture svCapture-demo2.yml --dry-run # check the configuration but don't do anything
mdi svCapture svCapture-demo2.yml
```

To submit the demo to your server cluster job scheduler, use:

```sh
mdi submit --dry-run svCapture-demo2.yml
mdi submit svCapture-demo2.yml # add options such as --account if needed on your server
```

It took us about 10 min total wall time to run the demo.

Depending on your HPC server, you might need to specify additional options,
which you can edit into the job file or add at the command line. For 
example, the following would specify a user account that takes precedence 
over any value found in the job file.

```sh
mdi submit svCapture-demo2.yml --account <userAcount>
```

### Examine the pipeline log

If you ran svCapture as a command above, the log stream will have printed to stdout.

If you submitted the demo to your job scheduler, you can view job status
and a log report using commands:

```sh
mdi status svCapture-demo2.yml
mdi report -j all svCapture-demo2.yml
```

[This repository file](https://github.com/wilsontelab/svx-mdi-tools/blob/main/docs/svCapture/svCapture-demo2.log) 
shows a complete log of the demo job sequence as it executed on our server.

### Access and visualize the results

The output of the demo pipeline will be in folder `./output`. 

```sh
ls -l output
```

You may examine the vcf files in the sample folders, etc.

Files 
`samples/samples.svCapture.find.mdi.package.zip` 
and 
`assembled_samples/assembled_samples.svCapture.assemble.mdi.package.zip`
are data packages ready to be uploaded into the svCapture R Shiny app. 

To quickly see the basics, you can view the results of the svCapture demo in our
[publicly accessible demo server](https://mdi-demo.wilsonte-umich.io/). 
Login using the passphrase `mdi-demo` and 
use the `Load from Server` button to find the svCapture demo packages and bookmarks.

The best way to run the app server yourself is by 
[installing the MDI Desktop](https://midataint.github.io/mdi-desktop-app/docs/installation).

Please be aware that the app will be functional but not very informative when looking at
the much reduced demo data set.
