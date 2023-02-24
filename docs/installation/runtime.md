---
title: Build environments
parent: Installation and usage
has_children: false
nav_order: 20
published: true
---

## {{page.title}}

Most MDI pipelines, including svCapture, depend
on third-party programs installed into an appropriate runtime
environment. You may use one of two methods to set up your environment.

### Singularity containers

[Singularity containers](https://sylabs.io/singularity/) are a platform
for code encapsulation that is especially useful on shared computing resources
such as many HPC servers (unlike Docker which is not allowed on most shared servers).

For servers and pipelines that support publicly shared Singularity containers you do not
need to do anything except to agree to downloading the required container if prompted. 

### Conda environments

Alternatively, all pipelines can be configured using [conda](https://docs.conda.io/en/latest/),
which must be installed and available on your system. You then call
the following MDI command(s) to build the required pipeline environement(s).

```sh
mdi download conda --create   # for downloading genome(s), if required
mdi svCapture conda --create  # or similar for another required pipeline
```

### Defaults and overrides

The MDI utility will recognize if you have Singularity on your server
and if the requested pipeline supports containers. If so, it will use containers
by default. Otherwise, you must build your own conda 
environment as described above.  You can force the use of locally built 
environments by setting the `--runtime` option, e.g.,

```sh
mdi svCapture align --runtime conda
```
