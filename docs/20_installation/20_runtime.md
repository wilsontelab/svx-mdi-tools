---
title: Build environments
parent: Installation and usage
has_children: false
nav_order: 20
published: true
---

## {{page.title}}

Most MDI pipelines, including those in svx-mdi-tools, depend
on third-party programs installed into an appropriate runtime
environment. You may use one of two methods to set up your environment.

### Singularity containers

[Singularity containers](https://sylabs.io/singularity/) are a platform
for code encapsulation that is especially useful on shared computing resources
(unlike Docker which is not allowed on many shared HPC servers).

For servers and pipelines that support publicly shared Singularity containers you do not
need to do anything except agree to downloading the required container when prompted. 

The following pipelines currently support singularity containers:
- genomex-mdi-tools/download

### Conda environments

Alternatively, all pipeline environments can be set up using [conda](https://docs.conda.io/en/latest/),
which must be installed and available on your system. You then call
the following MDI command(s) to create/build the required environment(s).

```sh
mdi download conda --create   # for downloading genome(s), if required
mdi svCapture conda --create  # or replace svCapture with another required pipeline
```

### Defaults and overrides

The MDI utility will recognize if you have Singularity on your server
and if the requested pipeline supports containers. If both are true, it will use containers
by default. Otherwise, you must build your own conda 
environment as described above.  

You can force the use of locally built 
environments by setting the `--runtime` option, e.g.,

```sh
mdi svCapture align --runtime conda
```
