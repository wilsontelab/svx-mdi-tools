# Wilson lab tools for structural variant analysis

The [Michigan Data Interface](https://midataint.github.io/) (MDI) 
is a framework for developing, installing and running 
Stage 1 HPC **pipelines** and Stage 2 interactive web applications 
(i.e., **apps**) in a standardized design interface.

This **svx-mdi-tools** repository contains pipelines and apps
for genome structural variant (SV) analysis by different library strategies
from the 
[Thomas Wilson laboratory](https://wilsonte-umich.github.io)
at the University of Michigan.

## Available tools and documentation

Pipelines and apps in stable release are:
- **svCapture** = find SVs in whole-genome capture sequencing

Pipelines and apps in beta are:
- **svWGS** = find SVs in whole-genome sequencing data (without capture)
- **scCNV** = find SVs in single-cell whole-genome amplification sequencing
- **svAmplicon** = find SVs in targeted amplicon sequences

Please see the [svx-mdi-tools documentation site](https://wilsontelab.github.io/svx-mdi-tools)
for detailed information on the usage of specific pipeline and apps.

---
## Quick Start Method 1: multi-suite installation (recommended)

You can install MDI tool suites in one of two ways: as a **multi-suite installation** that carries one or more distinct tool suites (recommended), or as a more contained **single-suite installation** dedicated to just the svx-mdi-tools suite.

In the recommended multi-suite mode, you will:
- clone and install the MDI framework
- add svx-mdi-tools (and potentially other suites) to your MDI installation
- call the _mdi_ utility to use tools from any installed suite

### Install the MDI framework

Please read the _install.sh_ menu options and the 
[MDI installer instructions](https://github.com/MiDataInt/mdi.git) to decide
which installation option is best for you. Briefly, choose option 1
if you will only run Stage 1 HPC pipelines from your installation.

```bash
git clone https://github.com/MiDataInt/mdi.git
cd mdi
./install.sh
```

### OPTIONAL: Add an _mdi_ alias to _.bashrc_

These commands will create a permanent named alias to the _mdi_
target script in your new installation.

```bash
./mdi alias --help
./mdi alias --alias mdi # change the alias name if you'd like 
`./mdi alias --alias mdi --get` # activate the alias in the current shell (or log out and back in)
mdi
```

Alternatively, you can add the MDI installation directory to your PATH variable,
or always change into the directory prior to calling _./mdi_.

### Add the svx-mdi-tools suite to your MDI installation

```bash
./mdi add --help
./mdi add -p -s wilsontelab/svx-mdi-tools 
```

Alternatively, you can perform the required suite addition steps one at a time:

```sh
nano config/suites.yml # or use any other text editor to edit suites.yml
```

```yml
# mdi/config/suites.yml
suites:
    - wilsontelab/svx-mdi-tools # add this tools suite to the config file
```

```sh
./install.sh # re-install to add the new tool suite
```



### Execute a Stage 1 pipeline from the command line

For help, call the _mdi_ utility with no arguments, which describes the format for pipeline calls. 

```bash
./mdi  # call the mdi utility directly without an alias, OR
mdi    # if you created an alias as described above
mdi svCapture # or change to one of the other pipelines...
# etc.
```

### Launch the Stage 2 web apps server

To launch the MDI web server, we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

---
## Quick Start Method 2: single-suite installation

In the alternative single-suite mode, you will install just the svx-mdi-tools suite by:
- cloning this tool suite repository
- running _install.sh_ to create a suite-specific MDI installation
- OPTIONAL: calling _alias.pl_ to create an alias to the suite's _run_ utility
- calling the _run_ utility to use a tool from the suite

### Install this tool suite

```bash
git clone https://github.com/wilsontelab/svx-mdi-tools.git
cd svx-mdi-tools
./install.sh
```

### OPTIONAL: Create an alias to the suite's _run_ utility

```bash
perl alias.pl svx # you can use a different alias name if you'd like
```

Alternatively, you can add the installation directory to your PATH variable,
or always change into the directory prior to calling _./run_.

### Execute a Stage 1 pipeline from the command line

For help, call the _run_ utility with no arguments, which describes the format for pipeline calls. 

```bash
./run  # call the run utility directly without an alias, OR
svx    # use the alias, if you created it as described above
```

### Launch the Stage 2 web apps server

To launch the MDI web server, we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.
