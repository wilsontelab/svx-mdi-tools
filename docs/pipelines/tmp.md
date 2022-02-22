---
title: Temporary
parent: Stage 1 Pipelines
has_children: false
nav_order: 0
---

{% include table-of-contents.md %}

### Pipeline actions

Each MDI pipeline must have one or more **actions** defined
in 'pipeline.yml' and organized into pipeline subfolders.
Actions are entered by users at the command line or in
'\<data\>.yml' configuration files.

```yml
# pipeline.yml (the pipeline's configuration file)
actions:
    actionName:
        ... # see other documentation for what goes here
```

```bash
# command line
mdi <pipeline> <actionName> ...
mdi myPipeline do ... # e.g., a single-action pipeline called 'myPipeline'
```

```yml
# <data>.yml (the specifications for analysis of a specific data set)
pipeline: myPipeline
options: 
    ... # see other documentation for what goes here
execute: # the list of actions to execute
    - do # 'do' is the standardized name for a single action
```

### Action steps

Each pipeline action might finally be executed as a series of 
sequential **steps**, organized into action subfolders.
Steps are executed in series from a single pipeline action call.
Their success can be independently monitored to allow a 
pipeline action to restart without repeating previously satisfied 
steps, e.g., if a server crashes during step 3, steps 1 and 2 
would not need to be repeated.

### Conda environments

All pipelines use conda to construct an appropriate execution
environment with proper versions of all required program
dependencies, for explicit version control, reproducibility
and portability. These might include any program called from 
the command line or a shell script to do data analysis work.

### Singularity containers

When implemented by a tool suite developer, any pipeline can
be wrapped into an optional Singularity container to provide complete
control over the entirety of the operating system and codebase
available to support the pipeline. The operating system, e.g.,
Ubuntu, and system libraries are specified in 'singularity.def',
while program dependencies are provided by conda environments
pre-installed into the container.



## Pipeline construction

Create one folder in '\<suite\>/pipelines' for each distinct data 
analysis pipeline carried in your suite. Each pipeline 
might define its own specific code and/or use code elements 
in the '\<suite\>/shared' folder.

Only one file is essential and must be present, called 'pipeline.yml'.

Optional files include (see the _template pipeline for usage):
- 'README.md', to document the pipeline 
- 'pipeline.pl', which can be used to set custom environment variables 
- 'singularity.def', if your pipeline will offer or require a container to run

### Pipeline structure and definition

Begin by editing file **pipeline.yml**, which is the configuration file that establishes your pipeline's identity, options, actions, etc. It dictates how users will provide information to your pipeline and where the pipeline will look for supporting scripts and definitions.

Edit file **README.md** to provide a summary of the purpose and usage of your new pipeline. You will usually want to link this file into your GitHub Pages documentation by adding Jekyll front matter.

#### Pipeline actions

Next, create a subfolder in your pipeline directory for each discrete **action**
defined in pipeline.yml. Many pipelines only require one action, which by convention is called 'do'. Alternatively, you might need multiple actions executed independently, e.g., a first action 'analyze' applied to individual samples followed by a second action 'compare' that integrates information from multiple samples.

By convention, the target script in an action folder is called '**Workflow.sh**' - 
it is the script that performs the work of the pipeline action.
There are no restrictions on exactly how Workflow.sh does its work. You may incorporate 
other code, make calls to programs, including calls to nested workflow managers such as snakemake, etc. Thus, one common pattern might be:

```bash
# pipelines/<pipeline>/<action>/Worflow.sh
TARGET_FILE=$DATA_NAME.XYZ
snakemake $SN_DRY_RUN $SN_FORCEALL \
    --cores $N_CPU \
    --snakefile $ACTION_DIR/Snakefile \
    --directory $TASK_DIR \
    $TARGET_FILE
checkPipe
```




### Output conventions

Data files written by a pipeline are always placed into a folder you specify using 
options '--output-dir' and '--data-name'. File names should always be prefixed with the value of option
'--data-name', such that pipeline output files follow the pattern:

```
<output-dir>/<data-name>/<data-name>.XXX
```

'--output-dir' and '--data-name' are thus universally required options, and their 
values must be the same between sequential actions applied to the same data. 



## Pipeline execution

Once the MDI is installed, pipelines can be called as executables from the command line, 
either by (i) a direct call to the pipeline or (ii) using the manager 'submit' command to queue
a pipeline job via a cluster server's job scheduler.

```bash
mdi <pipeline> ...
mdi submit <data.yml> ...
```

Use the --help options to learn about the different ways to configure pipeline calls.

```bash
mdi --help
mdi <command> --help
mdi <pipeline> --help
mdi <pipeline> <action> --help
```

### Job configuration files (specifying pipeline options)

You may provide all options through the command line 
as you would for any typical program.  Use '--help' to
see the available options. 

However, we recommend instead writing a '\<data\>.yml'
configuration file to set options, and then providing
the path to that file to the pipeline. This makes
it easy to assemble jobs and to keep a history of what
was done, especially if you use our job manager.

Config files are valid YAML files, although the interpreter
we use to read them only processes a subset of YAML features.
[Learn more about YAML on the internet](https://www.google.com/search?q=yaml+basics), 
or just proceed, it is intuitive and easy.

### Config file templates

To get a template to help you write your config file use:

```bash
mdi <pipelineName> template --help
mdi <pipelineName> template -a -c
```

In general, the syntax is:

```yml
# <data>.yml
---
pipeline: [suiteName/]pipelineName[:suiteVersion]
variables:
    VAR_NAME: value
pipelineAction:
    optionFamily:
        optionName1: $VAR_NAME # a single keyed option value
        optionName2:
            - valueA # an array of option values, executed in parallel
            - valueB
execute:
    - pipelineAction
```

The 'suiteName' and 'version' components of the pipeline declaration are optional, however, including suiteName can improve clarity and ensure that
you are always using the tool you intend. If provided, the version designation should be either:
- a suite release tag of the form 'v0.0.0'
- 'latest', to use the most recent release tag [the default]
- 'pre-release', to use the development code at the tip of the main branch
- for developers, the name of a code branch in the git repository

As a convenience for when you get tired of have many files
with the same option values (e.g., a shared data directory), you may
also create a file called 'pipeline.yml' or '\<pipelineName\>.yml'
in the same directory as '\<data\>.yml'. Options will be read
from 'pipeline.yml' first, then '\<data\>.yml', then finally
from any values you specify on the command line, with the last 
value that is read taking precedence, i.e., options specified on the 
command line have the highest precedence.

### Common workflow actions and options

Each pipeline supports parallel processing via options '--n-cpu' and '--ram-per-cpu'.

Option '--dry-run' allows a test of the action and options configuration prior to actual execution.

Many pipelines are designed to have staged execution with success monitoring. Commands 
'status' and 'rollback' provide support for monitoring and manipulating pipeline
execution status, repeating steps, etc. Pipelines that fail part-way through can be 
re-launched from the failure point by re-calling the initial action.




## Stage 1 versioning

### Pipeline versions

Individual pipeline versioning is optional but recommended as it will
help users to confidently access legacy versions of your code to analyze 
their data according to some previous standard, e.g., to ensure consistency 
between older and newer data sets.

Declaring pipeline versions is simple: just add a proper semantic version
declaration to pipeline.yml and update it prior to committing new code. 
It is not necessary to create Git tags for pipeline versions.

```yml
# pipelines/<pipeline>/pipeline.yml
pipeline:
    name: myPipeline
    description: "Description of myPipeline"
    version: v0.0.0
```

### External suite versions

If your pipeline uses code modules from external tool suites, you may
wish to specify the required versions of those external suites.
This is useful if you don't wish to adjust your pipeline to account for a
breaking change made in an external tool suite.  Declare such version
requirements as follows, replacing 'suiteName' with the name of the
external tool suite.

```yml
# pipelines/<pipeline>/pipeline.yml
suiteVersions:
    suiteName: v0.0.0
```

If you do not provide a version for an external tool suite,
the latest version of that suite will be used.

If you only use pipeline code from within your own tool suite, the 
suiteVersions dictionary can be omitted from pipeline.yml.
