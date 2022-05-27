---
published: false
---

## Code modules

**\<suite\>/shared/modules** defines chunks of execution code relevant 
to multiple pipelines. There are two types of modules - action modules 
and step modules. Each of module type might also be used from 
an external pipelines suite.

## Action modules

Action modules are fully encapsulated actions that can be added to a 
pipeline as follows:

```yml
# pipeline.yml
actions:
    actionName: # replace 'actionName' with the name of your action
        order: 1
        thread: name # optional
        module: example/path # relative folder path within 'shared/modules'
```

Thus, you assign a name to the module that is most informative in
the context of your pipeline. You can also declare an 'order' and
'thread' appropriate to your pipeline, as these are never set by
the module. Do not set any of the other typical config values
associated with actions (e.g., optionFamilies), they are defined
by the module. Similarly, _global definitions are not applied
to shared action modules, which are intended to be self-contained.

When writing an action module, use this format:

```yml
# shared/modules/example/path/module.yml
---
version: v0.0.0 # optional, for internal tracking
action: # required
    condaFamilies: 
        - shared-family
        - inline-family  
    optionFamilies:
        - shared-family
        - inline-family
    resources:
        required:
            total-ram: 4G
        recommended: 
            n-cpu: 1
            ram-per-cpu: 4G   
    job-manager:
        recommended:
            time-limit: 48:00:00   
    description: "generic description of the module's action"
condaFamilies: # if needed for declarations above
    inline-family: ...
optionFamilies:
    inline-family: ...
```

All condaFamilies and optionFamilies in module.yml are interpreted 
relative to the module's suite, not the calling suite, so that
suite developers have complete control over a module's action
even when it is called by another suite as an external module.

If inline component families are specified within module.yml,
they are appended at the end of the working pipeline.yml file 
during execution and therefore override families of the same name 
in the calling suite. It is the job of the calling pipeline to  
manage any collisions in family names between different actions. 

The calling pipeline and/or suite, not the module's tool suite, 
are responsible for _building_ the required conda environment or 
Singularity container using the definitions provided by an action module.

## Step modules

Step modules provide reusable code files that can be called
within a pipeline action that you configured yourself. Use of code in this 
way is more varied and flexible and hard to completely exemplify. One use
case would be a module that provides a Snakefile useful to multiple pipelines:

```bash
# Workflow.sh
snakemake --snakefile $MODULES_DIR/example/path/Snakefile
```

The environment variable $MODULES_DIR is available to all running pipelines 
to provide easy access to module files.

## External modules

A pipeline may also use a module of either type from a different pipelines suite, 
which must also be installed into the working MDI directory by setting 
'suite_dependencies' in the calling suite's _config.yml file.

To use an external action module in pipeline.yml, the syntax is:

```yml
# pipeline.yml
actions:
    actionName:
        module: <suite>//example/path
```

To use an external step module file in your pipeline scripts, the syntax is:

```bash
# Workflow.sh
$SUITES_DIR/<suite>/shared/modules/example/path/<file>
```

Care must be taken that file paths within external module scripts respect the fact 
that many environment variables, e.g., $MODULES_DIR, will be defined 
relative to the calling pipelines suite, not the external suite. Thus, the 
author of the external module must generally have anticipated it would be
used by other suites by ensuring that all file paths resolve properly.

## Module versioning

The version of a shared module is implicitly derived from the version of 
its parent suite, i.e., setting the version of a tool suite always yields 
the same, specific version of a module. 

For external shared modules, the module version can be set by requiring 
a specific version of the external suite at the pipeline level:

```yml
# pipeline.yml
suiteVersions: 
    suiteName: v0.0.0 
```

For internal shared modules, if two pipelines in your tool suite require different 
versions of a similar module they must have different names so that
they can be called differently by the two pipelines.
Alternatively, one pipeline could recreate the module in its pipeline.yml file.
