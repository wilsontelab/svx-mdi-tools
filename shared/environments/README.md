---
published: false
---

## Conda runtime environments

**\<suite\>/shared/environments** carries Anaconda/Miniconda run-time environment
configuration files that may be called by pipeline configs as follows:

```yml
# pipeline.yml
actions:
    actionName: # replace 'actionName' with the name of your action
        condaFamilies:
            - base
            - my-conda
condaFamilies:
    my-conda: # defines the 'my-conda' environment component
        channels: # optional, can often be omitted
            - abc
        dependencies: # load specific programs or versions
            - xyz=1.16.3
```

In the example above, 'base' must exist as a shared component, 
i.e., file 'shared/environments/base.yml' must exist. 'my-conda' might 
be fully private to the pipeline, or could also be a shared environment
for which the author needs to override a dependency version, etc.

The total set of condaFamilies is aggregated to create the final runtime
environment, i.e., installed programs, that are available to the action. 
The order in which the named conda families are listed under
"actions: actionName: condaFamilies:" is the order they are
loaded, i.e., the last one has highest precedence (e.g., to override
to a specific program version). However, conda entries for an action
in pipeline.yml will override entries in any shared environment file. 

## Creating shared environments

Shared environments are defined in YAML configuration files in 
'shared/environments' using the following syntax, where 
the name of the file is the name of the conda family. 

```yml
# shared/environments/NAME.yml = a single conda family called NAME
channels: ... # often omitted if defined upstream
dependencies: ...
```

## Available conda families

The following environment families are provided by the MDI suite template
as they are typical for many data analysis needs or support the 
demo, i.e., '_template', pipeline.

- **base** = establishes common conda channels and data tools
- **empty** = establishes common conda channels only
- **r-4.1** = adds R and a set of commonly used data packages

## Environment versioning

The version of a shared conda family is implicitly derived from the version of 
its parent suite, i.e., setting the version of a tool suite always yields 
the same, specific version of an environment config. 

When a conda family changes the version of a program in its environment,
it is important to advance the minor version of any pipelines that use it 
as well as the parent tool suite, to reflect the new set of program 
dependencies. Among other things, this ensures that container versions 
will be updated to reflect the new conda environments.

For external shared conda families, the environment version can be set by requiring 
a specific version of the external suite at the pipeline level:

```yml
# pipeline.yml
suiteVersions: 
    suiteName: v0.0.0 
```

For internal shared environments, if two pipelines in your tool suite require different 
versions of a similar conda family they must have different names so that
they can be called differently by the two pipelines.
Alternatively, one pipeline could override part of the conda family in its pipeline.yml file,
e.g., setting a specific version of one program in the environment.

For R, Python, and similar versioned languages or platforms, 
it is recommended to include the language version in the name of the 
conda family, e.g. 'R-4.1'. That way, a suite can easily maintain
environments for different language versions and the version in use
is clear to everyone.
