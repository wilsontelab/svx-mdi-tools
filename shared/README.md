---
published: false
---

## Shared components

MDI Stage 1 Pipelines suites support three types of reusable code 
components that can be shared between all pipelines in the suite,
which are defined in the '**\<suite\>/shared**' folder.

- **environments** are yml config files that create conda environments for job execution
- **modules** are script libraries that provide code for use by running pipelines
- **options** are yml config files that expose option families for job configuration

In all cases, the components are effectively placed inline into the 
pipeline.yml file that configures a specific pipeline.

### Private vs. shared components

Each of the component types listed above can be defined privately for a 
pipeline within its 'pipeline.yml' file, used from the shared folder, or both.

The framework first looks to see if a component is present in the shared
folder and loads that configuration first. It then further looks to see 
if there are pipeline specific definitions for the named component in 
pipeline.yml. If no shared component was found, the definition must 
exist in its entirety in pipeline.yml. If a shared component was found,
any further optional definitions in pipeline.yml will override the shared 
configuration.

Please see the environments, modules and options documentation
for specific syntax examples.

### External components

Component sharing can also extend beyond a single pipelines suite, such that
pipeline.yml may also attempt to load an environment, module, or option family from
a different pipelines suite, which of course must also be installed into 
the working MDI directory.

To use external components in pipeline.yml, prefix the component path
with '\<suite\>//', where \<suite\> is the name of the external suite from
which to load the component.

```yml
# pipeline.yml
actions:
    actionName:
        condaFamilies:
            - <suite>//shared-conda  
        module: <suite>//example/shared-module
        optionFamilies:
            - <suite>//shared-options
```

The pipelines framework will only look for external component files
in _definitive_ suite repositories; we cannot assume that an end 
user will have an active fork of any given external repository.

If your pipeline requires a specific version of the external tool suite,
you may override the default version of 'latest' at the pipeline level 
as follows:

```yml
# pipeline.yml
pipeline: ...
suiteVersions: # must come before 'actions'
    suiteName: v0.0.0 # use this version of a suite invoked as 'suite//module', etc. [latest]
actions: ...
```
