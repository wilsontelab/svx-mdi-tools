---
published: false
---

## Shared components

MDI Stage 1 Pipelines suites support three types of reusable code 
components that can be shared between all pipelines in the suite,
which are defined in the '**\<suite\>/shared**' folder.

- **environments** = yml config files that create conda environments for job execution
- **modules** = script libraries that provide code for use by running pipelines
- **options** = yml config files that expose option families for job configuration

In all cases, the components are effectively placed inline into the 
pipeline.yml file that configures a specific pipeline.

Please see the environments, modules, and options documentation
for more information.

### Private vs. shared components

Environments and options can be defined privately for a 
pipeline within its 'pipeline.yml' file, used from the shared folder, or both.
The framework first looks to see if a component is present in the shared
folder and loads that configuration first. It then further looks to see 
if there are pipeline-specific definitions for the named component in 
pipeline.yml. If no shared component was found, the definition must 
exist in its entirety in pipeline.yml. If a shared component was found,
any further definitions in pipeline.yml override the shared configuration.

Modules, by their nature, are encapsulated components that are 
always loaded from the 'shared' folder.

### External components

Component sharing can extend beyond a single tool suite, such that
pipeline.yml may also attempt to load an environment, module, or option 
family from a different tool suite, which must also be installed into 
the working MDI directory by setting 'suite_dependencies' in the calling suite's
_config.yml file.

```yml
# _config.yml
suite_dependencies:
    - <git_user>/<suite>
```

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

## Shared component versioning

Similar to pipelines, the version of a shared component is implicitly
derived from the version of its parent suite, i.e., setting the
version of a tool suite always yields the same, specific version of the component. 

If your pipeline requires a specific version of an external tool suite
(and therefore of its components), you may override the default version of 
'latest' at the pipeline level:

```yml
# pipeline.yml
pipeline: ...
suiteVersions: # must come before 'actions'
    suiteName: v0.0.0 # use this version of a suite invoked as 'suite//module', etc. [latest]
actions: ...
```
