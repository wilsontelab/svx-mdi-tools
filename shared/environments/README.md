---
published: false
---

## Conda runtime environments

**\<suite\>/shared/environments** carries Anaconda/Miniconda run-time environment
configuration files that may be called by pipeline configs. Create
one subfolder for each environment type, which will carry potentially multiple
versioned YAML files for that environment, to allow legacy environments to
be preserved even when current program versions advance.

Load an environment as follows:

```yml
# pipeline.yml
actions:
    actionName: # replace 'actionName' with the name of your action
        condaFamilies:
            - base-0.1 
            - my-conda
condaFamilies:
    my-conda: # defines the 'my-conda' environment component
        channels: # optional, can often be omitted
            - abc
        dependencies: # load specific programs or versions
            - xyz=1.16.3
```

In the example above, 'base-0.1' must exist as a shared component. 'my-conda'
might be fully private to the pipeline, or could also be a shared environment
for which the author needs to override a dependency version, etc.

The total set of condaFamilies are aggregated to create the final runtime
environment, i.e., installed programs, available to the action. The order 
in which the named conda families are listed under
"actions: actionName: condaFamilies:" are the order they are
loaded, i.e. the last one has highest precedence (e.g., to override
to a specific program version).

### Versioning

Conda family version numbers should be incremented whenever any of the
program versions within it are updated, or when the channels list
changes. It is not necessary to increment the family version when
a new program dependency (i.e., one that was previously omitted) is
added. That program may get added to the environment of existing
pipelines, but shouldn't interfere with other programs. Moreover,
if a pipeline had already declared that program as a dependency,
it's version will take precedence over that from the family.

R and Python family versions (and others for which it makes sense)
should match the major.minor version of R or Python, respectively. Other
families are incremented from 0.1 with each update.

If a version is not specified (e.g., 'base'), the latest version
is used. However, it is recommended to always declare
explicit versions in all pipeline configurations.

All older, i.e., now outdated, environment files must always be maintained
to allow users to run legacy versions of all code.

Finally, when a new conda family version is used by a pipeline,
it is important to remember to advance the minor version of the pipeline
and its parent tool suite to reflect the new set of program dependencies. 
Among other things, this ensures that container versions will be updated to 
include the new conda environments.

### Available conda families

The following environment families are provided by the template
as they are typical for many data analysis needs and support the 
demo, i.e., '_template', pipeline.

- **base** establishes common conda channels and data tools
- **r-essentials** adds R and a set of commonly used data packages
- **python-essentials** adds Python and a set of commonly used data packages
