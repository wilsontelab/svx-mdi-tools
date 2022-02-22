---
published: false
---

## Code modules

**\<suite\>/shared/modules** defines chunks of execution code that are relevant 
to multiple pipelines. There are two types of modules - action modules 
and step modules. Each of those module types might also be used from 
an external pipelines suite.

### Action modules

Action modules are fully encapsulated actions that can be added to a 
pipeline as follows:

```yml
# pipeline.yml
actions:
    actionName: # replace 'actionName' with the name of your action
        order: 1
        thread: name # optional
        module: example/path-0.1 # relative path within 'shared/modules'
```

Thus, you assign a name to the module that is most informative in
the context of your pipeline. You can also declare an 'order' and
'thread' appropriate to your pipeline, as these are never set by
the module. Do not set any of the other typical config values
associated with actions (e.g. optionFamilies), they are defined
by the module.

When writing an action module, use this format:

```yml
# shared/modules/example/path-0.1/module.yml
action:
    optionFamilies:
        - base
        - ...
    condaFamilies: 
        - base
        - ...        
    resources:
        required:
            total-ram: 4G
        recommended: 
            n-cpu: 1
            ram-per-cpu: 4G   
    description: "generic description of the module"
#optionFamilies: # similar to pipeline.yml if needed for families defined by the module
#condaFamilies: 
```

### Step modules

Step modules provide reusable code files that can be called
within a pipeline action that you configured yourself. Use of code in this 
way is more varied and flexible and hard to completely exemplify. One use
case would be a module that provides a Snakefile useful to multiple pipelines:

```bash
# Workflow.sh
checkWorkflowStep 1 exampleName example/path
snakemake --snakefile $MODULES_DIR/example/path-0.1/Snakefile
finishWorkflowStep
```

The environment variable $MODULES_DIR is available to all running pipelines 
to provide easy access to module files.

### External modules

A pipeline may also use a module of either type from a different pipelines suite,
which of course must also be installed into the working MDI directory.

To use an external action module in pipeline.yml, the syntax is:

```yml
# pipeline.yml
actions:
    actionName:
        module: <suite>//example/path-0.1 
```

To use an external step module file in your pipeline scripts, the syntax is:

```bash
# Workflow.sh
$SUITES_DIR/<suite>/shared/modules/example/path-0.1/<file>
```

Care must be taken that file paths within external module scripts respect the fact 
that many environment variables, e.g., $MODULES_DIR, will be defined 
relative to the calling pipelines suite, not the external suite. Thus, the 
author of the external module must generally have anticipated it would be
used by other pipelines by ensuring that all internal file paths will resolve
properly.
