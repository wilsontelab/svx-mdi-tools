---
published: false
---

## User options for configuring pipeline calls

**\<suite\>/shared/options** defines typical sets of pipeline configuration
options, i.e., values that may/must be specified by end users at the command
line or in their data.yml file. 

Options are organized into optionFamilies for clarity of presentation
and easier reading of configuration files. Option families can be invoked 
by pipeline configuration files as follows:

```yml
# pipeline.yml
actions:
    actionName: # replace 'actionName' with the name of your action
        optionFamilies:
            - shared-options # a shared option family
            - my-options
optionFamilies:
    my-options: # defines the 'my-options' option family
        order: 1
        options:
            my-option-name:
                order: 1
                short: x
                type: string
                required: true # set to false if a default is provided
                default: null
                description: "short description of the option's effect"
```

In the example above, 'shared-options' must exist as a shared component, 
i.e., file 'shared/options/shared-options.yml' must exist. 
'my-options' might be fully private to the pipeline, or could also be a 
shared option family for which the author needs to override some configuration
detail such as changing 'required' from true to false while providing
a 'default' value.

The 'order' key:value pairs allows you to control the order the families 
and their options are listed on help screens.

## Creating shared option families

Shared option families are defined in YAML configuration files as follows:

```yml
# shared/options/NAME.yml = a single option family called NAME
order: 1 # optional
options: # required
    option-1: ...
    option-2: ...
```
