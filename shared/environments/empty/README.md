**empty-0.1.yml** configures an empty runtime environment and shell.

Although uncommon, sometimes you might not wish to install any
specific program dependencies in the conda environment. In such cases, use:

```yml
actions: 
    actionName:
        condaFamilies:
            - empty-0.1
```
