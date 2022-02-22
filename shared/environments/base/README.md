
**base-xx.xx.yml** configures the basic runtime environment and shell.

Many pipelines load 'base' first, followed by:
- program families specified by shared conda environments
- any dependencies/versions peculiar to the pipeline (which override the above)
