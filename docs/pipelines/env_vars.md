---
title: Environment Variables
parent: Stage 1 Pipelines
has_children: false
nav_order: 5
---

## Environment variables

The job launcher sets a number of environment variables that are available
to your action script. 

First, all pipeline options are available, where an option named "--abc-def" 
sets environment variable "ABC_DEF", e.g., '--n-cpu' becomes N_CPU, etc.

Additionally, the launcher sets the following derivative environment variables,
which are useful for locating files and for other purposes:

| variable name | value | description |
|---------------|---------------|-------------|
| **TASK_DIR**          | $OUTPUT_DIR/$DATA_NAME | where output files should be placed |
| **DATA_FILE_PREFIX**  | $TASK_DIR/$DATA_NAME   | prefix to use for output file names |
| **PLOTS_DIR**         | $TASK_DIR/plots        | where output plots should be placed |
| **PLOT_PREFIX**       | $PLOTS_DIR/$DATA_NAME  | prefix to use for output plot names |
| **SUITE_NAME**        | | the name of the tool suite that carries the running pipeline |
| **PIPELINE_NAME**     | | the name of the running pipeline |
| **PIPELINE_ACTION**   | | the name of the running pipeline action being applied to $DATA_NAME |
| **TASK_PIPELINE_DIR** | $TASK_DIR/$PIPELINE_NAME            | status, code and log files specific to the running pipeline and task |
| **TASK_ACTION_DIR**   | $TASK_PIPELINE_DIR/$PIPELINE_ACTION | code and log files specific to the running action and task |
| **SUITES_DIR**        | $TASK_ACTION_DIR/suites                | the directory where working versions of accessible MDI code suites are found |
| **SUITE_DIR**         | $SUITES_DIR/$SUITE_NAME                | the working root directory of tool suite $SUITE_NAME |
| **PIPELINE_DIR**      | $SUITE_DIR/pipelines/$PIPELINE_NAME    | the root directory of pipeline $PIPELINE_NAME, which contains 'pipeline.yml' |
| **ACTION_DIR**        | $PIPELINE_DIR/$PIPELINE_ACTION         | the directory that contains the scripts for action $PIPELINE_ACTION, including 'Workflow.sh' |
| **SCRIPT_DIR**        | $ACTION_DIR                            | legacy synonym for above |
| **ACTION_SCRIPT**     |                                        | the primary action script, usually '$ACTION_DIR/Workflow.sh' (unless overridden) |
| **SCRIPT_TARGET**     | $ACTION_SCRIPT                         | legacy synonym for above |
| **MODULES_DIR**       | $SUITE_DIR/shared/modules              | the directory that contains all shared code modules in suite $SUITE_NAME |
| **LOGS_DIR**          | $TASK_ACTION_DIR/logs | where log files should be placed |
| **LOG_FILE_PREFIX**   | $LOGS_DIR/$DATA_NAME  | prefix to use for log file names |
| **TASK_LOG_FILE**     | $LOG_FILE_PREFIX.$PIPELINE_NAME.$PIPELINE_ACTION.task.log | the main log file for a running task; starts with job config YAML |
| **RAM_PER_CPU_INT**   | int($RAM_PER_CPU)                     | e.g., 1M becomes 1000000 |
| **TOTAL_RAM_INT**     | $RAM_PER_CPU_INT * $N_CPU             | RAM available to entire job, e.g., 1000000 * 4 = 4000000 |
| **TOTAL_RAM**         | $TOTAL_RAM_INT (as string)            | e.g., 4000000 becomes 4M |
