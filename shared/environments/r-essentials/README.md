
The packages contained in 'r-essentials' can be difficult to
find and are not well documented.  The following are known
to be present:

- jsonlite
- data.table

To check whether a package is available in a
conda environment, temporarily use these lines:

```bash
# Workflow.sh
Rscript -e 'library(foo); message("package foo OK\n")'
exit 1
```

and run:

````bash
mdi <pipeline> <action> <data.yml>
```

If the package is OK, add it to the list above. If
you get an error that R can't find the package, you
will need to add it to your conda dependencies.
