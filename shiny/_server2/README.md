
.../svCapture/_server

contains all code required to run a series of R Shiny
data visualization and retrieval servers.

Targets     filter SVs by their properties
            make triangle and other plots
            see and download tabular SV lists
            see molecule support and junctions of SVs

Peaks       make linear plots of capture targets
            to reveal SV peaks surrounding targets

Distribution   make distribution plots of various SV attributes
               while filtering for other attribute values

There are two main ways to launch a server:

1) ./_server/_shiny/run_app.sh
   with appropriate environment variables set
   (see run_app.sh for details)
   
2) source('./_server/<subserver_name>/local.R')
   from R, to run a local server instance
   after editing local.R with proper data paths

