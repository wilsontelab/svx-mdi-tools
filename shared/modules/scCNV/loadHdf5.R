#=====================================================================================
# extract required data elements from CellRanger-compatible hdf5 file
# discard some Cell Ranger outputs, e.g., normalization and cell clusters, if present
#-------------------------------------------------------------------------------------
message("loading bin data from HDF5 file")

# find the hdf5 file, can be nested in a subfolder of inputDir/inputName
inputDir <- file.path(env$INPUT_DIR, env$INPUT_NAME)
h5FileName <- "cnv_data.h5"
h5Files <- list.files(
    path = inputDir, 
    pattern = h5FileName, 
    full.names = TRUE, 
    recursive = TRUE
)
h5File <- h5Files[1]
if(h5File == "") stop(paste("file", h5FileName, "not found in directory", inputDir))

# extract the required data elements
h5 = H5Fopen(h5File)
constants <- h5read(h5, "constants", bit64conversion = "int")
genome_tracks <- h5read(h5, "genome_tracks", bit64conversion = "int")
genome_tracks$is_mappable <- NULL
metadata <- h5read(h5, "metadata", bit64conversion = "int")
per_cell_summary_metrics <- h5read(h5, "per_cell_summary_metrics", bit64conversion = "int")
raw_counts <- h5read(h5, "raw_counts", bit64conversion = "int")
H5Fclose(h5)
rm(h5FileName, h5Files, h5File, h5)
#=====================================================================================
