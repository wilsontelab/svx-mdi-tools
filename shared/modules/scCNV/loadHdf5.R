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

# # if possible, replace the svx pipeline's gc and mappability tracks with CellRanger's
# # unfortunately, the bin counts do not match between Cell Ranger and our genome builds
# genomeTracksDir <- file.path(env$MDI_DIR, "resources", "scCNV")
# genomeTracksGroup <- 'genome_tracks'
# genomeTrackFiles <- list(
#     GRCh38 = file.path(genomeTracksDir, paste("GRCh38", genomeTracksGroup, "h5", sep = ".")),
#     GRCm38 = file.path(genomeTracksDir, paste("GRCm38", genomeTracksGroup, "h5", sep = "."))
# )
# genomeTrackFiles$hg38 = genomeTrackFiles$GRCh38 # support both name types
# genomeTrackFiles$mm10 = genomeTrackFiles$GRCm38
# h5File <- genomeTrackFiles[[metadata$assembly]]
# if(is.null(genome_tracks$n_fraction)){ # TRUE for svx, this track not added during 'extract'
#     if(!is.null(h5File) && file.exists(h5File)){
#         message("overriding genome_tracks to Cell Ranger values")
#         message(h5File)
#         h5 = H5Fopen(h5File)
#         genome_tracks <- h5read(h5, "genome_tracks", bit64conversion = "int")
#         genome_tracks$is_mappable <- NULL
#         H5Fclose(h5)
#         rm(h5)
#     }
# } else { # only needed to export a CellRanger gene track file; uncomment as needed
#     if(!file.exists(h5File)){
#         dir.create(dirname(h5File), showWarnings = FALSE, recursive = TRUE)
#         h5createFile(h5File)
#         h5write(get(genomeTracksGroup), h5File, genomeTracksGroup)
#     }
# }
#=====================================================================================
