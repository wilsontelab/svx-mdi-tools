#----------------------------------------------------------------------
# define generic functions, i.e., methods, available to S3 classes,
# if the class declares method.class <- function()
#----------------------------------------------------------------------
# new generic declarations take the form:
#     genericName <- function(x, ...) {
#         UseMethod("genericName", x)
#     }
#----------------------------------------------------------------------
