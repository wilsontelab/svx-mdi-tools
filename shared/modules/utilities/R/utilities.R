
#----------------------------------------------------------------------
# miscellaneous generic support tools and functions
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# vector tools
#----------------------------------------------------------------------
collapseVector <- function(v, n) { # sum every n adjacent elements of a vector
    cv <- unname(tapply(v, (seq_along(v) - 1) %/% n, sum))
    tailLength <- length(v) %% n # number of input elements summed into incomplete last output element    
    if(tailLength != 0){
        cvLength <- length(cv) # expand incomplete last element to same scale as all others
        cv[cvLength] <- cv[cvLength] * n / tailLength          
    }
    cv
}

#----------------------------------------------------------------------
# data frame tools
#----------------------------------------------------------------------

# reduce a data frame to unique rows based on queried columns
uniqueRows <- function(df, cols) df[!duplicated(df[cols]), ] 

#----------------------------------------------------------------
# string functions
#----------------------------------------------------------------

# convert first character in word to upper case
ucFirst <- function(y) { 
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1, 1)), substring(c, 2),
      sep = "", collapse = " ")
}

#----------------------------------------------------------------
# miscellaneous functions
#----------------------------------------------------------------

# shortcut for the opposite of %in%
`%notin%` <- Negate(`%in%`)

# remove the first element of an vector
# return the altered vector (NOT the shifted value)
shift <- function(v) if(length(v) <= 1) c() else v[2:length(v)]
