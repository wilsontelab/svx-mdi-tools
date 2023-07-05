#-------------------------------------------------------------------
# define available heat map colors
#-------------------------------------------------------------------
# define available colors
cnvHues <- list( 
    red   = list(
        max = 0.8, # max intensity prevents excessive image brightness
        value = rgb(0.8, 0, 0, maxColorValue = 1),
        index = 1
    ),  
    green = list(
        max = 0.65,
        value = rgb(0, 0.65, 0, maxColorValue = 1),
        index = 2
    ), 
    blue  = list(
        max = 0.9,
        value = rgb(0, 0, 0.9, maxColorValue = 1),
        index = 3
    )
)
cnvLossColor = "blue" # consistent with prior work, e.g. blue = deletion
cnvGainColor = "red"
cnvHighlightColor = "green" # yields yellow for red+green, cyan for blue+green

#-------------------------------------------------------------------
# save intensity data as png heat map using ImageR
# much faster than heat map plotting in R (e.g. by rect or other)
#-------------------------------------------------------------------
# convert an intensity matrix to an array of RGB colors
parseHeatMapIntensity <- function(intensity, triangle, mirror){ # step 1: apply triangle, mirror and inversion
    intensity[is.na(intensity)]  <- 0 # ImageR requires actual values
    if(triangle) {
        intensity <- as.matrix(imager::imrotate(imager::as.cimg(intensity), -45)) # -45 is the angle that works!
        nc <- ncol(intensity)
        intensity <- intensity[, ceiling(nc / 2):nc]
    } else if(mirror) {
        intensity <- mirrorPyramid(intensity) # make rectangular view be symmetric
    }    
    t(intensity) # transpose for ImageR coordinate consistency
}
parseHeatMapHue <- function(intensity, hue){ # step 2: create the colors array
    intensity <- 1 - intensity # for color, 0 is max intensity
    max <- cnvHues[[hue]]$max
    hue_ <- max + (1 - max) * intensity
    switch(hue,
        "red"   = abind::abind(hue_, intensity, intensity, along = 3), # abind to merge matrices
        "green" = abind::abind(intensity, hue_, intensity, along = 3),
        "blue"  = abind::abind(intensity, intensity, hue_, along = 3)
    )    
}
printHeatMapImage <- function(img, file, xScaleFactor, yScaleFactor, transpose, decorate, ...){
    if(transpose){
        tmp <- xScaleFactor
        xScaleFactor <- yScaleFactor
        yScaleFactor <- tmp
        img <- aperm(img, c(2, 1, 3))
    }
    cimg <- imager::imrotate(suppressWarnings(imager::as.cimg(img)), -90) # rotate so 1,1 is bottom left of rectangle, or left of triangle
    if(xScaleFactor != 1 | yScaleFactor != 1) { # expand/contract the image when instructed
        dim <- dim(cimg)
        cimg <- imager::resize(cimg, size_x = dim[1] * xScaleFactor, size_y = dim[2] * yScaleFactor)
    }
    if(!is.null(decorate)) cimg <- decorate(cimg, xScaleFactor, yScaleFactor, ...) # decorate applied AFTER expansion
    imager::save.image(cimg, file)
}
# save intensity matrix [0,1=high] as a single-color heat map png
# triangle = TRUE yields a rotated triangle corresponding to the upper triangle of identity matrix
# mirror = TRUE makes a rectangular matrix be symmetrix from upper triangle
saveHeatMap_one_color <- function(intensity, hue, file, ..., # ... passed to decorateFN
                                  triangle = FALSE, mirror = TRUE,
                                  xScaleFactor = 1, yScaleFactor = 1, transpose = FALSE,
                                  decorate = NULL){   
    intensity <- parseHeatMapIntensity(intensity, triangle, mirror)
    cimg <- suppressWarnings(imager::as.cimg(parseHeatMapHue(intensity, hue)))
    printHeatMapImage(cimg, file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}
# same thing for two color intensity on a diverging axis (i.e. low color, high color)
# when input low intensity is >0, high should be 0, and vice versa
saveHeatMap_two_color <- function(intensityLow, intensityHigh, highlight, file, ..., # ... passed to decorateFN
                                  triangle = FALSE, mirror = TRUE,
                                  xScaleFactor = 1, yScaleFactor = 1, transpose = FALSE,
                                  decorate = NULL){
    intensityLow  <- parseHeatMapIntensity(intensityLow,  triangle, mirror)
    intensityHigh <- parseHeatMapIntensity(intensityHigh, triangle, mirror)
    highlight     <- parseHeatMapIntensity(highlight,     triangle, mirror)
    colorsLow     <- parseHeatMapHue(intensityLow,  cnvLossColor) 
    colorsHigh    <- parseHeatMapHue(intensityHigh, cnvGainColor)
    for(i in 1:3){ # execute the merge of the two, mutually exclusive, colors
        colorsHigh[,,i] <- ifelse(intensityLow > 0, colorsLow[,,i], colorsHigh[,,i])
    }
    colorsHigh[, , cnvHues[[cnvHighlightColor]]$index] <- pmax(
        highlight[, ],
        colorsHigh[, , cnvHues[[cnvHighlightColor]]$index]
    )
    printHeatMapImage(colorsHigh, file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}
# more complicated pseudo-coloring where low distribution = blue, middle = green, high = red
saveHeatMap_three_color <- function(intensityLow, intensityMid, intensityHigh,
                                    file, ..., # ... passed to decorateFN
                                    triangle = FALSE, mirror = TRUE,
                                    xScaleFactor = 1, yScaleFactor = 1, transpose = FALSE,
                                    decorate = NULL){
    printHeatMapImage(abind::abind(
        parseHeatMapIntensity(intensityHigh, triangle, mirror), # here, do not invert intensities
        parseHeatMapIntensity(intensityMid,  triangle, mirror),
        parseHeatMapIntensity(intensityLow,  triangle, mirror),
        along = 3
    ), file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}

#-------------------------------------------------------------------
# convert Z scores to color intensities in range [0,1]
#-------------------------------------------------------------------
# ramp up color signal only when Z >= 0
zToIntensity_pos <- function(Z, minZ = 1, maxZ = 4){
    Z[Z < 0] <- 0
    deltaZ <- pmin(pmax(Z, minZ), maxZ)
    (deltaZ - minZ) / (maxZ - minZ)  
}
# decay signal for negative Z, full intensity for Z > 0
zToIntensity_neg_decay <- function(Z, maxZ = 3){
    fullOn <- Z >= 0
    deltaZ <- pmin(pmax(-Z, 0), maxZ)
    ifelse(fullOn, 1, 1 - deltaZ / maxZ)
}
# decay signal on both sides of peak, full intensity for Z == 0
zToIntensity_symmetric <- function(Z, maxZ = 3){
    deltaZ <- pmin(abs(Z), maxZ)
    1 - deltaZ / maxZ
}
# 0.5 intensity for Z == 0, 0 for -maxZ, 1 for +maxZ
zToIntensity_midpoint <- function(Z, maxZ = 3){
    pmax(pmin(Z / maxZ / 2 + 0.5, 1), 0)
}
# ramp up signal on both sides of peak, no intensity for Z == 0
zToIntensity_inverted <- function(Z, minZ = 1, maxZ = 4){
    deltaZ <- pmin(pmax(abs(Z), minZ), maxZ)
    (deltaZ - minZ) / (maxZ - minZ)  
}

#-------------------------------------------------------------------
# functions to decorate ImageR images
#-------------------------------------------------------------------
adjustCoord <- function(pos, expand, max){
    pos <- (pos - 1) * expand # thus, all pixels of the region are VISIBLE (box masks adjacent pixels)
    min(max(pos, 1), max)
}
invertY <- function(y, ymax) ymax - y + 1
addRegionBox <- function(cimg, expand, triangle, box, col){ # box=c(x0, y0, x1, y1)
    if(is.null(box)) return(cimg)
    xmax <- dim(cimg)[2]
    ymax <- dim(cimg)[1]
    x0 <-         adjustCoord(box[1], expand, xmax) # remember, x are the columns
    y0 <- invertY(adjustCoord(box[2], expand, ymax), ymax)
    x1 <-         adjustCoord(box[3], expand, xmax) 
    y1 <- invertY(adjustCoord(box[4], expand, ymax), ymax)
    implot(cimg, rect(x0, y0, x1, y1, border=col, lwd=1))
}

#-------------------------------------------------------------------
# functions to resize images
#-------------------------------------------------------------------


expandImg <- function(img, h=1, v=1){
  if(h > 1){
    imgD <- dim(img)
    lrg <- as.cimg(array(0, dim=c(imgD[1]*h, imgD[2], 1, 3)))   
    for(x in 0:(h-1)) lrg[1:imgD[1] * h - x,,,] <- img[,,,]
    img <- lrg
  }
  if(v > 1){
    imgD <- dim(img)
    lrg <- as.cimg(array(0, dim=c(imgD[1], imgD[2]*v, 1, 3)))
    for(x in 0:(v-1)) lrg[,1:imgD[2] * v - x,,] <- img[,,,]
    img <- lrg
  }
  img
}