
#----------------------------------------------------------------------
# data_sources.R defines server-specific extensions to ../data_sources_session.R
#----------------------------------------------------------------------

# derived SV column data
addPlotValues <- function(dt, posPrefix, sample){
    reportProgress('addPlotValues')
    
    getCTIs <- function(side, posPrefix){ # indices of the capture targets
        chromCol <- paste("CHROM_",  side, sep="")
        posCol   <- paste(posPrefix, side, sep="")
        cts <- captureTargets[1:nCaptureTargets,]
        regionMatches <- mapply(function(chrom, pos){
            cts$chrom       == chrom &
            cts$paddedStart <= pos & # NB: padded
            cts$paddedEnd   >= pos
        }, dt[[chromCol]], dt[[posCol]])
        ctI <- rowSums(sapply(1:nCaptureTargets, function(ctI) ifelse(regionMatches[ctI,], ctI, 0)))
        ifelse(ctI==0, nCaptureTargets + 1, ctI)
    }
    getIs <- function(side, posPrefix, offset=TRUE){ # plot indices
        ctICol <- paste("ct",      side, sep="")
        posCol <- paste(posPrefix, side, sep="")
        ifelse(
            is.na(dt[[ctICol]]),
            NA,
            dt[[posCol]] - captureTargets[dt[[ctICol]], 'start'] + projectInfo$REGION_PADDING +
                if(offset) captureTargets[dt[[ctICol]], 'iOffset'] else 0
        )
    }    

    # get the _padded_ regions into which the ends fall
    dt[,':='(
        ct1 = getCTIs(1, posPrefix),
        ct2 = getCTIs(2, posPrefix)     
    )]

    # plot indices that increment over all regions (e.g. for triangle plot)
    dt[,':='(
        i1 = getIs(1, posPrefix), 
        i2 = getIs(2, posPrefix),
        
        # plot indices relative to the parent region (e.g. for linear plot)
        j1 = getIs(1, posPrefix, FALSE), 
        j2 = getIs(2, posPrefix, FALSE)
    )]

    # set values that aggregate over both SV ends
    dt[, ':='(
        size   = abs(i2 - i1), # as measured in plot coordinates NOT bp, so e.g. tt translocations DO have a size
        center = pmin(i1, i2) + abs(i2 - i1) / 2,
        
        # calculate the fraction of ends across all molecules that were shared with proper molecules
        # a few molecules are coming through find with dt$N_TOTAL>0 but (dt$N_GAPS + dt$N_SPLITS)==0
        # likely because of failed alignment to the junction?
        N_GAP_SPLIT = N_GAPS + N_SPLITS,
        FRAC_SHARED_PROPER = round(N_SHARED_PROPER / (N_GAPS + N_SPLITS), 2),
        
        # add a sample name column and return
        SAMPLE = sample
    )]

    # fix an oddity wherein sequence values 'NA' (a dinucleotide) are made NA (not avail) in R read.table
    dt[is.na(JXN_BASES) & !is.na(JXN_SEQ), JXN_BASES := 'NA']

    dt   
}

# specific SV filters
addSpecificFilters <- function(boo, sv){
    
    if(input$captureTarget != '-'){ # SVs with both ends mapped to a specific capture target
        region1 <- captureTargets[sv$ct1, 'region']
        region2 <- captureTargets[sv$ct2, 'region']
        boo <- boo & !is.na(sv$ct1) & !is.na(sv$ct2) & region1 == region2 & region1 == input$captureTarget
    }
    for(tgtType in allPairTypes){ # where the ends of an SV fell 
        if(tgtType %notin% input$targetTypeFilter) boo <- boo & sv$TARGET_CLASS != tgtType
    }        
    if(input$sequencedFilter != 'all'){ # whether we sequenced across the junction
        isSequenced <- isSequencedJunction(sv)
        boo <- boo & if(input$sequencedFilter == 'yes') isSequenced else !isSequenced
    }
    
    boo
}
