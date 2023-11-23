# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
filteredSvsTableUI <- function(ns, width){
    box(
        width = width,
        bufferedTableUI(ns('svsTable'))
    )  
}
filteredSvsTableServer <- function(id, input, svs, matchThreshold = NULL){
    bufferedTableServer(
        id = 'svsTable',
        parentId = id,
        parentInput = input,
        selection = 'single',
        tableData = reactive({
            svs <- svs()
            if(is.null(svs)) return( data.frame(message = "analysis not available for these data") )     
            svs[, .(
                svId = SV_ID,
                #---------------
                type = svx_jxnType_codeToX(edgeType, "name"),
                class = TARGET_CLASS,
                target = TARGET_REGION,
                #---------------
                nSmp = N_SAMPLES,
                nTot = N_TOTAL,
                nSplit = N_SPLITS,            
                nGap = N_GAPS,
                nClip = N_OUTER_CLIPS,
                nDpx = N_DUPLEX_GS,
                fPrp = round(SHARED_PROPER / 2, 2),
                #---------------
                count = STRAND_COUNT_GS,
                count1 = STRAND_COUNT1,
                count2 = STRAND_COUNT2,
                #---------------
                mapQ = MAX_MAPQ,
                # MATCH_TYPE = MATCH_TYPE,
                snv = if(is.null(matchThreshold)) "-" else ifelse(MATCH_TYPE >= matchThreshold(), "*", "-"),
                #---------------
                chr1 = CHROM_1,
                pos1 = POS_1,
                chr2 = CHROM_2,
                pos2 = POS_2,
                #---------------
                svSize = SV_SIZE,
                uHom = MICROHOM_LEN,
                jxnSeq = JXN_BASES
            )]
        })
        # 'SHARED_PROPER' =  'double',
        # 'SHARED_PROPER_GS'  =  'double',
        # #---------------
        # 'N_AMBIGUOUS'   =  'integer',
        # 'N_DOWNSAMPLED' =  'integer',
        # 'N_COLLAPSED'   =  'integer',
        #     lm[, ':='(
        #         FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
        #         Failed = failed,
        #         nReadPairs = commify(nReadPairs),
        #         nSourceMolecules = commify(nSourceMolecules),
        #         onTargetCoverage = commify(round(onTargetCoverage, 0)),
        #         offTargetCoverage = round(offTargetCoverage, 3),
        #         enrichment = commify(round(enrichment, 0)),
        #         efficiency = round(efficiency, 3)
        #     )]
        # editBoxes = list( # handle the libraryFailed checkbox
        #     libraryFailed = list(
        #         type = 'checkbox',
        #         boxColumn = 1,
        #         rawColumn = 2,
        #         handler = function(checked){ # enter the new failure value into our outcomes
        #             ss <- sampleSet$input$sampleSet
        #             failed <- outcomes[[ss]]
        #             failed[checked$selectedRow] <- checked$newValue
        #             outcomes[[ss]] <- failed # must replace the entire value for outcomes invalidation
        #             checked
        #         }
        #     )
        # )
    )
}
