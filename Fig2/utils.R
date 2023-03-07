scale_fill_MOCHA<- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('#F26E65', 
                            '#31B34A', 
                            '#89ACDA'
                            ), 
                          c('HOMER','MACS2',
                            'MOCHA'  
                            )),
    )
}


# scale_fill_MOCHA<- function(...){
#     ggplot2:::manual_scale(
#         'fill', 
#         values = setNames(c('#89ACDA', 
#                             '#F26E65', 
#                            '#31B34A'
#                             ), 
#                           c('Mocha','ArchR',
#                             'Signac'
#                             )),
#     )
# }


scale_col_MOCHA<- function(...){
    ggplot2:::manual_scale(
        'col', 
        values = setNames(c('#F26E65', 
                            '#31B34A', 
                            '#89ACDA'
                            ), 
                          c('HOMER','MACS2',
                            'MOCHA'  
                            )),
    )
}

scale_fill_cell<- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('#FF2C8C', 
                            '#D88000', 
                            '#A500A5',
                            '#00A1BC'
                            ), 
                          c('Monocytes','B',
                            'T','Myeloid'
                            )),
    )
}





sample_specific_rowRanges <- function(tileVec){
    openTiles <- tileVec[! is.na(tileVec) & tileVec==T]
    sampleGranges <- StringsToGRanges(names(openTiles))
    sampleGranges$tileID = names(openTiles)
    sampleGranges = as.data.table(sampleGranges)
    return(sampleGranges)
}


### reproducible_peaks by method 
getConsensusTiles <- function(tile_list, reproLimit=0.2){
            full_sc <- rbindlist(tile_list)
            tile_counts <- full_sc[, list(Prop = .N/length(tile_list)), 
                                   by=tileID]
            consensusTiles <- tile_counts[ Prop > reproLimit]

            consensusTiles$seqnames <- gsub(':[:0-9]*[-0-9]*','', consensusTiles$tileID)
            consensusTiles$range <- gsub('chr[0-9XY]*:','', consensusTiles$tileID)
            consensusTiles$start <- as.numeric(gsub('-[-0-9]+','', consensusTiles$range))
            consensusTiles$end <- as.numeric(gsub('[-0-9]+-','', consensusTiles$range))

            consensusTiles_gr <- makeGRangesFromDataFrame(consensusTiles,keep.extra.columns=T )


            return(consensusTiles_gr)
}




### impose color scheme
scale_fill_scmacs <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('red', 
                            'green', 
                            'blue', 
                            'yellow',
                            'purple',
                            'black',
                            'gray'), 
                          c('Homer','Macs2',
                            'MOCHA',  'Homer_mocha',
                            'Homer_Macs2',
                            'Macs2_mocha',
                            'Common')), 
        ...
    )
}



###############################################################################
###############################################################################

### additional open regions    
get_pairwise_ratios <- function(sortedFreqs,
                                tileCountsPerMethod,
                                melt_TSS_res,
                                melt_CTCF_res,
                                Venn_Res,
                                Intensity_res){
    

    a = tileCountsPerMethod[CellPop %in% abundant_cells$CellSubsets , 
                        mean(value), 
                        by=variable]                               
    a$Ratio = 1/(a$V1/a$V1[a$variable=='MOCHA'])
    a$Model = factor(a$variable,
                     levels=c('Homer','Macs2','MOCHA'))

    ### additional TSS regions                                                       
    b = melt_TSS_res[CellPop %in% abundant_cells$CellSubsets, 
                     mean(value), 
                     by=variable]
    b$Ratio = 1/(b$V1/b$V1[b$variable=='MOCHA'])
    b$Model = factor(b$variable,
                     levels=c('Homer','Macs2','MOCHA'))                               

    ### additional CTCF regions                                                       
    c = melt_CTCF_res[CellPop %in% abundant_cells$CellSubsets, 
                     mean(value), 
                     by=variable]
    c$Ratio = 1/(c$V1/c$V1[c$variable=='MOCHA'])
    c$Model = factor(c$variable,
                     levels=c('Homer','Macs2','MOCHA'))         

    setkey(a, Model);setkey(b,Model);setkey(c,Model)   

    ### Combine first three results 
    results <- data.frame(cbind(a[,1:3],b[,2:3],c[,2:3]))
    colnames(results) <- c('Model','Peaks',
                           'PeakRatio','TSS','TSS_Ratio',
                           'CTCF','CTCF_Ratio')  




    ### Venn additional regions 
    Venn_Res = dcast(Venn_Res, CellPop ~ variable, direction='wide')      
    Venn_Res = Venn_Res[CellPop %in% abundant_cells$CellSubsets]
    Venn_Res$MOCHA_Macs2_Diff <- Venn_Res$MOCHA - Venn_Res$Macs2
    Venn_Res$MOCHA_Homer_Diff <- Venn_Res$MOCHA - Venn_Res$Homer                               



    ### Add difference Unique Regions
    results$DiffsPeakPopulation <- c(mean( Venn_Res$MOCHA_Homer_Diff), 
                                     mean( Venn_Res$MOCHA_Macs2_Diff),
                                     0)                                    


    ### Intensity Ratios 
    avgIntPerCell <- Intensity_res[CellPop %in% abundant_cells$CellSubsets,
                   mean(Value),
                  by=list(CellPop, Model)]

    ratio_of_means <- avgIntPerCell[ , 
                          list(Model=Model, 
                               Ratio = 1/(V1 / first(V1))), 
                          by=list(CellPop)]

    results$IntensityRatios <- 
    ratio_of_means[, mean(Ratio), by=Model]$V1[c(3,2,1)]

    return(results)
                               
}
###############################################################################
###############################################################################


###############################################################################
###############################################################################


### generate cumulative vector
### provides cumulative 
### aggregate results for 
### summarizing figure 2 results 

get_cumulative_vector <- function(value_vector,model){
    
        grid = sort(unique(value_vector ))
        vec <- sapply(grid, 
               function(x)
                   sum(value_vector <= x)
               )

        res=data.frame(
            Pct0=grid,
            values=vec,
            Model=model)

        return(res)

}

### generate cumulative result
### provides results to summarize
### cumulative accessibility 
### results for Figure 2 

generate_cumulative_result <- function(mocha, macs2, homer, Metric, cell){
    cumul_res= rbind(get_cumulative_vector(mocha,'MOCHA'),
                 get_cumulative_vector(macs2,'Macs2'),
                 get_cumulative_vector(homer,'Homer'))
    
    cumul_res$Metric =Metric
    cumul_res$Cell = cell
    
    return(cumul_res)   
}

### generate histogram result
### provides results to summarize
### accessibility 
### results for Figure 2 

generate_histogram_result <- function(mocha, macs2, homer, Metric, cell){
    hist_res<- data.table(
        values = c(mocha, macs2, homer),
        model = c(rep('MOCHA', length(mocha)),
                  rep('MACS2', length(macs2)),
                  rep('HOMER', length(homer))),
        Metric=Metric
        )
    hist_res$Cell = cell
    return(hist_res)   
}


###############################################################################
###############################################################################

### draw venn diagram
### for each method 
drawVennDiagram <- function(mocha, macs2, homer, cell){

        venn.diagram(
          x = list(mocha,
                   macs2,
                   homer),
          category.names = c("MOCHA" , "MACS2" , "HOMER"),
          filename =  paste(cell, 'venn.png',sep='_'),

          output=TRUE,
            disable.logging=TRUE,


                # Output features
                imagetype="png" ,
                height = 480 , 
                width = 480 , 
                resolution = 300,
                compression = "lzw",

                # Circles
                lwd = 2,
                lty = 'blank',
                fill = c('Blue','Green','Red'),

                # Numbers
                cex = .6,
                fontface = "bold",
                fontfamily = "sans",

                # Set names
                cat.cex = 0.6,
                cat.fontface = "bold",
                cat.default.pos = "outer",
                cat.pos = c(-27, 27, 135),
                cat.dist = c(0.055, 0.055, 0.085),
                cat.fontfamily = "sans",
                rotation = 1
        )  
}
    
###############################################################################
###############################################################################


summarize_celltypes <- function(i, numCores=10){
    
    #### Extract MOCHA tiles 
    subsetTileResults <- tileResults[[cells[i]]]

    sampleTileMatrix <- RaggedExperiment::compactAssay(subsetTileResults, i='peak')
    
    MOCHA_tileList <- mclapply(1:ncol(sampleTileMatrix),
                               function(x)
                                   sample_specific_rowRanges(sampleTileMatrix[,x]),
                               mc.cores=numCores
                               )

    tsam <- RaggedExperiment::compactAssay(tileResults[[cells[i]]], 
                                           i='TotalIntensity')
    colnames(tsam) <- gsub('_','-', colnames(tsam))
    ################################################################
    ################################################################


    ################################################################
    ################################################################
    ### 

    #### Extract Macs2 Tiles 
    setwd('Macs2')
    setwd(cell_dirs[i])

    ### Macs2 Peaks
    fnames <- dir()
    fnames <- fnames[grep('broadPeak', fnames)]

    ## read summits list
    macs2_tile_list <- mclapply(fnames,
                                function(x) fread(x, data.table=T),
                                mc.cores=numCores
                                )

    ### convert broad 
    ### peaks into accessible
    ### regions (tiles)

    macs2_tile_list2 <- mclapply(1:length(macs2_tile_list),

                                  function(x) 
                                     peak_to_tile_macs2_bulk(x, 
                                                         macs2_tile_list),
                                  mc.cores=numCores

                                  )
    
    tmpNames <- gsub('_peaks.broadPeak','',fnames)
    tmpNames <- gsub('_','',tmpNames)
    tmpCell = gsub(' ','.', cells[i])
    tmpNames <- gsub(paste('..',tmpCell,sep=''), '', tmpNames)
    tmpNames <- gsub('\\.','-', tmpNames)
    tmpNames <- gsub('#','-', tmpNames)

    
    macs2_tile_list3 <- mclapply(1:length(macs2_tile_list2),
                                 function(x){
                                     
                                     ## get peak calls
                                     macs2Peak =macs2_tile_list2[[x]]
                                     
                                     ## identify non-empty tiles
                                     ## for that sample
                                     nonEmptyTiles <- data.table(
                                         values =tsam[,tmpNames[x]],
                                         tileID =row.names(tsam))
                                     
                                     nonEmptyIdx = which(!is.na(nonEmptyTiles[,1]))
                                     nonEmptyTiles=nonEmptyTiles[nonEmptyIdx,]
                                     
                                     ## get intersection of macs2 tiles 
                                     ## and non-empty tiles
                                     new_macs2Peak <- macs2Peak[macs2Peak$tileID 
                                                                %in% nonEmptyTiles$tileID,]
                                     new_macs2Peak
                                     
                                     },
                                 mc.cores=10
                                 )
    macs2_tile_list2 <- macs2_tile_list3                             

    ################################################################
    ################################################################


    ################################################################
    ################################################################
    #### Extract Homer Tiles 
    setwd('../../Homer')
    setwd(cell_dirs[i])

    ### Homer Peaks
    fnames <- dir()
    peakNames <- paste(fnames,'/peaks.txt', sep='')
    peakNames <- peakNames[!peakNames %in% c("ArchRLogs/peaks.txt",'-style/peaks.txt')]

    ### extract tiles 
    homer_peak_list <- mclapply(peakNames,
                                        function(x)
                                            fread(x),
                                        mc.cores=numCores
                                        )

    ### Transform to tiles 
    homer_peak_list2 <- mclapply(1:length(homer_peak_list),
                                  function(x) peak_to_tile_homer_bulk(x, 
                                                                      homer_peak_list),
                                  mc.cores=15
                                  )
    
    ### translate names 
    tmpHomerName <- gsub('_|/peaks.txt','',peakNames)
    tmpCell = gsub(' ','.', cells[i])
    tmpHomerName <- gsub(paste('..',tmpCell,sep=''), '', tmpHomerName)
    tmpHomerName <- gsub('\\.','-',tmpHomerName)
    tmpHomerName <- gsub('#','-', tmpHomerName)
    
    homer_peak_list3 <- mclapply(1:length(homer_peak_list2),
                                 function(x){
                                     
                                     ## get peak calls
                                     homerPeak =homer_peak_list2[[x]]
                                     
                                     ## identify non-empty tiles
                                     ## for that sample
                                     nonEmptyTiles <- data.table(
                                         values =tsam[,tmpHomerName[x]],
                                         tileID =row.names(tsam))
                                     
                                     nonEmptyIdx = which(!is.na(nonEmptyTiles[,1]))
                                     nonEmptyTiles=nonEmptyTiles[nonEmptyIdx,]
                                     
                                     ## get intersection of macs2 tiles 
                                     ## and non-empty tiles
                                     new_homerPeak <- homerPeak[homerPeak$tileID 
                                                                %in% nonEmptyTiles$tileID,]
                                     new_homerPeak
                                     
                                     },
                                 mc.cores=10
                                 )
    homer_peak_list2 <- homer_peak_list3

    ################################################################
    ################################################################
    ### Panel B: Covid-19

    N =  ncol(sampleTileMatrix)

    tilesPerMethod<- data.table(
                Macs2 = sapply(macs2_tile_list2, function(x) nrow(x)
                            ),
                Homer = sapply(homer_peak_list2, function(x) nrow(x)
                            ),
                MOCHA = colSums(sampleTileMatrix, na.rm=T)
            )

    tilesPerMethod =  melt(tilesPerMethod)
    tilesPerMethod$CellPop = cells[i]
    
    homer_population <- getConsensusTiles(homer_peak_list2)
    macs2_population <- getConsensusTiles(macs2_tile_list2)
    mocha_population <- getConsensusTiles(MOCHA_tileList)

    ################################################################
    ################################################################
    
    mocha_tsam <- tsam[row.names(tsam) %in% mocha_population$tileID,]
    macs2_tsam <- tsam[row.names(tsam) %in% macs2_population$tileID,]
    homer_tsam <- tsam[row.names(tsam) %in% homer_population$tileID,]

    ####################################################
    ####################################################
    # calculate zeros 

    mocha_zeroes <- rowMeans(is.na(mocha_tsam))
    macs2_zeroes <- rowMeans(is.na(macs2_tsam))
    homer_zeroes <- rowMeans(is.na(homer_tsam))

    cumul_tiles= generate_cumulative_result(mocha_zeroes,
                               macs2_zeroes,
                               homer_zeroes,
                               'Tiles',
                               cells[i])
    
    hist_tiles= generate_histogram_result(mocha_zeroes,
                               macs2_zeroes,
                               homer_zeroes,
                               'Tiles',
                               cells[i])

    
    ##############################################################
    ### CTCF Per method 

    mocha_ctcf <- subsetByOverlaps(mocha_population, ctcf)
    macs2_ctcf <- subsetByOverlaps(macs2_population, ctcf)
    homer_ctcf <- subsetByOverlaps(homer_population, ctcf)

    mocha_ctcf_zeroes <- mocha_zeroes[names(mocha_zeroes) %in% mocha_ctcf$tileID]
    macs2_ctcf_zeroes <- macs2_zeroes[names(macs2_zeroes) %in% macs2_ctcf$tileID]
    homer_ctcf_zeroes <- homer_zeroes[names(homer_zeroes) %in% homer_ctcf$tileID]

    cumul_ctcf= generate_cumulative_result(mocha_ctcf_zeroes,
                               macs2_ctcf_zeroes,
                               homer_ctcf_zeroes,
                               'CTCF',
                               cells[i])
    
    hist_ctcf= generate_histogram_result(mocha_ctcf_zeroes,
                               macs2_ctcf_zeroes,
                               homer_ctcf_zeroes,
                               'CTCF',
                               cells[i])

    ##############################################################
    ### TSS Per method 

    mocha_tss <- subsetByOverlaps(mocha_population, new_tss)
    macs2_tss <- subsetByOverlaps(macs2_population, new_tss)
    homer_tss <- subsetByOverlaps(homer_population, new_tss)

    mocha_tss_zeroes <- mocha_zeroes[names(mocha_zeroes) %in% mocha_tss$tileID]
    macs2_tss_zeroes <- macs2_zeroes[names(macs2_zeroes) %in% macs2_tss$tileID]
    homer_tss_zeroes <- homer_zeroes[names(homer_zeroes) %in% homer_tss$tileID]

    cumul_tss= generate_cumulative_result(mocha_tss_zeroes,
                               macs2_tss_zeroes,
                               homer_tss_zeroes,
                               'TSS',
                               cells[i])
    
    hist_tss= generate_histogram_result(mocha_tss_zeroes,
                               macs2_tss_zeroes,
                               homer_tss_zeroes,
                               'TSS',
                               cells[i])
    ####################################################
    ####################################################
    ## calculate promoters
    ##
    require(tidyverse)
    require(magrittr)

    library(TxDb.Hsapiens.UCSC.hg38.refGene)
    library(org.Hs.eg.db)
    TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
    Org <- org.Hs.eg.db

    mocha_gr <- annotateTiles(mocha_population, TxDb, Org)
    macs2_gr <- annotateTiles(macs2_population, TxDb, Org)
    homer_gr <- annotateTiles(homer_population, TxDb, Org)
    
    macs2_homer_gr <- unique(c(macs2_gr, homer_gr))

    mocha_promoter_zeroes <- mocha_zeroes[names(mocha_zeroes) %in% 
                                          mocha_gr$tileID[mocha_gr$tileType=='Promoter']]
    macs2_promoter_zeroes <- macs2_zeroes[names(macs2_zeroes) %in% 
                                          macs2_gr$tileID[macs2_gr$tileType=='Promoter']]
    homer_promoter_zeroes <- homer_zeroes[names(homer_zeroes) %in% 
                                          homer_gr$tileID[homer_gr$tileType=='Promoter']]
  
    cumul_promoter= generate_cumulative_result(mocha_promoter_zeroes,
                               macs2_promoter_zeroes,
                               homer_promoter_zeroes,
                               'Promoter',
                               cells[i])
    
    hist_promoter= generate_histogram_result(mocha_promoter_zeroes,
                               macs2_promoter_zeroes,
                               homer_promoter_zeroes,
                               'Promoter',
                               cells[i])    
    
    mocha_unique_idx <- which(!mocha_gr$tileID %in% union(macs2_gr$tileID,
                                                         homer_gr$tileID))
    mocha_missed_idx <- which(!macs2_homer_gr$tileID %in%
                             mocha_gr$tileID)
    common_all_three <- which(mocha_gr$tileID %in% intersect(macs2_gr$tileID,
                                                         homer_gr$tileID))   
    
    mocha_unique_tileType =  as.data.table(mocha_gr[mocha_unique_idx])
    mocha_missed_tileType =  as.data.table(macs2_homer_gr[mocha_missed_idx])
    common_three_tileType =  as.data.table(mocha_gr[common_all_three])
    
    ### Get Counts 
    mocha_unique_tileType =  mocha_unique_tileType[, list(value=.N), by=tileType]
    mocha_missed_tileType =  mocha_missed_tileType[, list(value=.N), by=tileType]
    common_three_tileType =  common_three_tileType[, list(value=.N), by=tileType]
        
    ### capture intensity 
    tsam2 = as.data.table(tsam)
    tsam2$tiles <- row.names(tsam)
    
    ### unique lambda
    tsam = log2(tsam+1)
    unique_lambda= rowMeans(tsam[row.names(tsam) %in% mocha_gr$tileID[mocha_unique_idx],],
                            na.rm=T)
    
    ### missed lambda
    missedLambda= rowMeans(tsam[row.names(tsam) %in% macs2_homer_gr$tileID[mocha_missed_idx],],
                           na.rm=T)
    ### all three 
    all3Lambda= rowMeans(tsam[row.names(tsam) %in% mocha_gr$tileID[common_all_three],],
                           na.rm=T)    
    ####################################################
    ####################################################    
    require(VennDiagram)

    ### Draw venn diagram
    ### for all 3 methods 
    setwd(datasetDir)
    setwd('results')
    VennList <- list(mocha_population$tileID,
                    macs2_population$tileID,
                    homer_population$tileID)
    
    ### 
    results_hist <- rbind(hist_tiles,
                          hist_tss,
                          hist_ctcf,
                          hist_promoter)
    
    results_cumul <- rbind(cumul_tiles,
                           cumul_ctcf,
                           cumul_tss,
                           cumul_promoter)

    setwd(datasetDir)
    return(list(Hist= results_hist,
                Cumul=results_cumul,
                Tiles=tilesPerMethod,
                VennList=VennList,
                TileTypesUniqueMOCHA=mocha_unique_tileType,
                TileTypesAllThree = common_three_tileType,
                TileTypesMissedMOCHA=mocha_missed_tileType,
                MochaUnique_lambda=unique_lambda,
                MochaMissed_lambda= missedLambda,
                All3Lambda=all3Lambda
                
               ))
        
}





summarize_upset_plots <- function(i){
    
    venn = venns[[i]]
    counts <- venn[, .N, by=Method]
    
    mocha_tiles <- venn[Method=='MOCHA']$Tiles
    macs2_tiles <- venn[Method=='MACS2']$Tiles
    homer_tiles <- venn[Method=='HOMER']$Tiles
    
    mocha_unique <- sum(!mocha_tiles %in% union(macs2_tiles, homer_tiles))
    mocha_recall_macs2 <- sum(macs2_tiles%in% mocha_tiles)/length(macs2_tiles) 
    mocha_recall_homer <- sum(homer_tiles%in% mocha_tiles)/length(homer_tiles) 
    
    
    results <- data.frame(
        MOCHA = length(mocha_tiles),
        MACS2 = length(macs2_tiles),
        HOMER = length(homer_tiles),
        MOCHA_Unique = mocha_unique, 
        Union = length(Reduce(intersect, list(mocha_tiles, macs2_tiles, homer_tiles))),
        MACS2_Detected_C= sum(macs2_tiles%in% mocha_tiles),
        MACS2_Detected_P = mocha_recall_macs2,
        HOMER_Detected_C = mocha_recall_homer,
        HOMER_Detected_P = sum(homer_tiles%in% mocha_tiles),
        Cell = cells[i]
    )
        
        
    return(results)
        
}

    
blank_theme <- theme_minimal()+
           theme(
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.border = element_blank(),
           panel.grid=element_blank(),
           axis.ticks = element_blank(),
           plot.title=element_text(size=14, face="bold")
           )
    