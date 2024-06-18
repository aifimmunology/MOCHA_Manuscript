# -*- coding: utf-8 -*-
# ###########################################################
# ###########################################################
#
# Author: Samir Rachid Zaim
# Date: 11/06/2021
#
# ###########################################################
# ###########################################################

## Load Libraries
require(MOCHA)
require(data.table)
require(ArchR)
require(ggpubr)

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

## Load the ArchR Project
## and extract the metadata object
ArchRProj <- loadArchRProject("/home/jupyter/FullCovid")
metadf <- as.data.table(ArchRProj@cellColData)
studySignal = median(metadf$nFrags)


## Subset to visit one 
lookup_table <- unique(metadf[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

lookup_table <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms),]



# ##########################################################
# ##########################################################

## subset ArchR Project
idxSample <- BiocGenerics::which(ArchRProj$Sample %in% lookup_table$Sample)
cellsSample <- ArchRProj$cellNames[idxSample]
ArchRProj <- ArchRProj[cellsSample, ]


###########################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
###########################################################


tileResults <- callOpenTiles( 
    ArchRProj,
    cellPopLabel = "CellSubsets" ,
    cellPopulations = "CD16 Mono",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
    Org = "org.Hs.eg.db",
    numCores = 20,
    outDir =NULL,
    studySignal = studySignal
)


###########################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
###########################################################

SampleTileMatrices <- getSampleTileMatrix( 
    tileResults,
    cellPopulations = "CD16 Mono",
    groupColumn = "COVID_status",
    threshold = 0.2
)

sampleTileMatrix <- MOCHA::getCellPopMatrix(SampleTileMatrices, 'CD16 Mono')

## subset ArchR Project
group = ifelse(colnames(sampleTileMatrix)[1:39] %in% 
               lookup_table$Sample[lookup_table$COVID_status=='Positive'],
               1,
               0)
###########################################################
###########################################################


###########################################################
###########################################################
peaks_to_test <- c(3, 229,
                   4210, 20966,
                   55251, 110009,
                   162531,190779,
                   204103)

estimate_runtime <- function(tilesTested){
       
    idx = sample(nrow(sampleTileMatrix),tilesTested)
    tmp_sampleTileMatrix = SampleTileMatrices[idx]
      
    time <- system.time(
        
        
        differentials <- getDifferentialAccessibleTiles(
            SampleTileObj = tmp_sampleTileMatrix,
            cellPopulation = "CD16 Mono",
            groupColumn = "COVID_status",
            foreground =  "Positive",
            background =  "Negative",
            fdrToDisplay = 0.2,
            outputGRanges = FALSE,
            numCores = 60
        )
        
    )
    
    time         
    
}

timeRuns <- lapply(peaks_to_test,
       function(x){
           estimate_runtime(x)
     
           }
       )

timeRuns = data.table(do.call(rbind, timeRuns))
timeRuns$PeaksTested = peaks_to_test
timeRuns$RunTime = timeRuns$elapsed
timeRuns$Model = 'MOCHA'
timeRuns =  timeRuns[, c('PeaksTested','RunTime','Model')]


write.csv(timeRuns,       
          file='/home/jupyter/MOCHA_Manuscript/Fig3/panelC_Runtime/MOCHA_runtime.csv',
          row.names=F
          )
          
    

