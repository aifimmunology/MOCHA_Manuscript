# ###########################################################
# ###########################################################
#
# Author: Samir Rachid Zaim
# Date: 09/26/2022
#
# ###########################################################
# ###########################################################

## Load Libraries
require(MOCHA)
require(ArchR)

## Load the ArchR Project
## and extract the metadata object
ArchRProj <- loadArchRProject("/home/jupyter/FullCovid")
metadata = as.data.table(getCellColData(ArchRProj))
studySignal = median(metadata$nFrags)

# Define your annotation package for TxDb object(s)
# and genome-wide annotation
# Here our samples are human using hg38 as a reference.
# For more info: https://bioconductor.org/packages/3.15/data/annotation/

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

###########################################################
###########################################################

## Get metadata information
## at the sample level
lookup_table <- unique(metadata[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table$Sample[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]

## subset ArchR Project
idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
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
    numCores = 50,
    outDir = NULL,
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
    threshold = 0.2,
    verbose=T
)

###########################################################
# 4. Get differential accessibility for specific 
#    cell populations. Here we are comparing MAIT  
#    cells between samples where our groupColumn 
#    "COVID_status" is Positive (our foreground) 
#    to Negative samples (our background).
###########################################################

differentials <- getDifferentialAccessibleTiles(
    SampleTileObj = SampleTileMatrices,
    cellPopulation = "CD16 Mono",
    groupColumn = "COVID_status",
    foreground =  "Positive",
    background =  "Negative",
    fdrToDisplay = 1,
    outputGRanges = TRUE,
    numCores = 40
)

cd16_differentials <- differentials[FDR <= 0.2]
write.csv(cd16_differentials,
          file='cd16_mocha.csv')



unFilteredSummaries <- getDifferentialAccessibleTiles(
    SampleTileObj = SampleTileMatrices,
    cellPopulation = "CD16 Mono",
    groupColumn = "COVID_status",
    foreground =  "Positive",
    background =  "Negative",
    fdrToDisplay = 0.2,
    outputGRanges = FALSE,
    numCores = 40,
    noiseThreshold = 0,
    minZeroDiff = 0
)

write.csv(unFilteredSummaries,
          file='cd16_allTiles.csv')

###########################################################
# 5. Save output from tile Matrix
#    to run co-accessibility comparisons
###########################################################
sampleTileMatrix <- as.data.table(scMACS::getCellPopMatrix(SampleTileMatrices, 'CD16 Mono'))

sampleTileMatrix$tileID <- rownames(SampleTileMatrices)

write.csv(sampleTileMatrix,
          'tileMatrix.csv',
         row.names=F)