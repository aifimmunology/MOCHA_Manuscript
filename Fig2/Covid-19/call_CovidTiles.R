####################################################
# 0. Load libraries, ArchR project, and annotation
#    databases. Optionally filter the ArchR project
#    to a subset of samples.
####################################################

library(MOCHA)
library(ArchR)

# You should substitute this with your own ArchR project.
# You must have completed cell labeling with your ArchR project.
ArchRProj <- ArchR::loadArchRProject("/home/jupyter/FullCovid")
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


####################################################
# 1. Setting Parameters
#    These should be set according to your ArchR 
#    project and investgative question.
#
#    For more details on each of these parameters, 
#    view the help pages for each function using 
#    ?callOpenTiles and ?getSampleTileMatrix
####################################################

# Parameters for calling open tiles.
cellPopLabel <- "CellSubsets" 
cellPopulations <- c('CD16 Mono','B naive','CD4 CTL TEM')
numCores <- 30

# Parameters for generating the sample-tile matrices
threshold <- 0.2
groupColumn <- "COVID_status"

####################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
####################################################

tileResults <- callOpenTiles( 
    ArchRProj,
    cellPopLabel = "CellSubsets" ,
    cellPopulations = cellPopulations,
    TxDb = 'TxDb.Hsapiens.UCSC.hg38.refGene',
    Org = 'org.Hs.eg.db',
    numCores = 20,
    studySignal = studySignal,
    outDir=NULL
)

saveRDS(tileResults,
     file='/home/jupyter/MOCHA_Manuscript/Fig2/Covid-19/MOCHA.RDS')
