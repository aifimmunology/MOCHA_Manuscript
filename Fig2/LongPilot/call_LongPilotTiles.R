# ###################################################
# 0. Load libraries, ArchR project, and annotation
#    databases. Optionally filter the ArchR project
#    to a subset of samples.
# ###################################################

library(scMACS)
library(ArchR)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)


# You should substitute this with your own ArchR project.
# You must have completed cell labeling with your ArchR project.
ArchRProj <- ArchR::loadArchRProject("/home/jupyter/longPilot")
metadata = as.data.table(getCellColData(ArchRProj))
studySignal = median(metadata$nFrags)

# Define your annotation package for TxDb object(s)
# and genome-wide annotation
# Here our samples are human using hg38 as a reference.
# For more info: https://bioconductor.org/packages/3.15/data/annotation/
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

# ###################################################
# 1. Setting Parameters
#    These should be set according to your ArchR 
#    project and investgative question.
#
#    For more details on each of these parameters, 
#    view the help pages for each function using 
#    ?callOpenTiles and ?getSampleTileMatrix
# ###################################################

# Parameters for calling open tiles.
cellPopLabel <- "predictedGroup_Col2.5" 
cellPopulations <- c('CD16 Mono','B naive','CD8 TEM','CD14 Mono')
numCores <- 30

# Parameters for generating the sample-tile matrices
threshold <- 0.2

# ###################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
# ###################################################

tileResults <- callOpenTiles( 
    ArchRProj,
    cellPopLabel = cellPopLabel,
    cellPopulations = cellPopulations,
    TxDb = TxDb,
    Org = Org,
    numCores = 20,
    studySignal = studySignal,
    outDir=NULL
)

saveRDS(tileResults,
     file='/home/jupyter/MOCHA_Manuscript/Fig2/LongPilot/MOCHA.RDS')
