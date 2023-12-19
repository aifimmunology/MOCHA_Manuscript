library(ArchR)
library(MOCHA)
library(stringr)
library(doParallel)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(tidyverse)

FullCovid <- loadArchRProject('FullCovidf')

FullCovid$predictedCellType <- FullCovid$predictedGroup_Co2
FullCovid$predictedCellType[grepl('ASDC|cDC|pDC', FullCovid$predictedCellType)] = 'DC'
FullCovid$predictedCellType[grepl('gdT|dnT', FullCovid$predictedCellType)] = 'OtherT'
FullCovid$predictedCellType[grepl('Treg', FullCovid$predictedCellType)] = 'Treg'
FullCovid$predictedCellType[grepl('CD4 CTL|CD4 TEM|CD4 TCM|CD4 Proliferating', FullCovid$predictedCellType)] = 'CD4 Effector'
FullCovid$predictedCellType[grepl('CD8 TEM|CD8 TCM|CD8 TEMRAS', FullCovid$predictedCellType)] = 'CD8 Effector'
saveArchRProject(FullCovid)


cellTypes <- grep('Plasma|ILC|HSPC', names(table(FullCovid$predictedCellType)), invert = TRUE, value=TRUE)

studySignal= median(FullCovid$nFrags)

## filter down to just COVID- for correlations?
tileResults <- MOCHA::callOpenTiles(
    FullCovid,
    cellPopLabel = "predictedCellType",
    cellPopulations= cellTypes,
    TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
    Org = org.Hs.eg.db,
    outDir = getOutputDirectory(FullCovid),
    studySignal= studySignal,
    verbose = TRUE,
    fast = TRUE, 
    numCores = 35
)


saveRDS(tileResults, 'tileResults_AllCellTypes.rds')

tileResults <- readRDS('tileResults_AllCellTypes.rds')


STM <- getSampleTileMatrix( 
    tileResults,
    groupColumn = 'Stage',
    threshold = 0.2,
    numCores = 40
)

saveRDS(STM, 'SampleTileMatrix_AllCellTypes.rds')

