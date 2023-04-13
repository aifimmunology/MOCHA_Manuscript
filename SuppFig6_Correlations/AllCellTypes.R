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

STM2 <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList = 'Positive')

accMat <- assays(STM2)[['CD16 Mono']][mcols(rowRanges(STM))[,5],]

metaData <- colData(STM)

startTime <- Sys.time()
lmemList <- randomEffectModeling(metaData, accMat[c(1:100),], c('Age', 'Sex', 'days_since_symptoms'),
                                 numCores = 15)
stopTime <- Sys.time()
saveRDS(lmemList, 'CD16_Mono_ModelingList.rds')
varDecomp <- getVarDecomp(lmemList, numCores= 15)

write.csv(varDecomp, 'CD16_Mono_VarDecomp.csv')

 df <-  data.frame(exp = as.numeric(accMat[1,]), 
                metaData, stringsAsFactors = FALSE)
            list(df, formula1)

lme4:lmer(formula = formula1, data = df)
a
runModeling(list(df, formula1))
            
    varForm <- paste0(unlist(lapply(c('Age', 'Sex', 'days_since_symptoms'), function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
    formula1 <- as.formula(paste('exp ~ ',varForm, sep = ''))

    formula1 <- as.formula(paste('exp ~ ',varForm, sep = ''))

    MetaDF <- dplyr::filter(as.data.frame(MetaDF), Sample %in% colnames(CountDF))
    CountDF <- CountDF[,match(colnames(CountDF), MetaDF$Sample)]


startTime <- Sys.time()
lmemList <- randomEffectModeling(metaData, accMat[c(1:100),], c('Age', 'Sex', 'days_since_symptoms'),
                                 numCores = 15)
stopTime <- Sys.time()

varComp1 <- getVarDecomp(lmemList)
varComp2 <- getVarDecomp(lmemList)

tmp <- runDecomposition(metaData, accMat[c(1:100),], c('Age', 'Sex', 'days_since_symptoms'),
                                 numCores = 15)
