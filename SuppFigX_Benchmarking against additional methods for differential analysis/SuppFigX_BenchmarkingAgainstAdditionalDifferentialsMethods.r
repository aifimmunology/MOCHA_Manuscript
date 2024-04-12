## Generate raw pseudobulk counts for other methods

library(ArchR)

tmp <- loadArchRProject('CD16s_EarlyInfectionComp2')

allPeaks <- read.csv('cd16_signac.csv')
allPeaks <- sub("-", ":", allPeaks[,1])

tmp <- addFeatureMatrix(tmp, features = MOCHA::StringsToGRanges(allPeaks), ceiling = 10^20, matrixName = 'RawMatrix',
                        binarize = FALSE)

summarizedCounts <- getGroupSE(
  ArchRProj = tmp,
  useMatrix = 'RawMatrix',
  groupBy = "Sample",
  divideN = FALSE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

rowRanges(summarizedCounts) <- MOCHA::StringsToGRanges(allPeaks)

saveRDS(summarizedCounts, 'SummarizedCounts_Raw_CD16Monocytes.rds')

### DESeq2
library(DESeq2)

## Followed tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## Subset down to early infection
summarizedCounts$COVID_status = summarizedCounts$DaysPSO < 0

dds <- DESeqDataSet(summarizedCounts, design = ~ COVID_status)
dds <- DESeq(dds)
res <- results(dds)

saveRDS(dds, 'DATS_DESeq2_Object.rds')

peakSet <- rowRanges(summarizedCounts)
mcols(peakSet) <- res

write.csv(peakSet, 'DATs_DESeq2_EarlyInfection.csv')

mat_deo <- counts(dds, normalized=TRUE)
rownames(mat_deo) <- rownames(summarizedCounts)
write.csv(mat_deo, 'NormMatrix_DESeq2_EarlyInfection.csv')


### Now with DiffChiPL
devtools::install_github("yancychy/DiffChIPL")
devtools::install_github("onofriAndreaPG/aomisc")

##Set up design
str1 = "EarlyInfection"
group= summarizedCounts$COVID_status
ctrName = "Uninfected"
treatName = "COVID19"
groupName = ifelse(group, ctrName, treatName)
design0 <- cbind(rep(1, length(summarizedCounts$COVID_status)), summarizedCounts$COVID_status)
colnames(design0) <- c(ctrName, treatName)

## Normalize
cpmD = cpmNorm(assays(summarizedCounts)[[1]])
rownames(cpmD) <- MOCHA::GRangesToString(rowRanges(summarizedCounts))

write.csv(cpmD, 'NorMatrix_DiffChipL_EarlyInfection.csv')

#Run Differential
resA = DiffChIPL(cpmD, design0, group0 = group)
saveRDS(resA, 'DATs_DiffChipL_FullDEObject.rds')
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE

write.csv(rtRlimm3, 'DATS_DiffChipL_Differentials.csv')

#####



### Now with DiffChiPL
#conda install bioconda::bioconductor-biocinstaller
#conda install bioconda::bioconductor-limma bioconda::bioconductor-sgseq bioconda::bioconductor-bamsignals bioconda::bioconductor-edger
#conda install r-devtools
#conda install bioconda::bioconductor-
devtools::install_github("yancychy/DiffChIPL")
devtools::install_github("onofriAndreaPG/aomisc")
library(DiffChIPL)

##Set up design
str1 = "EarlyInfection"
group= summarizedCounts$COVID_status
ctrName = "Uninfected"
treatName = "COVID19"
groupName = ifelse(group, ctrName, treatName)
design0 <- cbind(rep(1, length(summarizedCounts$COVID_status)), summarizedCounts$COVID_status)
colnames(design0) <- c(ctrName, treatName)

## Normalize
cpmD = cpmNorm(assays(summarizedCounts)[[1]])
rownames(cpmD) <- MOCHA::GRangesToString(rowRanges(summarizedCounts))

write.csv(cpmD, 'NorMatrix_DiffChipL_EarlyInfection.csv')

#Run Differential
resA = DiffChIPL(cpmD, design0, group0 = group)
saveRDS(resA, 'DATs_DiffChipL_FullDEObject.rds')
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE

write.csv(rtRlimm3, 'DATS_DiffChipL_Differentials.csv')

######### Now let's run permutations for both methods. 
permut_samples <- readRDS('permutation.RDS')

library(DESeq2)

dds <- readRDS('DATS_DESeq2_Object.rds')

diffList = list()

for(perSamp in permut_samples){
    tmpMeta <- perSamp$metadata
    dds$NewStatus <- factor(tmpMeta[colnames(dds),]$PermutedLabel)
    design(dds) = ~NewStatus
    newDDS <- nbinomWaldTest(dds)
    res <- results(newDDS)
    diffList = append(diffList, list(res))
}
saveRDS(diffList, 'DESeq2_Permutations.rds')


### This requires a different environment. Packages are not mutually compatible.
## conda activate "./env/MOCHA_Manuscript"
library(DiffChIPL)
cpmD <- read.csv('NorMatrix_DiffChipL_EarlyInfection.csv', row.names = 1)

ChIPList = lapply(permut_samples, function(perSamp){
     tmpMeta <- perSamp$metadata
    
    str1 = "EarlyInfection"
    group =  tmpMeta[gsub("\\.","-",colnames(cpmD)),]$PermutedLabel == 'Positive'
    ctrName = "Positive"
    treatName = "Negative"
    groupName = ifelse(group, ctrName, treatName)
    design0 <- cbind(rep(1, length(group)), group)
    colnames(design0) <- c(ctrName, treatName)

    DiffChIPL(cpmD, design0, group0 = group)
})
saveRDS(ChIPList, 'DiffChipL_Permutations.rds')  














