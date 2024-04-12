

library(ArchR)
library(tidyverse)


setwd('scMACS_Analysis')

#### Let's run the ArchR permutations. 
cd16Mono <- loadArchRProject('CD16s_EarlyInfectionComp')
STM <- readRDS('CD16_SampleTileObj.rds')

addArchRThreads(15)
cd16Mono <- addPeakSet(cd16Mono,rowRanges(STM), force = TRUE)
saveArchRProject(cd16Mono)
cd16Mono <- addPeakMatrix(cd16Mono,  force = TRUE)
saveArchRProject(cd16Mono)

permut_samples <- readRDS('permutation.RDS')

diffObjList = list()
for(perSamp in permut_samples){
    
     newStatus = perSamp$metadata[cd16Mono$Sample,]$PermutedLabel
     cd16Mono <- addCellColData(cd16Mono, data = newStatus, name = 'NewStatus', 
                                cells = getCellNames(cd16Mono), force = TRUE)
     saveArchRProject(cd16Mono)
     tmpMark <- getMarkerFeatures(cd16Mono, 
                             groupBy = "NewStatus", 
                              useGroups = 'Positive',
                              bgdGroups = 'Negative',
                              useMatrix = 'PeakMatrix',
                              maxCells = 24744)
    diffObjList = append(diffObjList, list(tmpMark))
    gc()
}
saveRDS(diffObjList, 'ArchR_Differentials_PermutationTests.rds')


##########################################################################
###### Now for the same for Signac
###########################################################################
##
library(Seurat)
library(Signac)

seurObj <- readRDS('CD16_Mono_SignacObject.RDS')

permut_samples <- readRDS('permutation.RDS')

diffSigList = list()
for(perSamp in permut_samples){
    
     newStatus = perSamp$metadata[seurObj@meta.data$Sample,]$PermutedLabel
     seurObj@meta.data$NewStatus=newStatus
    Idents(seurObj) <- newStatus
     SeurDaps <- FindMarkers(
                  object = seurObj,
                  ident.1 = "Positive",
                  ident.2 = "Negative",
                  test.use = 'LR',
                  latent.vars = 'nFrags'
                )
    diffSigList  = append(diffSigList , list(SeurDaps))
    gc()
}
saveRDS(diffSigList, 'Seurat_Differentials_PermutationTests.rds')


##########################################################################
###### Now for the same for DESeq2 and DiffChipL
###########################################################################

######### Now let's run permutations for both methods. 


setwd('scMACS_Analysis')

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


##########################################################################
###### Now summarize across all methods
###########################################################################


## Load in all permutation results
ChIPList <- readRDS('DiffChipL_Permutations.rds')  
diffChIP <- lapply(ChIPList, function(XX){
        dplyr::filter(XX$resDE, adj.P.Val < 0.1)
    })

diffDE <- readRDS('DESeq2_Permutations.rds')  
diffDESeq <- lapply(diffDE, function(XX){
        XX <- as.data.frame(XX)
        XX$Feature = rownames(ChIPList[[1]])
        dplyr::filter(XX, padj < 0.1)
    })

diffArchR <- readRDS('ArchR_Differentials_PermutationTests.rds')
diffArchR <- lapply(diffArchR, function(XX){
         ArchR::getMarkers(XX,cutOff = 'FDR < 0.05')[[1]]
    })

diffSigList <- readRDS('Seurat_Differentials_PermutationTests.rds')
diffSig <- lapply(diffSigList, function(XX){
        dplyr::filter(XX, p_val_adj < 0.05)
    })

allDiffs <- list(diffChIP, diffDESeq, diffArchR, diffSig)
names(allDiffs) <- c("DiffChipL", "DESeq2", "ArchR", "Seurat")

allDiffDF <- do.call('rbind', lapply(seq_along(allDiffs), function(XX){
    
        tmpMat <- data.frame(NumberOfDATs = unlist(lapply(allDiffs[[XX]], function(ZZ) dim(ZZ)[1])))
        tmpMat$Method = names(allDiffs)[XX]
        tmpMat
    
    }))
                             
                         
mocha_list <- lapply(readRDS('permutation.RDS'), function(XX) 
        sum(XX[[1]]$FDR < 0.1))
allDiffDF <- rbind(allDiffDF, data.frame(
                    NumberOfDATs = unlist(mocha_list),
                    Method = rep('MOCHA', length(mocha_list)))
                   )
                     
write.csv(allDiffDF, 'Permutations_Across_DifferentialMethods.csv')
                    
allDiffDF <- read.csv('Permutations_Across_DifferentialMethods.csv')[,-1]
                     
###add in totals differentials
totalArchr <- sum(read.csv('cd16_ArchR.csv')$FDR < 0.05)
totalSignac <- sum(read.csv('cd16_signac.csv')$p_val_adj < 0.05)
totalDiffChip <- sum(read.csv('DATS_DiffChipL_Differentials.csv')$adj.P.Val < 0.05)
                     
totalDESEq <- sum(read.csv('DATs_DESeq2_EarlyInfection.csv')$padj <0.05, na.rm = TRUE) 
totalMOCHA <- sum(read.csv('cd16_mocha.csv')$FDR < 0.2)

df2 <- data.frame(Totals = c(totalArchr, totalSignac, totalDiffChip,totalDESEq, totalMOCHA), 
                  Method = c('ArchR', 'Signac', 'DiffChipL',
                             'DESeq2', 'MOCHA'))

allDiffDF$Method[allDiffDF$Method =='Seurat'] = 'Signac'                 
allDiffDF <- dplyr::full_join(allDiffDF, df2, by ='Method') %>%
                     dplyr::mutate(Rate = NumberOfDATs/Totals)
dplyr::group_by(allDiffDF, Method) %>% dplyr::summarize(AvgRate = median(Rate), MaxRate = max(Rate), MinRate = min(Rate))
                     
write.csv(allDiffDF, 'Permutations_Across_DifferentialMethods_withRate.csv')
  
allDiffDF <- read.csv('Permutations_Across_DifferentialMethods_withRate.csv')[,-1]
                     
pdf('Permutations_AcrossMethods.pdf')        
ggplot(allDiffDF, aes(x = Method, y = NumberOfDATs +1, color = Method)) + geom_jitter() + 
                             theme_bw() + xlab('Differential Method') + 
                             ylab('Number of Significant DATs') + 
                             theme(legend.position = 'none') +
                             scale_y_continuous(trans = 'log2')
ggplot(allDiffDF, aes(x = Method, y = Rate, color = Method)) + geom_jitter() + 
                             theme_bw() + xlab('Differential Method') + 
                             ylab('False Positive Rate')+
                             theme(legend.position = 'none') 
dev.off()
                     
cd16 <- read.csv('cd16 mono res.csv')[,-1]
cd4 <- read.csv('cd4 naive res.csv')[,-1]
bnaive <- read.csv('b naive res.csv')[,-1]
               
cd16$NB1_Statistic[is.na(cd16$NB1_Statistic)] = 1
cd16$Group = cd16$NB1_Pvalue < 0.05
                     
cd16Sum <- dplyr::mutate(cd16 ,Group =
        ifelse(NB1_Statistic>1, 'Overinflated', 'Underinflated')) %>%
        dplyr::group_by(Group) %>%
        dplyr::filter(NB1_Pvalue < 0.05) %>%
        dplyr::summarize(TileNum = dplyr::n())
                     
ggplot(dplyr::filter(cd16, NB1_Pvalue <0.05), aes(x =  NB1_Statistic>1)) + geom_bar() + theme_bw()