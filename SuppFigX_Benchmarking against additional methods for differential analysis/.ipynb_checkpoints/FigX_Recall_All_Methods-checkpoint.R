

library(ArchR)
library(tidyverse)


setwd('scMACS_Analysis')



#### Let's run the ArchR Recall. 
cd16Mono <- loadArchRProject('CD16s_EarlyInfectionComp')

down_samples <- readRDS('sampleList_downsampling.RDS')

archr_baseline <- getMarkerFeatures(cd16Mono, 
                         groupBy = "EarlyComparison", 
                          useGroups = 'CD16.Mono.Positive',
                          bgdGroups = 'CD16.Mono.Negative',
                          useMatrix = 'PeakMatrix',
                          maxCells = 24744)
saveRDS(archr_baseline, 'ArchR_baseline.rds')

diffObjList <- lapply(down_samples, function(perSamp){
       gc()
    lapply(perSamp, function(XX){

     cd16Mono2 <-subsetCells(cd16Mono, cellNames = 
                        getCellNames(cd16Mono)[cd16Mono$Sample %in% XX])
     getMarkerFeatures(cd16Mono2, 
                         groupBy = "EarlyComparison", 
                          useGroups = 'CD16.Mono.Positive',
                          bgdGroups = 'CD16.Mono.Negative',
                          useMatrix = 'PeakMatrix',
                          maxCells = 24744)
    
   
    })
})
saveRDS(diffObjList, 'ArchR_Differentials_Recall.rds')
  

##########################################################################
###### Now for the same for Signac
###########################################################################
##
library(Seurat)
library(Signac)

seurObj <- readRDS('CD16_Mono_SignacObject.RDS')

down_samples <- readRDS('sampleList_downsampling.RDS')
baseLine <- FindMarkers(object = seurObj,
                      ident.1 = "CD16.Mono.Positive",
                      ident.2 = "CD16.Mono.Negative",
                      test.use = 'LR',
                      latent.vars = 'nFrags'
                    )
saveRDS(baseLine, 'Seurat_baseline.rds')
diffSigList <- lapply(down_samples, function(perSamp){
       gc()
    pbapply::pblapply(cl = 10, perSamp, function(XX){
         subSeur <- subset(seurObj, subset = Sample %in% XX)
         FindMarkers(object = subSeur,
                      ident.1 = "CD16.Mono.Positive",
                      ident.2 = "CD16.Mono.Negative",
                      min.pct = 0.001,
                      logfc.threshold = 0.05,
                      test.use = 'LR',
                      latent.vars = 'nFrags'
                    )
        })

})
saveRDS(diffSigList, 'Seurat_Differentials_Recall.rds')


##########################################################################
###### Now for the same for DESeq2 and DiffChipL
###########################################################################

######### Now let's run Recall for both methods. 


setwd('scMACS_Analysis')


down_samples <- readRDS('sampleList_downsampling.RDS')

library(DESeq2)

dds <- readRDS('DATS_DESeq2_Object.rds')

diffList <- lapply(down_samples, function(perSamp){
        gc()
    pbapply::pblapply(cl = 10, perSamp, function(XX){
        newDDS <- nbinomWaldTest(dds[,XX])
        results(newDDS)
        })
})
saveRDS(diffList, 'DESeq2_Recall.rds')


### This requires a different environment. Packages are not mutually compatible.
## conda activate "./env/MOCHA_Manuscript"
library(DiffChIPL)
cpmD <- read.csv('NorMatrix_DiffChipL_EarlyInfection.csv', row.names = 1)

summarizedCounts <- readRDS('SummarizedCounts_Raw_CD16Monocytes.rds')

ChIPList = lapply(down_samples, function(perSamp){

    lapply(perSamp, function(XX){

        str1 = "EarlyInfection"
        ctrName = "Positive"
        treatName = "Negative"
        group= summarizedCounts[,
                        colnames(summarizedCounts) %in% XX]$COVID_status
        cpmD_tmp = cpmD[,colnames(summarizedCounts) %in% XX]
        
        groupName = ifelse(group, ctrName, treatName)
        design0 <- cbind(rep(1, length(group)), group)
        colnames(design0) <- c(ctrName, treatName)
        
        DiffChIPL(cpmD_tmp, design0, group0 = group)
    })
})
saveRDS(ChIPList, 'DiffChipL_Recall.rds')  



##########################################################################
###### Now summarize across all methods
###########################################################################
down_samples <- readRDS('sampleList_downsampling.RDS')

## Load in all permutation results
ChIP_base <- read.csv('DATS_DiffChipL_Differentials.csv', row.names = 1)
ChIP_base2 <- ChIP_base[grepl('chr4', rownames(ChIP_base)),]
ChIP_base2$adj.P.Val <- p.adjust(ChIP_base2$P.Value, 'fdr')
sigChIP <- which(ChIP_base2$adj.P.Val < 0.1)
ChIPList <- readRDS('DiffChipL_Recall.rds') 
ChIPRecall <- do.call('rbind', 
    lapply(seq_along(down_samples), function(YY){
        perSamp <- ChIPList[[YY]]
        downSampleNumber = unique(lengths(down_samples[[YY]]))
        print(YY)
        do.call('rbind', lapply(X = seq_along(perSamp), function(XX){
            print(XX)
                subDAT <- perSamp[[XX]]$resDE
                subDAT <- subDAT[grepl('chr4', rownames(subDAT)),]
                subDAT$adj.P.Val <- p.adjust(subDAT$P.Value, 'fdr')
                recall1 <- sum(which(subDAT$adj.P.Val < 0.1) %in% sigChIP)/length(sigChIP)
                newDAT1 <- sum(!which(subDAT$adj.P.Val < 0.1) %in% sigChIP)/sum(subDAT$adj.P.Val < 0.1)
                data.frame(SampleNumber = downSampleNumber, 
                           Iteration = XX, Method = 'DiffChIPL',
                           Recall = recall1, 
                           NewDifferentials = newDAT1)
            }))
    }))

library(DESeq2)
DE_base <- read.csv('DATs_DESeq2_EarlyInfection.csv', row.names = 1)
DE_base2 <- DE_base[grepl('chr4', DE_base$seqnames),]
DE_base2$padj <- p.adjust(DE_base2$pvalue, 'fdr')
sigDE <- which(DE_base2$padj < 0.1)
diffDE <- readRDS('DESeq2_Recall.rds')  
deRecall <- do.call('rbind', 
    lapply(seq_along(down_samples), function(YY){
        perSamp <- diffDE[[YY]]
        downSampleNumber = unique(lengths(down_samples[[YY]]))
        do.call('rbind', lapply(X = seq_along(perSamp), function(XX){
                subPerSamp = perSamp[[XX]][grepl('chr4', DE_base$seqnames),]
                subPerSamp$padj <- p.adjust(subPerSamp$pvalue, 'fdr')
                recall1 <- sum(which(subPerSamp$padj < 0.1) %in% sigDE)/length(sigDE)
                newDAT1 <- sum(! which(subPerSamp$padj < 0.1) %in% sigDE)/sum(subPerSamp$padj < 0.1, na.rm = TRUE)
                data.frame(SampleNumber = downSampleNumber, 
                            Method = 'DESeq2',
                           Iteration = XX, Recall = recall1, 
                           NewDifferentials = newDAT1)
            }))
    }))


library(ArchR)
diffArchR <- readRDS('ArchR_Differentials_Recall.rds')
archrBase1 <- readRDS('ArchR_baseline.rds')
archrBase1f <- archrBase1[rowData(archrBase1)$seqnames == 'chr4',]
assays(archrBase1f)[['FDR']][,1] = p.adjust(assays(archrBase1f)[['Pval']][,1], 'fdr')
archrBase <- ArchR::getMarkers(archrBase1f, cutOff = 'FDR < 0.05')[[1]]
sigArchR <- archrBase$idx

ArchRRecall <- do.call('rbind', 
    lapply(seq_along(down_samples), function(YY){
        perSamp <- diffArchR[[YY]]
        print(YY)
        downSampleNumber = unique(lengths(down_samples[[YY]]))
        gc()
        do.call('rbind', pbapply::pblapply(X = seq_along(perSamp), function(XX){
                tmpMarkers <- perSamp[[XX]]
                tmp1 <- rowData( tmpMarkers)
                tmpMarkersf <- tmpMarkers[tmp1$seqnames == 'chr4',]
                assays(tmpMarkersf)[['FDR']][,1] = p.adjust(assays(tmpMarkersf)[['Pval']][,1], 'fdr')
                tmpArchR <- ArchR::getMarkers(tmpMarkersf, cutOff = 'FDR < 0.05')[[1]]
                print(XX)
                recall1 <- sum(tmpArchR$idx %in% sigArchR)/length(sigArchR)
                newDAT1 <- sum(! tmpArchR$idx %in% sigArchR)/length(tmpArchR$idx)
                data.frame(SampleNumber = downSampleNumber, 
                            Method = 'ArchR',
                           Iteration = XX, Recall = recall1, 
                           NewDifferentials = newDAT1)
            }, cl = 3))
    }))

diffSigList <- readRDS('Seurat_Differentials_Recall.rds')
signacBase <- readRDS('Seurat_baseline.rds')
sigSeurat <- rownames(signacBase)[signacBase$p_val_adj < 0.05]
SeurRecall <- do.call('rbind', 
    lapply(seq_along(down_samples), function(YY){
        perSamp <- diffSigList[[YY]]
        downSampleNumber = unique(lengths(down_samples[[YY]]))
        do.call('rbind', lapply(X = seq_along(perSamp), function(XX){
                recall1 <- sum(rownames(perSamp[[XX]])[perSamp[[XX]]$p_val_adj < 0.05] %in% 
                               sigSeurat)/length(sigSeurat)
                newDAT1 <- sum(!rownames(perSamp[[XX]])[perSamp[[XX]]$p_val_adj < 0.05] %in% 
                               sigSeurat)/sum(perSamp[[XX]]$p_val_adj < 0.05)
                data.frame(SampleNumber = downSampleNumber, 
                            Method = 'Seurat',
                           Iteration = XX, Recall = recall1, 
                           NewDifferentials = newDAT1)
            }))
    }))

mocha_recall <- read.csv('downsamplingRes.csv')[,-1][-211,]
mocha_recall$Iteration = ArchRRecall$Iteration
mocha_recall$Method = 'MOCHA'
mocha_recall <- dplyr::rename(mocha_recall[,c('N', 'Method', 'Iteration', 'Recall', 'FalsePostive')],
                              SampleNumber = N, NewDifferentials = FalsePostive)

allRecall <- do.call('rbind', list(ArchRRecall, ChIPRecall, deRecall, SeurRecall, mocha_recall))
write.csv(allRecall, 'Recall_Summary.csv')
allRecall <- read.csv('Recall_Summary.csv')[,-1]
allRecall$Method[allRecall$Method == 'Seurat'] = 'Signac'

allRecall <- rbind(data.frame(
        SampleNumber = rep(39, 5),
        Method = c('Signac', 'MOCHA', 'DiffChIPL', 'ArchR', "DESeq2"),
        Iteration = rep(1,5),
        Recall = rep(1,5),
        NewDifferentials = rep(0,5)
        ), 
      allRecall)

pdf('Method_Recall.pdf')
ggplot(allRecall, aes(x = SampleNumber, y = Recall, color = Method, group = Method)) + 
    geom_smooth() + geom_jitter()+theme_bw() + xlab('Number of Samples') + ylab('Recall')

ggplot(allRecall, aes(x = SampleNumber, y = Recall, color = Method, group = Method)) + 
    geom_smooth() + geom_jitter()+theme_bw() + xlab('Number of Samples') + ylab('Recall') +
facet_wrap(~Method) + theme(legend.position = 'none')

dev.off()