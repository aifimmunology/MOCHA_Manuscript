
library(MOCHA)
library(SummarizedExperiment)

setwd('scMACS_Analysis')

Spearman1 <- read.csv('CoAcc_Spearman.csv')
ZISpearman1 <- read.csv('CoAcc_ZI_Spearman.csv')

CorrDF <- read.csv('Correlation_Dataframe.csv')

STM <- readRDS('MOCHA_All_SampleTileMatrix.RDS')


gene1 <- rowRanges(STM)$Gene[match(CorrDF$Tile1, rownames(STM))] 
gene1[rowRanges(STM)$tileType != 'Promoter'] = NA

gene2 <- rowRanges(STM)$Gene[match(CorrDF$Tile2, rownames(STM))] 
gene2[rowRanges(STM)$tileType != 'Promoter'] = NA

CorrDF$Gene1 = gene1
CorrDF$Gene2 = gene2

subsetCorr <- CorrDF[(!is.na(CorrDF$Gene1) | !is.na(CorrDF$Gene2)) & (CorrDF$ZI_FDR < 0.05 | CorrDF$Norm_FDR < 0.05),]

subsetCorr[grepl('CD6', subsetCorr$Gene1) | grepl('CD6', subsetCorr$Gene2), ] 

totalObj <- combineSampleTileMatrix(STM)

totalMat <- assays(totalObj)[[1]]

## Identify strongest correlations that are different between ZI and Norm. 
ZISpecific <- arrange(CorrDF, desc(ZI_Correlation)) %>% filter(Norm_FDR >= 0.05 & ZI_FDR < 0.05) %>% dplyr::slice_head( n = 9)  #2358 total
NormSpecific <- arrange(CorrDF, desc(Norm_Correlation)) %>% filter(ZI_FDR >= 0.05 & Norm_FDR < 0.05) %>% dplyr::slice_head( n = 9)  #174 total


ZIPlots <- lapply(1:dim(ZISpecific)[1], function(x){
    
    plot_coaccessibility(ZISpecific, x, totalMat,fname = 'fileName', addLabel = TRUE)
    
    subMat <- data.frame(Tile1 = totalMat[ZISpecific$Tile1[x],], Tile2 = totalMat[ZISpecific$Tile2[x],])
    
    plot_coaccessibility


})

NormPlots <- lapply(1:dim(NormSpecific)[1], function(x){
    
    subMat <- data.frame(Tile1 = totalMat[NormSpecific$Tile1[x],], Tile2 = totalMat[NormSpecific$Tile2[x],])
    
    ggplot() + geom_point(data = subMat,aes(x = Tile1, y = Tile2)) +  
        geom_smooth(data = subMat, aes(x = Tile1, y = Tile2), method=lm, color = 'blue') +
        geom_smooth(data = subMat[rowSums(subMat == 0) == 0,], 
                    aes(x = Tile1, y = Tile2, color = 'red'), method=lm) +
    theme_bw() +theme(legend.position = "none") 


})

pdf('CorrelationAnalysis_StrongestCorrelations_Unique_ZISpearman_Spearman.pdf')

annotate_figure(ggpubr::ggarrange(plotlist = ZIPlots, ncol = 3, nrow = 3), 
                top = text_grob("Significant in Only ZI Spearman"))

annotate_figure(ggpubr::ggarrange(plotlist = NormPlots, ncol = 3, nrow = 3), 
                top = text_grob("Significant in Only Normal Spearman"))
dev.off()

write.csv(subsetCorr[grepl('CD6', subsetCorr$Gene1) | grepl('CD6', subsetCorr$Gene2), ], 'Correlations_CD6_TCell_Regulation.csv')

regAvg <- extractRegion(STM, region ='chr11:60969000- 61023499', cellPopulations = 'CD8 Naive', sampleSpecific = FALSE, numCores = 35)

saveRDS(regAvg, 'CD4Naive_region.rds')

plotRegion(regAvg, collapseGenes = 'longestTx', whichGene = 'CD6')

