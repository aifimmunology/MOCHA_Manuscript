library(ArchR)
library(MOCHA)
library(stringr)
library(doParallel)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(tidyverse)

setwd('scMACS_Analysis')


MonoDCE <- loadArchRProject('MonoDC_Edits')
studySignal = 3628

####################################################
# 2. Call open tiles (main peak calling step)
#    Done once for all specified cell populations
####################################################


tileResults <- MOCHA::callOpenTiles(
    MonoDCE,
    cellPopLabel = "predictedGroup_Co2",
    cellPopulations= 'CD16 Mono',
    TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
    Org = org.Hs.eg.db,
    studySignal = studySignal,
    numCores = 25
)

saveRDS(tileResults, 'scMACS_tileResults_CD16.RDS')

SampleTileObj <-  MOCHA::getSampleTileMatrix( 
    tileResults,
    groupColumn= 'InfectionStages',
    threshold = 0.2,
    numCores = 35,
)
SampleTileObj <- annotateTiles(SampleTileObj)

saveRDS(SampleTileObj, 'scMACS_SampleTileMatrix.RDS')


###############################
### Add just time-related peaks
    
    
varPeaks <- read.csv('variance_decomposition_peaks.csv') 

subVarPeaks <- varPeaks[,c('PTID','Time','Age','Sex','Residual')]
maxComp <- colnames(subVarPeaks)[apply(subVarPeaks, 1, which.max)]
varPeaks$MaxVar <- maxComp
                                
Time_Regs <- filter(varPeaks, MaxVar == 'Time')
TimeRegs <- Time_Regs$Gene %>% sub("\\.",":",.) %>%
            sub("\\.","-", .) %>% StringsToGRanges()    

CD16s <- subsetArchRProject(MonoDCE, 
                            cells = getCellNames(MonoDCE)[
                                MonoDCE$predictedGroup_Co2 == 'CD16 Mono'],
                            outputDirectory = 'CD16sLong', dropCells = TRUE)
saveArchRProject(CD16s)
CD16s <- loadArchRProject('CD16sLong')

CD16s <- addPeakSet(CD16s, TimeRegs, force = TRUE)
CD16s <- addPeakMatrix(CD16s, force = TRUE)
addArchRThreads(50)
CD16s <- addIterativeLSI(CD16s,  useMatrix = "PeakMatrix", name = "TimeLSI",   iterations = 5, force = TRUE)
CD16s <- addUMAP(CD16s,  reducedDims ="TimeLSI",
  name = "Time_UMAP", force = TRUE)
CD16sUMAP <- getEmbedding(CD16s, 'Time_UMAP')

## Run LSI on tile matrix without specifically choosing time-related peaks
CD16s <- addIterativeLSI(CD16s, name = "normLSI", force = TRUE)
CD16s <- addUMAP(CD16s,  reducedDims ="normLSI",
  name = "normUMAP", force = TRUE)

kClustList <- mclapply(c(3:10), function(x){
    
                kmeans(CD16sUMAP, centers =x, iter.max = 1000, nstart = 100)
    
    }, mc.cores = 25)
kClustdf <- lapply(kClustList, function(x) x$cluster) %>% do.call('cbind',.)
colnames(kClustdf) <- paste('KCluster_',3:10,sep="")
                   
all(rownames(kClustdf) == rownames(CD16sUMAP))
colnames(CD16sUMAP) <- c('UMAP1', 'UMAP2')
UMAPdf <- cbind(CD16sUMAP, kClustdf)
                   
library(ggrastr)
pdf('CD16s_KClusterLongitudinal.pdf')
for(i in 3:10){
                   
               
    p1 <- ggplot(UMAPdf, aes_string(x = 'UMAP1', y ='UMAP2',
                            color = paste('KCluster_',i,sep=""))) + geom_point_rast() 
    print(p1)
                   
}
dev.off()
    
CD16s <- addClusters(CD16s,  reducedDims ="TimeLSI",
  name = "Clusters3", force = TRUE)
saveArchRProject(CD16s)

daysPSO <- CD16s$DaysPSO
daysPSO[daysPSO == -50] = -1


CD16s <- addCellColData(CD16s, data = daysPSO, name = "DaysPSO", cells = getCellNames(CD16s), force= TRUE)
CD16s <- addCellColData(CD16s, data = kClustdf[,1], name = "kClusters", cells = rownames(kClustdf), force= TRUE)
saveArchRProject(CD16s)

pdf('All_CD16s_LongitudinalUMAP.pdf')                   
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "DaysPSO")
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "InfectionStages")
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "Clusters3")    
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "kClusters")    
dev.off()
                   
MonocleTraj <- getMonocleTrajectories(
                  CD16s,
                  name = "Trajectory",
                useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  embedding = 'Time_UMAP')
                    
CD16s <- addMonocleTrajectory(CD16s, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                 monocleCDS = MonocleTraj,
                                  groupBy = 'InfectionStages')
                   
                  
saveArchRProject(CD16s)
                   
trajCells <- getCellColData(CD16s) %>% as.data.frame() %>% mutate(Cells = rownames(.)) %>%
                   dplyr::select(Cells, Trajectory, InfectionStages) %>% arrange(Trajectory)
ggplot(trajCells, aes(x = Trajectory, y =1, fill = InfectionStages)) + geom_tile()
                   
trajGS <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory', useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
    
trajDev <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
    
    

pdf('CD16s_Fig5_Trajectories.pdf')

plotTrajectory(CD16s, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS,  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev,  pal = paletteContinuous(set = "horizonExtra"))

dev.off()
            
                  
### Find Monocle trajectories and run Over-representation analysis
GS_Traj <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory', useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
ChromVar_Traj <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory', useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE) 

allGenes <- getGenes(CD16s)
allMotifs <- getPositions(CD16s, 'CD16LongCISBPMotif')
enrichDataBaseList <- WebGestaltR::listGeneSet()
                   
GeneEnrichment <- lapply(enrichDataBaseList$name[c(2,7,8,9,10)], function(x){
    
        simplifiedORA(x, unique(rowData(GS_Traj)$name), unique(allGenes$symbol))
    
    })
                   
specMotifList <- unique(rowData(ChromVar_Traj)$name)
                   
TF_Enrichment <- lapply(enrichDataBaseList$name[c(2,7,8,9,10)], function(x){
    
        simplifiedORA(x, gsub("_.*", "", specMotifList), 
                      gsub("_.*", "", names(allMotifs)))
    
    })                   

pdf('CD16_Fig5_Panther_Database.pdf')


ggplot(GeneEnrichment[[3]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Panther Pathway Enrichment on Genescore Trajectory')

ggplot(TF_Enrichment[[3]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Panther Pathway Enrichment on ChromVar Trajectory')

ggplot(GeneEnrichment[[4]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Reactome Pathway Enrichment on GeneScore Trajectory')

                   
dev.off()


############### Now do the same as above, but with the standard UMAP
library(monocle3)
MonocleTraj2 <- getMonocleTrajectories(
                  CD16s,
                  name = "Trajectory2",
                useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  embedding = "normUMAP")
                    
CD16s <- addMonocleTrajectory(CD16s, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                 name = 'Trajectory2',
                                 monocleCDS = MonocleTraj2,
                                  groupBy = 'InfectionStages', force = TRUE)
                   
                  
saveArchRProject(CD16s)
         
trajectoryAnnotation = data.frame(DaysPSO = CD16s$days_since_symptoms, Trajectory = CD16s$Trajectory2) %>%
                   dplyr::filter(!is.na(DaysPSO)) %>% dplyr::mutate(y = 1)
                   
pdf('CD16s_Fig5_DaysPSO_Over_Trajectory.pdf')
ggplot(trajectoryAnnotation, aes(x = Trajectory, y = y, fill = DaysPSO)) + geom_tile() + theme_minimal()
dev.off()
                   
       
### Pull out trajectories for GeneScores and Motif Deviations
                   
trajGS2 <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory2', 
                         useMatrix = "GeneScoreMatrix", log2Norm = TRUE)

trajDev2 <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory2',
                         useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
    

pdf('CD16s_Fig5_Default_Trajectories.pdf')

plotTrajectory(CD16s, name = "Trajectory2", trajectory = "Trajectory2", 
                embedding = "normUMAP")[[1]]
plotTrajectoryHeatmap(trajGS2,  pal = paletteContinuous(set = "horizonExtra"))
plotTrajectoryHeatmap(trajDev2,  pal = paletteContinuous(set = "horizonExtra"))

dev.off()

### For default UMAP trajectory,  Over-representation analysis 
allGenes <- getGenes(CD16s)
allMotifs <- getPositions(CD16s, 'CD16LongCISBPMotif')
enrichDataBaseList <- WebGestaltR::listGeneSet()
                   
GeneEnrichment <- lapply(enrichDataBaseList$name[c(2,7,8,9,10)], function(x){
    
        simplifiedORA(x, unique(rowData(trajGS2)$name), unique(allGenes$symbol))
    
    })
                   
specMotifList <- unique(rowData(trajDev2)$name)
                   
TF_Enrichment <- lapply(enrichDataBaseList$name[c(2,7,8,9,10)], function(x){
    
        simplifiedORA(x, gsub("_.*", "", specMotifList), 
                      gsub("_.*", "", names(allMotifs)))
    
    })                   

simplifiedORA(enrichDataBaseList$name[8], gsub("_.*", "", specMotifList), 
                      gsub("_.*", "", names(allMotifs)))
                   
pdf('CD16_Fig5_Panther_Database_DefaultUMAP.pdf')


ggplot(GeneEnrichment[[3]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Panther Pathway Enrichment on Genescore Trajectory')

ggplot(TF_Enrichment[[3]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Panther Pathway Enrichment on ChromVar Trajectory')

ggplot(GeneEnrichment[[4]], aes(x = description,
                           y = FDR, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('FDR') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Reactome Pathway Enrichment on GeneScore Trajectory')

                   
dev.off()




 


##### Let's try a different type of density plot

   
allDen <- plotDensity(ArchRProj= CD16s, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
           
NormalDen <- plotDensity(ArchRProj= CD16s, embeddingName= "normUMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)


library('ggpubr')
    
plotAll <- ggpubr::ggarrange(plotlist = allDen, nrow = 2, ncol = 2)
plotNormAll <- ggpubr::ggarrange(plotlist = NormalDen, nrow = 2, ncol = 2)




    
pdf('CD16_LongitudinalUMAP_Density.pdf')
annotate_figure(plotAll, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))
           
annotate_figure(plotNormAll, top = text_grob("CD16s in All Donors (Unbiased UMAP)", 
              face = "bold", size = 14))    

dev.off()
    

################################################################################


#### Now for pseudobulk analysis

#################################################################################

STM <- readRDS('scMACS_SampleTileMatrix.RDS')

sampleMat <- assays(STM)[['CD16 Mono']]
compMat <- cbind(as.data.frame(StringsToGRanges(results$SampleTileMatrix$tileID)), sampleMat)
name(compMat) <- 'counts'

## getMotifs

data(human_pwms_v2)
STM <- addMotifSet(STM, pwms = human_pwms_v2, motifSetName = 'Cisbp', w = 7)
saveRDS(STM,'scMACS_SampleTileMatrix.RDS')

MonoDCE <- loadArchRProject('MonoDC_Edits')
metadf <- getCellColData(MonoDCE) %>% as.data.frame()
subMat <- dplyr::select(metadf, c(Sample, PTID, DaysPSO, InfectionStages, PASC, Age, AgeGroups)) %>% distinct()
#subMat$Sample <- gsub("-",".", subMat$Sample)
subMat <- subMat %>% group_by(PTID) %>%
            mutate(PASC = case_when(any(grepl('Non',PASC)) ~ 'Noninflamed PASC', 
                                    any(grepl('Infl',PASC)) ~ 'Inflamed PASC',
                                      TRUE ~ 'Recovered')) %>%
            ungroup() %>%
            mutate(PASC_Cat = ifelse(grepl("PASC",PASC), 
                                 'PASC', 'Recovered')) %>% 
            ungroup() 

## Pull in Variance decomposition results
varPeaks <- read.csv('variance_decomposition_peaks.csv') 
subVarPeaks <- varPeaks[,c('PTID','Time','Age','Sex','Residual')]
maxComp <- colnames(subVarPeaks)[apply(subVarPeaks, 1, which.max)]
varPeaks$MaxVar <- maxComp
Time_Regs <- filter(varPeaks, MaxVar == 'Time')
TimeRegs <- Time_Regs$Gene %>% sub("\\.",":",.) %>%
            sub("\\.","-", .) %>% MOCHA::StringsToGRanges()

samplePeakMat <- sampleMat
colnames(samplePeakMat) <- gsub("__.*","", colnames(samplePeakMat)
AllPeaks <- MOCHA::StringsToGRanges(samplePeakMat$tileID)
                        
samplePeakMat[is.na(samplePeakMat)] <- 0
all_pca <- prcomp(t(samplePeakMat[rownames(samplePeakMat) %in% GRangesToString(TimeRegs),])) 
pca_plot <- all_pca$x  %>% as.data.frame() %>% dplyr::mutate(Sample = rownames(.)) %>%
                inner_join(subMat, by = 'Sample') %>% arrange(DaysPSO)
                                
pdf('Psuedobulk_PCA_Longitudinal_CD16s.pdf')
         
ggplot(pca_plot) + 
     geom_point(aes(x = PC1, y = PC2, shape = InfectionStages),
               size = 2.5) +
     theme_bw()
                                
ggplot(pca_plot) + 
    geom_point(aes(x = PC2, y = PC3, shape = InfectionStages),
               size = 2.5) +
     theme_bw()
                                
ggplot(pca_plot) + 
     geom_point(aes(x = PC3, y = PC4, shape = InfectionStages),
               size = 2.5) +
     theme_bw()  
            
ggplot(pca_plot) + 
     geom_point(aes(x = PC3, y = PC4, shape = InfectionStages),
               size = 2.5) +
     theme_bw()  

dev.off()
                                
TimeMatv2 <- samplePeakMat[samplePeakMat$tileID %in% GRangesToString(TimeRegs), -'tileID'] 
colnames(TimeMatv2) <- gsub("__.*","",colnames(TimeMatv2))
TimeMatv2 <- TimeMatv2[,colnames(TimeMatv2) %in% 
                     subMat$Sample[subMat$InfectionStages != "Uninfected"], with = FALSE]
set.seed(1)
TimeUMAPv2 <- as.data.frame(uwot::umap(t(TimeMatv2)))
colnames(TimeUMAPv2) = c('UMAP1', 'UMAP2')
TimeUMAPv2$Sample = rownames(TimeUMAPv2)
TimeUMAP2v2 <- inner_join(TimeUMAPv2, subMat, by = 'Sample') %>% arrange(DaysPSO)
TimeUMAP2v2$KMeans <- kmeans(TimeUMAP2v2[,c('UMAP1', 'UMAP2')], centers =3, nstart = 100)$cluster
TimeUMAP2v2 <- TimeUMAP2v2 %>% group_by(PTID) %>%
                    dplyr::mutate(Last = ifelse(DaysPSO == max(DaysPSO),
                                                'Last', 'NotLast')) %>%
                    dplyr::mutate(TrajectoryGroup = 
                                  ifelse(any(KMeans == 1 & Last == 'Last'),
                                         'Group 1', 'Group 2'))
                                
write.csv(TimeUMAP2v2, 'Pseudobulk_UMAP_Longitudinal_CD16s.csv')
                              
TimeUMAP2v2 <- read.csv('Pseudobulk_UMAP_Longitudinal_CD16s.csv')  %>% 
                                group_by(PTID) %>% arrange(PTID, DaysPSO) %>%
                     dplyr::mutate(Last = ifelse(DaysPSO == min(DaysPSO), 'First', Last))

allPeakUMAP <- uwot::umap(t(sampleMat)) %>% as.data.frame() %>% mutate(Sample = rownames(.)) %>%
                    dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
                    inner_join(subMat, by = 'Sample') %>% arrange(DaysPSO)
                                
pdf('Pseudobulk_TimeUMAP_v3.pdf')
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, color = PASC_status, shape = InfectionStages),
               size = 2) + theme_bw()
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, color = PASC, shape = InfectionStages),
               size = 2) + theme_bw()
                                
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, color = DaysPSO, shape = InfectionStages),
               size = 2) + theme_bw() +
            scale_colour_gradient(limits = c(0,30))
                                
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = PASC),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color =  PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()
                                
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()

ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ AgeGroups)
                                
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ PTID)

ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ PTID)
                                
ggplot(TimeUMAP2v2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ TrajectoryGroup)    
                                
ggplot(TimeUMAP2v2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = as.factor(PTID)),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
                                 geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages),
               size = 2.5) 
                                
                                
ggplot(TimeUMAP2v2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID
                  ),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5) 
                                
                                
ggplot(TimeUMAP2v2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID
                  color = batch_id),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5) 

dev.off()
      
                                
#### All Peaks UMAP
                                
STM <- readRDS('scMACS_SampleTileMatrix.RDS')

STM2 <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList = 'Positive')
sampleMat <- assays(STM2)[['CD16 Mono']]
sampleMat[is.na(sampleMat)] = 0

set.seed(20)
                                
umapAll <- as.data.frame(uwot::umap(2^t(sampleMat)-1))
colnames(umapAll ) = c('UMAP1', 'UMAP2')
umapAll$Sample <- rownames(umapAll)

umapAll<- inner_join(umapAll, as.data.frame(colData(STM2)), by = 'Sample') %>% arrange(DaysPSO)

umapAll <- umapAll%>% group_by(PTID) %>%
                    dplyr::mutate(Last = case_when(DaysPSO == max(DaysPSO) ~ 'Last',
  						    DaysPSO == min(DaysPSO) ~ 'First',
                                 TRUE ~ 'Other'))
umapAll$KMeans <- kmeans(umapAll [,c('UMAP1', 'UMAP2')], centers =3, nstart = 100)$cluster
                                
umapAll2 <- as.data.frame(uwot::umap(t(sampleMat)))
colnames(umapAll2 ) = c('UMAP1', 'UMAP2')
umapAll2$Sample <- rownames(umapAll2)
                                


umapAll2<- inner_join(umapAll2, as.data.frame(colData(STM2)), by = 'Sample') %>% arrange(DaysPSO)

umapAll2 <- umapAll2 %>% group_by(PTID) %>%
                    dplyr::mutate(Last = case_when(DaysPSO == max(DaysPSO) ~ 'Last',
  						    DaysPSO == min(DaysPSO) ~ 'First',
                                 TRUE ~ 'Other'))
umapAll2$KMeans <- kmeans(umapAll2 [,c('UMAP1', 'UMAP2')], centers =3, nstart = 100)$cluster

                                
pdf('Pseudobulk_UMAP2.pdf')

ggplot(umapAll, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = Last),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() + ggtitle('Raw data, all tiles')

ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5)  + ggtitle('Raw data, all tiles')
                                
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color =Last),
               size = 2.5)  + ggtitle('Raw data, all tiles')

ggplot(umapAll2, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = Last),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() + ggtitle('Log2 data, all tiles')

ggplot(umapAll2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5) + ggtitle('Log2 data, all tiles')
                                
ggplot(umapAll2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color =Last),
               size = 2.5)  + ggtitle('Log2 data, all tiles')
                                
dev.off()

### All peaks outside of the Time-related peakset
     
varPeaks <- read.csv('variance_decomposition_peaks.csv')                        
maxComp <- colnames(subVarPeaks)[apply(subVarPeaks, 1, which.max)]
varPeaks$MaxVar <- maxComp
Time_Regs <- filter(varPeaks, MaxVar == 'Time')
TimeRegs <- Time_Regs$Gene %>% sub("\\.",":",.) %>%
            sub("\\.","-", .) %>% MOCHA::StringsToGRanges()

sampleMat2 <- sampleMat[!rownames(sampleMat) %in% MOCHA::GRangesToString(TimeRegs),]

set.seed(1)
                                
umapNonTime <- as.data.frame(uwot::umap(t(sampleMat2)))
colnames(umapNonTime ) = c('UMAP1', 'UMAP2')
umapNonTime$Sample <- rownames(umapNonTime)

umapNonTime<- inner_join(umapNonTime, as.data.frame(colData(STM2)), by = 'Sample') %>% arrange(DaysPSO)

umapNonTime <- umapNonTime%>% group_by(PTID) %>%
                    dplyr::mutate(Last = case_when(DaysPSO == max(DaysPSO) ~ 'Last',
  						    DaysPSO == min(DaysPSO) ~ 'First',
                                 TRUE ~ 'Other'))
umapNonTime$KMeans <- kmeans(umapNonTime [,c('UMAP1', 'UMAP2')], centers =3, nstart = 100)$cluster
                                
pdf('Pseudobulk_NonTimeUMAP.pdf')

ggplot(umapNonTime, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = Last),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()

ggplot(umapNonTime) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5) 
                                
                                
dev.off()                                
                                
    
################################################################

##Identify groups of donors based on Normal and Abnormal Response

Abnormal = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)
Normal = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)


STM <- readRDS('scMACS_SampleTileMatrix.RDS')

meta1 <- colData(STM)
meta1$ResponseType =  dplyr::case_when(meta1$PTID %in% Normal ~ 'Normal',
                                       meta1$PTID %in% Abnormal ~ 'Abnormal',
                                      TRUE ~ 'Uninfected')
firstVisit <- as.data.frame(meta1) %>% dplyr::group_by(PTID) %>%
    dplyr::arrange(DaysPSO) %>%
    dplyr::slice_head(n=1) %>% 
    dplyr::filter(InfectionStages == 'Early Infection')

LastVisit <- as.data.frame(meta1) %>% dplyr::group_by(PTID) %>%
    dplyr::arrange(DaysPSO) %>%
    dplyr::slice_tail(n=1) %>% 
    dplyr::filter(InfectionStages == 'Recovery Phase')

meta1$VisitType <- dplyr::case_when(meta1$Sample %in% firstVisit$Sample ~ 'First',
                       meta1$Sample %in% LastVisit$Sample ~ 'Last',
                        meta1$InfectionStages == 'Uninfected' ~ 'Only',
                       TRUE ~ 'Other')

meta1$ResponseType_Visit <- 
    paste(meta1$ResponseType, meta1$VisitType, sep = '_')

colData(STM) <- meta1
saveRDS(STM, 'scMACS_SampleTileMatrix.RDS')

compGroups <- list(c('Normal_First', 'Uninfected_Only'),
                   c('Normal_Last', 'Uninfected_Only'),
                   c('Abnormal_First', 'Uninfected_Only'),
                   c('Abnormal_Last', 'Uninfected_Only'))
names(compGroups) <- c('Normal_First', 'Normal_Last',
                       'Abnormal_First','Abnormal_Last')
                                
timeTiles <- sub("\\.",":",Time_Regs$Gene) %>% sub("\\.", "-", .)

allDATsList <- lapply(compGroups, function(x){
    
                     getDifferentialAccessibleTiles(STM[rownames(STM) %in% timeTiles,] ,
                                            cellPopulation = 'CD16 Mono',
                                           groupColumn = 'ResponseType_Visit',
                                           foreground = x[1],
                                           background = x[2],
                                           signalThreshold = 0,
                                           minZeroDiff = 0,
                                           fdrToDisplay = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 15)
    
    })
#869 DATs - Normal First Visit
#181 DATs - Normal Last Visit
#0 DATs - Abnormal First Visit
#0 DATs - Abnormal Last Visit

hist(allDATsList[[1]]$P_value, breaks = 100)
hist(allDATsList[[2]]$P_value, breaks = 100)
hist(allDATsList[[3]]$P_value, breaks = 100)
hist(allDATsList[[4]]$P_value, breaks = 100)

hist(allDATsList[[1]]$Avg_Intensity_Control, breaks = 100, xlim = c(5,18))
hist(allDATsList[[1]]$Avg_Intensity_Case, breaks = 100, xlim = c(5,18))
hist(allDATsList[[2]]$Avg_Intensity_Control, breaks = 100, xlim = c(5,18))
hist(allDATsList[[2]]$Avg_Intensity_Case, breaks = 100, xlim = c(5,18))
hist(allDATsList[[3]]$Avg_Intensity_Control, breaks = 100, xlim = c(5,18))
hist(allDATsList[[3]]$Avg_Intensity_Case, breaks = 100, xlim = c(5,18))
hist(allDATsList[[4]]$Avg_Intensity_Control, breaks = 100, xlim = c(5,18))
hist(allDATsList[[4]]$Avg_Intensity_Case, breaks = 100, xlim = c(5,18))

hist(allDATsList[[1]]$Pct0_Control, breaks = 100)
hist(allDATsList[[1]]$Pct0_Case, breaks = 100)
hist(allDATsList[[2]]$Pct0_Control, breaks = 100)
hist(allDATsList[[2]]$Pct0_Case, breaks = 100)
hist(allDATsList[[3]]$Pct0_Control, breaks = 100)
hist(allDATsList[[3]]$Pct0_Case, breaks = 100)
hist(allDATsList[[4]]$Pct0_Control, breaks = 100)
hist(allDATsList[[4]]$Pct0_Case, breaks = 100)

dev.off()

allEnrichments <- lapply(allDATsList, function(x){
    
        foreList <- plyranges::filter(x, FDR <= 0.2)
        backList <- plyranges::filter(x, FDR > 0.2)
    
        MotifEnrichment(foreList, backList,
                        STM@metadata$Cisbp, numCores = 15)
    
    })
                                
                                
#################### Try again with just:
#  Acute response first vs Abnormal Response last 
#  Acute response first vs Normal Response last
                                
AcuteID = c(32220, 31874, 31924, 32251, 32196, 32209, 32140, 32054, 32255, 32415, 
            32131, 31945, 32416)
Abnormal_Last = c(32220, 32038, 31874, 32124, 32251, 32196, 32415)
Normal_Last = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)

meta1 <- colData(STM)

firstVisit <- as.data.frame(meta1) %>% dplyr::group_by(PTID) %>%
    dplyr::arrange(DaysPSO) %>%
    dplyr::slice_head(n=1) %>% 
    dplyr::filter(InfectionStages == 'Early Infection')

LastVisit <- as.data.frame(meta1) %>% dplyr::group_by(PTID) %>%
    dplyr::arrange(DaysPSO) %>%
    dplyr::slice_tail(n=1) %>% 
    dplyr::filter(InfectionStages == 'Recovery Phase')
                                
meta1$SampleType =  dplyr::case_when(meta1$PTID %in% Normal_Last & 
                                     meta1$Sample %in% LastVisit$Sample ~ 'NormalLast',
                                       meta1$PTID %in% Abnormal_Last &
                                     meta1$Sample %in% LastVisit$Sample ~ 'AbnormalLast',
                                       meta1$PTID %in% AcuteID &
                                     meta1$Sample %in% firstVisit$Sample ~ 'AcuteResponse',
                                     meta1$COVID_status == 'Negative' ~ 'Uninfected',
                                      TRUE ~ 'Other')

colData(STM) <- meta1
saveRDS(STM, 'scMACS_SampleTileMatrix.RDS')       


compGroups <- list(c('NormalLast',  'AcuteResponse'),
                  c('AbnormalLast',  'AcuteResponse'),
                   c('AbnormalLast', 'NormalLast'))

names(compGroups) <- c('Normal_v_Acute', 'Abnormal_v_Acute',
                       'Abnormal_v_Normal')
                                
timeTiles <- sub("\\.",":",Time_Regs$Gene) %>% sub("\\.", "-", .)

allDATsList <- lapply(compGroups, function(x){
    
                     getDifferentialAccessibleTiles(STM[rownames(STM) %in% timeTiles,] ,
                                            cellPopulation = 'CD16 Mono',
                                           groupColumn = 'SampleType',
                                           foreground = x[1],
                                           background = x[2],
                                           signalThreshold = 12,
                                           minZeroDiff = 0.5,
                                           fdrToDisplay = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 15)
    
    })
#1392 differential regions found at FDR 0.2
#744 differential regions found at FDR 0.2
#229 differential regions found at FDR 0.2        
                            
saveRDS(allDATsList, 'CD16_Trajectory_Part2_DATs.rds')




                                
##################################################################################################
                                
############### Re-running modelings
                                
##################################################################################################
                                
STM <- readRDS('scMACS_SampleTileMatrix.RDS')

STM2 <- subsetMOCHAObject(STM, subsetBy =  'COVID_status', groupList = 'Positive')
             
allAcc <- getCellPopMatrix(STM2, 'CD16 Mono')
set.seed(1)
allUMAP <- as.data.frame(uwot::umap(t(allAcc)))
colnames(allUMAP) <- c('UMAP1','UMAP2')
allUMAP$Sample <- rownames(allUMAP)
                                
subMat <-  colData(STM2) %>% as.data.frame()
                                
allUMAP2 <- dplyr::inner_join(allUMAP, subMat, by = 'Sample') %>% arrange(DaysPSO)
allUMAP2$KMeans <- kmeans(allUMAP2[,c('UMAP1', 'UMAP2')], centers =3, nstart = 100)$cluster
allUMAP2 <- allUMAP2  %>% dplyr::group_by(PTID) %>%
                    dplyr::mutate(Last = ifelse(DaysPSO == max(DaysPSO),
                                                'Last', 'NotLast')) %>%
                    dplyr::mutate(TrajectoryGroup = 
                                  ifelse(any(KMeans == 1 & Last == 'Last'),
                                         'Group 1', 'Group 2'))
pdf('UMAP_AllTiles.pdf')
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, color = PASC_status, shape = InfectionStages),
               size = 2) + theme_bw()

                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, color = DaysPSO, shape = InfectionStages),
               size = 2) + theme_bw() +
            scale_colour_gradient(limits = c(0,30))
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color =  PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw()

ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ AgeGroups)
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ PTID)

ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ PTID)
                                
ggplot(allUMAP2) + 
    geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages,
                  color = PASC_status),
               size = 2.5) +
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = PASC_status),
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    facet_wrap(~ TrajectoryGroup)    
                                
ggplot(allUMAP2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, 
                  color = as.factor(PTID)),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
                                 geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages),
               size = 2.5) 
                                
                                
ggplot(allUMAP2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID
                  ),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5) 
                                
                                
ggplot(allUMAP2) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID
                  color = batch_id),
              alpha = 0.30, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "closed")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5)          
                                
                                
dev.off()
                                
### Generate UMAP 
                                
32416
42409
32415
32255
31207
6       1 32245
7       1 32209
8       1 32054
9       1 32131
10      1 32140
11      1 31945
## Older groups                               
#Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
#Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)
                                
#Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
#Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)
                                
STM_G1 <- subsetMOCHAObject(STM2[rowRanges(STM2)$tileType ==  'Promoter',], 
                            subsetBy =  'PTID', groupList = Group1)
STM_G2 <- subsetMOCHAObject(STM2[rowRanges(STM2)$tileType ==  'Promoter',],
                            subsetBy =  'PTID', groupList = Group2)
                                
G1_lmem <- linearModeling(STM_G1,formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID),
                          CellType = 'CD16 Mono', numCores = 30)
                                
saveRDS(G1_lmem, 'Group1_LinearModel.rds')                             
G1_lmem <- readRDS('Group1_LinearModel.rds')
                                
#Group 2
G2_lmem <- linearModeling(STM_G2,
                          formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID),
                          CellType = 'CD16 Mono', numCores = 30)
                               
saveRDS(G2_lmem, 'Group2_LinearModel.rds')
G2_lmem <- readRDS('Group2_LinearModel.rds')
                                
G1_pval <- parallel::mclapply(seq_along(G1_lmem), function(x) {
    tryCatch(expr = {summary(G1_lmem[[x]])$coefficients[,5]},
             error = function(cond){NA})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
saveRDS(G1_pval, 'Group1_LinearModel_PValues.rds')  
G1_pval <- readRDS('Group1_LinearModel_PValues.rds')  
                                
G2_pval <- parallel::mclapply(seq_along(G2_lmem), function(x) {
    tryCatch(expr = {summary(G2_lmem[[x]])$coefficients[,5]},
             error = function(cond){NA})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
saveRDS(G2_pval, 'Group2_LinearModel_PValues.rds')                                
G2_pval <- readRDS('Group2_LinearModel_PValues.rds') 

### Pull out slopes
G1_slopes <- parallel::mclapply(seq_along(G1_lmem), function(x) {
    tryCatch(expr = {summary(G1_lmem[[x]])$coefficients[,1]},
             error = function(cond){NA})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
G2_slopes <- parallel::mclapply(seq_along(G2_lmem), function(x) {
    tryCatch(expr = {summary(G2_lmem[[x]])$coefficients[,1]},
             error = function(cond){NA})
}, mc.cores = 20) %>% do.call('rbind',.)

saveRDS(G1_slopes, 'Group1_LinearModel_Slopes.rds') 
saveRDS(G2_slopes, 'Group2_LinearModel_Slopes.rds') 
                                
G1_slopes <- readRDS('Group1_LinearModel_Slopes.rds') 
G2_slopes <- readRDS('Group2_LinearModel_Slopes.rds') 
                                
G1_combined <- data.frame(Tiles = MOCHA::GRangesToString(rowRanges(STM_G1)),
                          P_values = as.data.frame(G1_pval)$days_since_symptoms, 
                          Slopes = as.data.frame(G1_slopes)$days_since_symptoms,
                         Group = rep('Group1', length(rowRanges(STM_G1))))  
                          
G2_combined <- data.frame(Tiles = MOCHA::GRangesToString(rowRanges(STM_G2)),
                          P_values = as.data.frame(G2_pval)$days_since_symptoms, 
                          Slopes = as.data.frame(G2_slopes)$days_since_symptoms,
                         Group = rep('Group2', length(rowRanges(STM_G2))))                          
volcTiles <- rbind(G1_combined, G2_combined) %>% 
                dplyr::mutate(Significant = ifelse(P_values < 0.05, 'Significant', 'Not')) %>% 
                dplyr::mutate(Significant = ifelse(is.na(Significant),'Not' ,Significant ))
                                
pdf('VolcanoPlots.pdf')
         
ggplot(G1_combined, aes(x = Slopes, y = -log10(P_values))) + 
                                geom_point() + theme_minimal()+
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')

ggplot(G2_combined, aes(x = Slopes, y = -log10(P_values))) + 
                                geom_point() + theme_minimal()+
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')
                                
ggplot(volcTiles, aes(x = Slopes, y = -log10(P_values), 
                     color = Significant)) + 
            geom_point() + theme_minimal() +
            facet_wrap(~ Group, ncol = 1 ) +
            scale_color_manual(values=c('Significant' = 'red', 'Not' = 'black') ) +
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')
            
                                
ggplot(volcTiles, aes(x = Slopes, y = -log10(P_values),)) + 
            geom_point() + theme_minimal() +
            facet_wrap(~ Group, ncol = 2)+
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')
 
                                
dev.off()
                                

PromoGR_G1 <-  rowRanges(STM_G1)
mcols(PromoGR_G1) <- cbind(mcols(PromoGR_G1), G1_pval)
                                
saveRDS(PromoGR_G1, 'Group1_CD16Mono_TimeTiles.rds')
PromoGR_G1 <- readRDS('Group1_CD16Mono_TimeTiles.rds')
                                
                                
PromoGR_G2 <-  rowRanges(STM_G2)
G2_pval1 <-  G2_pval
G1_pval1 <-  G1_pval
colnames(G1_pval1) <- paste('Pvalue', colnames(G1_pval1), sep = '_')
colnames(G2_pval1) <- paste('Pvalue', colnames(G2_pval1), sep = '_')
                                
mcols(PromoGR_G1) <- cbind(mcols(PromoGR_G1), G1_pval1)
mcols(PromoGR_G2) <- cbind(mcols(PromoGR_G2), G2_pval1)
                                
colnames(G1_slopes) <- paste('Coefficient', colnames(G1_slopes), sep = '_')
colnames(G2_slopes) <- paste('Coefficient', colnames(G2_slopes), sep = '_')
                                
mcols(PromoGR_G1) <- cbind(mcols(PromoGR_G1), G1_slopes)
mcols(PromoGR_G2) <- cbind(mcols(PromoGR_G2), G2_slopes)
                                
saveRDS(PromoGR_G2, 'Group2_CD16Mono_TimeTiles.rds')
PromoGR_G2 <- readRDS('Group2_CD16Mono_TimeTiles.rds')    
                                
write.csv(as.data.frame(plyranges::filter(PromoGR_G2, Pvalue_days_since_symptoms < 0.05)),
          'SupplementalTable_Group2_Significant_PromoterTiles.csv')
write.csv(as.data.frame(plyranges::filter(PromoGR_G1,  Pvalue_days_since_symptoms < 0.05)),
          'SupplementalTable_Group1_Significant_PromoterTiles.csv')
                                

              
accMat2 <- getCellPopMatrix(STM_G2, 'CD16 Mono')
subAccMat2 <- accMat2[rownames(accMat2)  %in% names(PromoGR_G2)[PromoGR_G2$days_since_symptoms < 0.05],]    
column_ha2 = HeatmapAnnotation(Days =  anno_barplot(colData(STM_G2)$days_since_symptoms[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)]),
                        Donor = as.character(colData(STM_G2)$PTID[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)]))
                                

     
    
accMat1 <- getCellPopMatrix(STM_G1, 'CD16 Mono')
subAccMat1 <- accMat1[rownames(accMat1)  %in% names(PromoGR_G1)[PromoGR_G1$days_since_symptoms < 0.05],]   
                                
column_ha1 = HeatmapAnnotation(Days =  anno_barplot(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)]),
                        Donor = as.character(colData(STM_G1)$PTID[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)]))
                                
col_fun = circlize::colorRamp2(c(0, 12, 20), c("blue", "white", "red"))
col_fun2 = circlize::colorRamp2(c(0, 5, 20), c("blue", "white", "red"))

pdf('SubGroup_CD16Mono_Heatmaps.pdf')
                                
ComplexHeatmap::Heatmap(subAccMat1, top_annotation = column_ha1, 
                        column_order = order(as.numeric(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)])), 
                        column_title = 'Group 1 Infection Time Course',
                        #col = col_fun,
                        heatmap_legend_param = list(col_fun = col_fun),
                          show_column_names = FALSE,
                        show_row_names= FALSE)   
                                
ComplexHeatmap::Heatmap(subAccMat2, top_annotation = column_ha2, 
                        column_order = order(as.numeric(colData(STM_G2)$days_since_symptoms[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)])), 
                            column_title = 'Group 2 Infection Time Course',
                        #col = col_fun,
                        show_column_names = FALSE,
                        show_row_names= FALSE)   

ComplexHeatmap::Heatmap(subAccMat1, top_annotation = column_ha1, 
                        column_order = order(as.numeric(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)])), 
                        column_title = 'Group 1 Infection Time Course',
                        col = col_fun,
                        heatmap_legend_param = list(col_fun = col_fun),
                          show_column_names = FALSE,
                        show_row_names= FALSE)   
                                
ComplexHeatmap::Heatmap(subAccMat2, top_annotation = column_ha2, 
                        column_order = order(as.numeric(colData(STM_G2)$days_since_symptoms[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)])), 
                            column_title = 'Group 2 Infection Time Course',
                        col = col_fun,
                        show_column_names = FALSE,
                        show_row_names= FALSE)   
                                
ComplexHeatmap::Heatmap(subAccMat1, top_annotation = column_ha1, 
                        column_order = order(as.numeric(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)])), 
                        column_title = 'Group 1 Infection Time Course',
                        col = col_fun2,
                        heatmap_legend_param = list(col_fun = col_fun),
                          show_column_names = FALSE,
                        show_row_names= FALSE)   
                                
ComplexHeatmap::Heatmap(subAccMat2, top_annotation = column_ha2, 
                        column_order = order(as.numeric(colData(STM_G2)$days_since_symptoms[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)])), 
                            column_title = 'Group 2 Infection Time Course',
                        col = col_fun2,
                        show_column_names = FALSE,
                        show_row_names= FALSE)   
                                

dev.off()
                                
                                
## Pathway enrichment:
                                
G1_geneset <- PromoGR_G1$Gene[PromoGR_G1$Pvalue_days_since_symptoms < 0.05 & !is.na(PromoGR_G1$Pvalue_days_since_symptoms)]  %>% strsplit(., split = ", ") %>% unlist()   %>%
                                unique()
 
G2_geneset <- PromoGR_G2$Gene[PromoGR_G2$Pvalue_days_since_symptoms < 0.05 & !is.na(PromoGR_G2$Pvalue_days_since_symptoms)]  %>% strsplit(., split = ", ") %>% unlist()   %>%
                                unique()
write.table(G1_geneset, 'Group1_Geneset.txt', row.names= F)
write.table(G2_geneset, 'Group2_Geneset.txt', row.names= F)
allGenes <- rowRanges(STM)$Gene %>% strsplit(., split = ", ") %>% unlist()   %>%
                                unique()
#allMotifs <- getPositions(CD16s, 'CD16LongCISBPMotif')
enrichDataBaseList <- WebGestaltR::listGeneSet()
                                
                   
GE_G1 <- lapply(enrichDataBaseList$name[c(2,7,8,9,10, 59, 60, 61)], function(x){
    
        simplifiedORA(x, G1_geneset , allGenes)
    
    })
                                
GE_G2 <- lapply(enrichDataBaseList$name[c(2,7,8,9,10, 59, 60, 61)], function(x){
    
        simplifiedORA(x, G2_geneset, allGenes)
    
    })
                                

    
CD16s <- loadArchRProject('CD16sLong')

allGenes2 <- getGenes(CD16s)
#allMotifs <- getPositions(CD16s, 'CD16LongCISBPMotif')
                                
trajGS2 <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory2', 
                         useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
CD16s <- addPeakSet(CD16s, rowRanges(STM), force = TRUE)
CD16s <- addPeakMatrix(CD16s)
saveArchRProject(CD16s)
trajPeak2 <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory2', 
                         useMatrix = "PeakMatrix")

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)                              
trajGS2Mat <- plotTrajectoryHeatmap(trajGS2, returnMatrix = TRUE)
trajPeakMat <- plotTrajectoryHeatmap(trajPeak2, returnMatrix = TRUE)
                                
trajTiles <- annotateTiles(MOCHA::StringsToGRanges(gsub("_","-",rownames(trajPeakMat))),
                          TxDb =TxDb.Hsapiens.UCSC.hg38.refGene, Org = org.Hs.eg.db)
trajPromoterTiles <- unique(unlist(
                                lapply(plyranges::filter(trajTiles, tileType == 'Promoter')$Gene,
                                       function(x){strsplit(x, split= ', ')})))      

GE_Trajectory <- lapply(enrichDataBaseList$name[c(2,7,8,9,10, 59, 60, 61)], function(x){
    
          simplifiedORA(x, unique(gsub('.*:','', rownames(trajGS2Mat))),  
                        unique(as.character(allGenes2$symbol)))
    
    })
                                
GE_TileTrajectory <- lapply(enrichDataBaseList$name[c(2,7,8,9,10, 59, 60, 61)], function(x){
    
          simplifiedORA(x,trajPromoterTiles, allGenes)
    
    })
            
                   
saveRDS(list(GE_G1, GE_G2, GE_Trajectory, GE_TileTrajectory),'Fig5_AllPathwayEnrichments.rds')

                   
                                
GE_Group1 <- WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = enrichDataBaseList$name[9],
                         interestGene = G1_geneset, 
                        interestGeneType ="genesymbol",
                        referenceGene =   allGenes, 
                                    fdrThr = 1,
                            referenceGeneType= "genesymbol")     %>% as.data.frame() %>%
                        dplyr::mutate(GroupType = 'Group1') 
                                
GE_Group2 <- WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = enrichDataBaseList$name[9],
                         interestGene = G2_geneset, 
                        interestGeneType ="genesymbol",
                        referenceGene =   allGenes, 
                                    fdrThr = 1,
                            referenceGeneType= "genesymbol")   %>% as.data.frame() %>%
                        dplyr::mutate(GroupType = 'Group2')
## Plot out Monocle3 trajectory, with a brief comparison to other groups
GE_Traj <- WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = enrichDataBaseList$name[9],
                         interestGene = unique(gsub('.*:','', rownames(trajGS2Mat))),
                        interestGeneType ="genesymbol",
                        referenceGene = unique(allGenes2$symbol), 
                                    fdrThr = 1,
                            referenceGeneType= "genesymbol")  %>% as.data.frame() %>%
                        dplyr::mutate(GroupType = 'Monocle3')
write.csv(GE_Traj[GE_Traj$FDR < 0.05, ], 'SupplementalFigure_Monocle3_Trajectory_GeneScores.csv')
                                  
GE_Traj2 <- WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = enrichDataBaseList$name[9],
                         interestGene =trajPromoterTiles,
                        interestGeneType ="genesymbol",
                        referenceGene = allGenes, 
                                    fdrThr = 1,
                            referenceGeneType= "genesymbol")  %>% as.data.frame() %>%
                        dplyr::mutate(GroupType = 'Monocle3')
                                
write.csv(GE_Traj2[GE_Traj2$FDR < 0.05, ], 'SupplementalFigure_Monocle3_Trajectory_AlteredPromoters.csv')
                                
library(UpSetR)
                        
upsetInput <- list('Group2' = GE_Group2$description[GE_Group2$FDR < 0.05],
                'Group1' = GE_Group1$description[GE_Group1$FDR < 0.05],
                'Monocle3' = GE_Traj$description[GE_Traj$FDR < 0.05])
                                
upsetInput2 <- list('Group2' = GE_Group2$description[GE_Group2$FDR < 0.05],
                'Group1' = GE_Group1$description[GE_Group1$FDR < 0.05],
                'Monocle3' = GE_Traj2$description[GE_Traj2$FDR < 0.05])
                                
pdf('Figure5_MonocleTrajectory.pdf')
                                
plotTrajectory(CD16s, embedding = "normUMAP", trajectory = 'Trajectory2', name = 'Trajectory2')
plotTrajectoryHeatmap(trajGS2)
                                
ggplot(as.data.frame(GE_Traj[GE_Traj$FDR < 0.05,]), aes(x = description, y = -log10(FDR))) +
        geom_col() + 
       xlab('-log10 of FDR') + ylab('Reactime Pathway') + theme_minimal() +
                                theme(axis.text.x = element_text(angle=-45,vjust = 0.5, hjust=0.5 ))
                                
UpSetR::upset(fromList(upsetInput), order.by = "freq", sets.bar.color=c("blue","orange","red"),
            main.bar.color = c("blue","orange","red", 'black') )
                                
plotTrajectoryHeatmap(trajPeak2)
                                
ggplot(as.data.frame(GE_Traj2[GE_Traj2$FDR < 0.05,]), aes(x = description, y = -log10(FDR))) +
        geom_col() + 
       xlab('-log10 of FDR') + ylab('Reactime Pathway') + theme_minimal() +
                                theme(axis.text.x = element_text(angle=-45,vjust = 0.5, hjust=0.5 ))
                                
UpSetR::upset(fromList(upsetInput2), order.by = "freq", sets.bar.color=c("blue","orange","red"),
            main.bar.color = c("blue",'black' ,"purple", "orange",'brown', "red") )
                                
dev.off()
                                
                                
                                
## Putting Monocle3 together with others is a mess. Do it separately.             
toKeep <- c(GE_Group2$description[GE_Group2$FDR < 0.05],
            GE_Group1$description[GE_Group1$FDR < 0.05])
   
## Merge results into one heatmap
allPathwaysEnrich <- rbind(dplyr::filter(as.data.frame(GE_Group1), description %in% toKeep), 
                     dplyr::filter(as.data.frame(GE_Group2), description %in% toKeep)) %>%
                tidyr::pivot_wider(id_cols = 'description', names_from = 'GroupType',
                                   values_from = "enrichmentRatio")
allPathwaysFDR <- rbind(dplyr::filter(as.data.frame(GE_Group1), description %in% toKeep), 
                     dplyr::filter(as.data.frame(GE_Group2), description %in% toKeep)) %>%
                 tidyr::pivot_wider(id_cols = 'description', names_from = 'GroupType',
                                   values_from = "FDR")      
                                
allPathways <- rbind(dplyr::filter(as.data.frame(GE_Group1), description %in% toKeep), 
                     dplyr::filter(as.data.frame(GE_Group2), description %in% toKeep)) 

                             
allPathEnMat <- as.matrix(allPathwaysEnrich[,c('Group1','Group2')])
allPathFDRMat <- as.matrix(allPathwaysFDR[,c('Group1','Group2')] < 0.05)
rownames(allPathEnMat) <- allPathwaysEnrich$description
rownames(allPathFDRMat) <- allPathwaysFDR$description
all(rownames(allPathEnMat) == rownames(allPathFDRMat))
                                
allPathEnMat[is.na(allPathEnMat)] = 0
allPathFDRMat[is.na(allPathFDRMat)] = 0
#TRUE
                                
rowA <- rowAnnotation(
    foo = anno_text(stringr::str_wrap(allPathwaysEnrich$description, 35),
                    just = "left",  gp = gpar(fontsize =5),
                    location = unit(1, "npc")
    ), legend_param = list(direction = "horizontal"))
           
                
                                
pathwayAnnot = dplyr::case_when(
    grepl('TLR|MyD88|TGF|NTRK1|stimuli|NLR|Interleukin|catenin|NOD',
          rownames(allPathFDRMat)) ~ 'Immune Response',
    grepl('UPR|DNA|Apoptot|Repair|Ub|Senescence|Death|Deub', 
          rownames(allPathFDRMat)) ~ 'DNA Damage & Apoptosis',
    grepl('Chromatin|pigenetic',  rownames(allPathFDRMat)) ~ 'Epigenetics',
    grepl('Mitot|G2|Separation|Cell Cycle',  rownames(allPathFDRMat)) ~ 'Cell Cycle',
    TRUE ~ "Other")
                                
pathwayAnnot2 = dplyr::case_when(
       grepl('TLR|MyD88|TGF|NTRK1|stimuli|NLR|Interleukin|catenin|NOD',
             allPathways$description) ~ 'Immune Response',
    grepl('UPR|DNA|Apoptot|Repair|Ub|Senescence|Death|Deub',
          allPathways$description) ~ 'DNA Damage & Apoptosis',
    grepl('Chromatin|pigenetic',   allPathways$description) ~ 'Epigenetics',
    grepl('Signaling|GTPase',   allPathways$description) ~ 'General Signaling',
    grepl('Neuro',   allPathways$description) ~ 'Neurodegeneration',
    grepl('Mitot|G2|Separation|Cell Cycle',   allPathways$description) ~ 'Cell Cycle',
    TRUE ~ "Other")
                                
                                
allPathwaysFDR$pathwayAnnot <- pathwayAnnot
allPathwaysEnrich$pathwayAnnot <- pathwayAnnot  
allPathways$pathwayAnnot <- pathwayAnnot2

col_fun1 = colorRamp2(c(1, 2.5), c("white", "red"))
                                
pathTypeColors <- c('green','brown', 'purple', 'yellow', 'black')
names(pathTypeColors) <- unique(pathwayAnnot)
                                
rowA2 <- rowAnnotation(
    Group1 = as.data.frame(allPathwaysFDR[,c('Group1')]) < 0.05,
    Group2 = as.data.frame(allPathwaysFDR[,c('Group2')]) < 0.05,
    Type = pathwayAnnot,
    col = list(Group1 = c('FALSE' = 'white', 'TRUE' = 'blue'), 
               Group2 = c('FALSE' = 'white', 'TRUE' = 'blue'), 
               Type = pathTypeColors) ,
    
    legend_param = list(direction = "horizontal",direction = "horizontal",direction = "horizontal" ))
                                
                                
#FDRH <- ComplexHeatmap::Heatmap(allPathFDRMat, name = '-log FDR',
#                                width = unit(2, "cm"),
#                               cluster_columns=  FALSE,
#                                show_row_dend = FALSE,
#                                #col = col_fun2,
#                                heatmap_legend_param = list(labels =c('False', 'True'), 
#                                                             title = "FDR < 0.05", 
#                                                             legend_gp = gpar(fill = 1:2),
#                                                            direction = "horizontal"),
#                                right_annotation = rowAnnotation('Type' = pathwayAnnot),
#                                row_names_gp = gpar(fontsize = 12),
#                                  column_title= 'Significance')
     
#H_list <- EnH + FDRH 
                                
EnH <- ComplexHeatmap::Heatmap(allPathEnMat, name = 'Enrichment Ratio',
                               cluster_columns=  FALSE,
                                width = unit(2, "cm"),
                                #show_row_names = FALSE, 
                               #right_annotation = rowA,
                              show_row_dend = FALSE,
                               col = col_fun1,
                            right_annotation = rowA2,
                              heatmap_legend_param = list(direction = "horizontal"),
                             row_names_gp = gpar(fontsize = 12),
                               column_title= 'Enrichment')
                                
                                
pdf('SubGroup_CD16Mono_Pathway_Heatmaps.pdf')
 
EnH
#FDRH
                                
draw(EnH, heatmap_legend_side = "bottom")
                   
ggplot(allPathways, aes(x = GroupType, y = str_wrap(description, 30), 
                        color = enrichmentRatio,
                       size = -log10(FDR))) +
                   geom_point() +  
                   scale_size_continuous(range = c(0, 4)) + theme_minimal()
                   
                   
ggplot(allPathways, aes(x = -log10(FDR), y = enrichmentRatio, 
                        color = pathwayAnnot)) +
                   geom_point() +  
                   scale_size_continuous(range = c(0, 4)) + theme_minimal() +
                   facet_wrap(~GroupType)


                                
dev.off()
                                
sumSigPathways <- allPathwaysFDR %>% 
                    tidyr::pivot_longer(cols = c('Group1', 'Group2'), names_to = 'Group',
                                        values_to = 'PValues') %>%
                        dplyr::filter(PValues < 0.05) %>% #dplyr::mutate(PValues = PValues < 0.05) %>%
                    dplyr::group_by(Group, pathwayAnnot) %>% dplyr::summarize(PathNumb = dplyr::n()) %>%
                    dplyr::mutate(pathwayAnnot = factor(pathwayAnnot, 
                            levels = c('Immune Response', 'Cell Cycle', 'Epigenetics', 
                                       'DNA Damage & Apoptosis', 'Other')))
            
                                
pdf('SubGroup_CD16Mono_Pathway_BarChart.pdf')
         
ggplot(sumSigPathways, aes(x = Group, y = PathNumb,
                           fill = pathwayAnnot)) +
                        scale_fill_manual(values = pathTypeColors) + geom_col(position = 'dodge') +
                    theme_minimal() + xlab('Response Group') + ylab('Number of Enriched Pathways')                  
                                
dev.off()
                   
#################################################################
                   
#### 
library(chromVAR)
library(chromVARmotifs)
data(human_pwms_v2)

STM <- addMotifSet(STM, pwms = human_pwms_v2,
                      w = 7, returnObj = TRUE, motifSetName = "Motifs")
                   
saveRDS(STM, 'scMACS_SampleTileMatrix.RDS')
                   

### Generate a heatmap:
#Each row is one gene, each column is one sample, ranked by days since symptoms. 
## Donor ID information as bar at the top. 
                                
### Then generate a heatmap of pathway enrichment for all pathways detected in Group 1 and Group 2. 
                                
        
###### Functions used. 
                        
simplifiedORA <- function(database, foreground, background){
    
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                            referenceGeneType= "genesymbol")
    
    }
