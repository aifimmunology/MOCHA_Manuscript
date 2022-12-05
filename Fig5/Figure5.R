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

CD16s <- subsetArchRProject(MonoDCE, 
                            cells = getCellNames(MonoDCE)[
                                MonoDCE$predictedGroup_Co2 == 'CD16 Mono'],
                            outputDirectory = 'CD16sLong', dropCells = TRUE)
saveArchRProject(CD16s)
CD16s <- loadArchRProject('CD16sLong')

CD16s <- addPeakSet(CD16s, rowRanges(STM), force = TRUE)
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

#### All Peaks UMAP
                                
STM <- readRDS('scMACS_SampleTileMatrix.RDS')
                   
varPeaks <- read.csv('variance_decomposition_peaks.csv') 

subVarPeaks <- varPeaks[,c('PTID','Time','Age','Sex','Residual')]
maxComp <- colnames(subVarPeaks)[apply(subVarPeaks, 1, which.max)]
varPeaks$MaxVar <- maxComp
                                
Time_Regs <- filter(varPeaks, MaxVar == 'Time')
TimeRegs <- Time_Regs$Gene %>% sub("\\.",":",.) %>%
            sub("\\.","-", .) %>% StringsToGRanges() 

STM2 <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList = 'Positive')
sampleMat <- assays(STM2)[['CD16 Mono']]
sampleMat[is.na(sampleMat)] = 0

set.seed(2)
                                
umapAll <- as.data.frame(uwot::umap(2^t(sampleMat)-1), min_dist = 0.025)
colnames(umapAll ) = c('UMAP1', 'UMAP2')
umapAll$Sample <- rownames(umapAll)

umapAll<- inner_join(umapAll, as.data.frame(colData(STM2)), by = 'Sample') %>% arrange(DaysPSO)

umapAll <- umapAll%>% group_by(PTID) %>%
                    dplyr::mutate(Last = case_when(DaysPSO == max(DaysPSO) ~ 'Last',
  						    DaysPSO == min(DaysPSO) ~ 'First',
                                 TRUE ~ 'Other'))
umapAll$KMeans <- kmeans(umapAll [,c('UMAP1', 'UMAP2')], centers =3, nstart = 200)$cluster
                                
write.csv(umapAll, 'Pseudobulk_UMAP_alltiles.csv')
umapAll <- read.csv('Pseudobulk_UMAP_alltiles.csv')
                   
### Identify groups:
                   
#Group 1 --> starts in Cluster 3 and ends in Cluster 1
Group1 <- umapAll %>% dplyr::group_by(PTID) %>% 
    filter(any(Last == 'Last' & KMeans == 1)) %>%
                   select(PTID) %>% unlist() %>% unique()
#42409 32416 32054 31207 32255 32415 32131 32209 32140 31945 32245
Group2 <- umapAll %>% dplyr::group_by(PTID) %>% 
    filter(any(Last == 'Last' & KMeans != 1)) %>%
                   select(PTID) %>% unlist() %>% unique()
#32038 32220 32251 31924 31874 32196 32124

## Older groups     #32415 is different                          
Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)

#newer groups
#Group1 = c(42409, 32416, 32054, 31207, 32255, 32415, 32131, 32209, 32140, 31945, 32245)
#Group2 = c(32038, 32220, 32251, 31924, 31874, 32196, 32124)
                   

              
umapAll <- dplyr::mutate(umapAll, Groups = ifelse(PTID %in% Group1, 'Group1', 'Group2'))
write.csv(umapAll, 'Pseudobulk_UMAP_alltiles.csv')
umapAll <- read.csv('Pseudobulk_UMAP_alltiles.csv')                   
                   
pdf('Pseudobulk_UMAP.pdf')

ggplot(umapAll, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = Last),
               size = 2.5) +
           scale_color_manual(values = c('First' = '#2CB04A', 'Last' = '#E9282A', 'Other' = 'black'))+
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID), alpha = 0.40,
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() + ggtitle('Raw data, all tiles') 


ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = as.factor(KMeans)),
               size = 2.5)  + ggtitle('Raw data, all tiles') 
                                
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color =Last),
               size = 2.5)  + ggtitle('Raw data, all tiles') +
      scale_color_manual(values = c('First' = '#2CB04A', 'Last' = '#E9282A', 'Other' = 'black'))
                            
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, color = PASC_status),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages, color = PASC_status),
               size = 2.5)  + ggtitle('Raw data, all tiles') 
                   
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, color = PASC_status),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages),
               size = 2.5)  + ggtitle('Raw data, all tiles') +
        facet_wrap(~ PTID)
                   
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID, color = PASC_status),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
            geom_point(aes(x = UMAP1, y= UMAP2, shape = InfectionStages),
               size = 2.5)  + ggtitle('Raw data, all tiles') +
        facet_wrap(~ Groups)
                   
ggplot(umapAll) + 
    geom_path(aes(x = UMAP1, y = UMAP2, group = PTID),
              alpha = 0.40, linejoin = "mitre",
             arrow = arrow(angle = 30, 
                           length = unit(0.25, "inches"),
      ends = "last", type = "open")) + theme_bw() +
    geom_point(aes(x = UMAP1, y= UMAP2,
            size = factor(InfectionStages, 
                           levels = c('Recovery Phase', 'Late Infection', 'Early Infection')), 
                color = PASC_status, shape = InfectionStages))  + 
                   ggtitle('Raw data, all tiles') +
        facet_wrap(~ Groups) + theme(legend.position = 'bottom') +
                scale_size_discrete(range = c(2,4,6))
                                
                                
dev.off()
                                

                                            

    
################################################################

##Identify groups of donors based on Normal and Abnormal Response

STM <- readRDS('scMACS_SampleTileMatrix.RDS')
                   
Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)

meta1 <- colData(STM)
meta1$ResponseType =  dplyr::case_when(meta1$PTID %in% Group1 ~ 'Group1',
                                       meta1$PTID %in% Group2 ~ 'Group2',
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

compGroups <- list(c('Group1_First', 'Group1_Last'),
                   c('Group2_First', 'Group2_Last'),
                   c('Group1_Last', 'Group2_Last'),
                   c('Group1_First', 'Group2_First'))
names(compGroups) <- c('Group 1 Trajectory', 'Group 2 Trajectory',
                       'Last Comparison','First Comparison')
                            
allDAPsList <- lapply(compGroups, function(x){
    
                     getDifferentialAccessibleTiles(STM[rowRanges(STM)$tileType ==  'Promoter',],
                                            cellPopulation = 'CD16 Mono',
                                           groupColumn = 'ResponseType_Visit',
                                           foreground = x[1],
                                           background = x[2],
                                           fdrToDisplay = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 15)
    
    })
                   
saveRDS(allDAPsList, 'CD16_Fig5_DifferentialTiles.rds')
                   
allDAPsList <- readRDS('CD16_Fig5_DifferentialTiles.rds')
                   
#For all tiles with old groups:
#1828 DATs - 'Group 1 Trajectory'
#0 DATs - 'Group 2 Trajectory'
#2791 DATs - 'Last Comparison'
#0  DATs - 'First Comparison'
                   
#For all tiles:
#2132 DATs - 'Group 1 Trajectory'
#0 DATs - 'Group 2 Trajectory'
#0 DATs - 'Last Comparison'
#0  DATs - 'First Comparison'
    
#For promoters: 
#2080 DATs - 'Group 1 Trajectory'
#0 DATs - 'Group 2 Trajectory'
#3132 DATs - 'Last Comparison'
#0  DATs - 'First Comparison'                   

allGenes <- rowRanges(STM)$Gene %>% strsplit(., split = ", ") %>% unlist()   %>%
                                unique()
#allMotifs <- getPositions(CD16s, 'CD16LongCISBPMotif')
enrichDataBaseList <- WebGestaltR::listGeneSet()
                                
                   
Group1_TrajectoryPath <- lapply(enrichDataBaseList$name[c(2,7,8,9,10, 59, 60, 61)], function(x){
    
        simplifiedORA(x, G1_geneset , allGenes)
    
    })
                   
                   
##Identify co-accessible tiles to the Differentially Accessible Promoters. 
STM2 <- subsetMOCHAObject(STM, subsetBy = 'PTID', groupList = append(Group1, Group2))
allCoAccList <- lapply(allDAPsList[c(1,3)], function(x){
    
                    regions <- plyranges::filter(x, FDR <= 0.2)
                    tmpCoAcc <- getCoAccessibleLinks(STM, cellPopulation = 'CD16 Mono', regions = regions, 
                                 windowSize = 1 * 10^6, numCores = 20, 
                                 ZI = TRUE, verbose = FALSE)
                    filterCoAccessibleLinks(tmpCoAcc)
    })
                   
saveRDS(allCoAccList , 'CD16_Fig5_coAccTiles.rds')
                   
allDAPs <- unlist(lapply(allDAPsList[c(1,3)], function(y){ plyranges::filter(y, FDR <= 0.2)})) %>% 
                   as(., 'GRangesList') %>% unlist()
backPromoters <- plyranges::filter(rowRanges(STM), tileType == 'Promoter') %>%
                   plyranges::filter_by_non_overlaps(allDAPs)
                   
backgroundInter <- getCoAccessibleLinks(STM2, cellPopulation = 'CD16 Mono', 
                        regions = backPromoters, 
                                 windowSize = 1 * 10^6, numCores = 10, 
                                 ZI = TRUE, verbose = FALSE)
saveRDS(backgroundInter, 'CD16_Fig5_background.rds')                   

allEnrichments <- lapply(allDATsList, function(x){
    
        foreList <- plyranges::filter(x, FDR <= 0.2)
        backList <- plyranges::filter(x, FDR > 0.2)
    
        MotifEnrichment(foreList, backList,
                        STM@metadata$Cisbp, numCores = 15)
    
    })
###############################################################################################
                   
## Summarize PASC by group                   
dataSummary <- colData(STM) %>% as.data.frame() %>% dplyr::filter(PASC_status != 'Controls') %>%
    dplyr::mutate(Group = ifelse(PTID %in% Group1, 'Group1','Group2')) %>%
    dplyr::mutate(PASC = ifelse(grepl('PASC', PASC_status), 'PASC', 'Recovered')) %>% 
                   group_by(PTID) %>% slice_head(n=1) %>%
                   dplyr::group_by(Group, PASC) %>% dplyr::summarize(Num = dplyr::n()) 
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
                   
pdf('Fig5_PieCharts.pdf')
ggplot(dataSummary[dataSummary$Group =='Group1',], aes(x = "", y = Num, fill = PASC)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste(round(Num/10,2)*100,'%',sep='')),
            position = position_stack(vjust = 0.5),
            size=12) +
  coord_polar(theta = "y") +blank_theme+
  theme(legend.position='none',
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_brewer()+ggtitle('Group1: Pasc Distribution')

ggplot(dataSummary[dataSummary$Group =='Group2',], aes(x = "", y = Num, fill = PASC)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste(round(Num/8,2)*100,'%',sep='')),
            position = position_stack(vjust = 0.5),
            size=12) +
  coord_polar(theta = "y") +blank_theme+
  theme(legend.position='none',
        title=element_text(size=22, hjust = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_brewer()+ggtitle('Group2: Pasc Distribution')                   
dev.off()               
    
                                
#################################################### Linear Modeling
                           
Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)                              

STM <- readRDS('scMACS_SampleTileMatrix.RDS')
                                      
STM_G1 <- subsetMOCHAObject(STM[rowRanges(STM)$tileType ==  'Promoter',], 
                            subsetBy =  'PTID', groupList = Group1)
STM_G2 <- subsetMOCHAObject(STM[rowRanges(STM)$tileType ==  'Promoter',],
                            subsetBy =  'PTID', groupList = Group2)

mat1 <- MOCHA::getCellPopMatrix(STM,'CD16 Mono',NAtoZero = FALSE)
hist(rowSums(!is.na(mat1))/dim(mat1)[2] )
                   
library(lmerTest)
G1_lmem <- linearModeling(STM_G1,formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID), NAtoZero = FALSE, 
                          CellType = 'CD16 Mono',  numCores = 20)
G1_lmem2 <- linearModeling(STM_G1,formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID), NAtoZero = TRUE, 
                          CellType = 'CD16 Mono',  numCores = 20)                   
                   
saveRDS(G1_lmem, 'Group1_LinearModel.rds')
saveRDS(G1_lmem2, 'With_Zeros/Group1_LinearModel.rds')
G1_lmem <- readRDS('Group1_LinearModel.rds')
                   
                         
#Group 2: With and without zeros. 
G2_lmem <- linearModeling(STM_G2,
                formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID),
                          CellType = 'CD16 Mono',  NAtoZero = FALSE, numCores = 20)
G2_lmem2 <- linearModeling(STM_G2,
                formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID),
                          CellType = 'CD16 Mono',  NAtoZero = TRUE, numCores = 20                          
saveRDS(G2_lmem, 'Group2_LinearModel.rds')
saveRDS(G2_lmem2, 'With_Zeros/Group2_LinearModel.rds')     
                           
G2_lmem <- readRDS('Group2_LinearModel.rds')
G2_lmem2 <- readRDS('With_Zeros/Group2_LinearModel.rds')
                   

#If there were too many NAs for a given region, then the linear model will return an error.
#Let's give it a dummy variable to preserve information. 
NullBindP <- matrix(nrow = 1, ncol = 4, data = 1)                    
G1_pval <- parallel::mclapply(seq_along(G1_lmem), function(x) {
    tryCatch(expr = {summary(G1_lmem[[x]])$coefficients[,5]},
             error = function(cond){NullBindP})
}, mc.cores = 20) %>% do.call('rbind',.)
                           
G1_pval2 <- parallel::mclapply(seq_along(G1_lmem2), function(x) {
    tryCatch(expr = {summary(G1_lmem2[[x]])$coefficients[,5]},
             error = function(cond){NullBindP})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
saveRDS(G1_pval, 'Group1_LinearModel_PValues.rds')  
saveRDS(G1_pval2, 'With_Zeros/Group1_LinearModel_PValues.rds')  
                           
G1_pval <- readRDS('Group1_LinearModel_PValues.rds')  
G1_pval2 <- readRDS('With_Zeros/Group1_LinearModel_PValues.rds')  
                                
G2_pval <- parallel::mclapply(seq_along(G2_lmem), function(x) {
    tryCatch(expr = {summary(G2_lmem[[x]])$coefficients[,5]},
             error = function(cond){NullBindP})
}, mc.cores = 20) %>% do.call('rbind',.)
G2_pval2 <- parallel::mclapply(seq_along(G2_lmem2), function(x) {
    tryCatch(expr = {summary(G2_lmem2[[x]])$coefficients[,5]},
             error = function(cond){NullBindP})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
saveRDS(G2_pval, 'Group2_LinearModel_PValues.rds')  
saveRDS(G2_pval2, 'With_Zeros/Group2_LinearModel_PValues.rds')  
                           
G2_pval <- readRDS('Group2_LinearModel_PValues.rds')  
G2_pval2 <- readRDS('With_Zeros/Group2_LinearModel_PValues.rds')  

### Pull out slopes. Same dummy variable for regions that were too sparse. 
NullBindSlope <- matrix(nrow = 1, ncol = 4, data = 0)   
G1_slopes <- parallel::mclapply(seq_along(G1_lmem), function(x) {
    tryCatch(expr = {summary(G1_lmem[[x]])$coefficients[,1]},
             error = function(cond){NullBindSlope})
}, mc.cores = 20) %>% do.call('rbind',.)
                           
G1_slopes2 <- parallel::mclapply(seq_along(G1_lmem2), function(x) {
    tryCatch(expr = {summary(G1_lmem2[[x]])$coefficients[,1]},
             error = function(cond){NullBindSlope})
}, mc.cores = 20) %>% do.call('rbind',.)
                                
G2_slopes <- parallel::mclapply(seq_along(G2_lmem), function(x) {
    tryCatch(expr = {summary(G2_lmem[[x]])$coefficients[,1]},
             error = function(cond){NullBindSlope})
}, mc.cores = 20) %>% do.call('rbind',.)
G2_Slopes2 <- parallel::mclapply(seq_along(G2_lmem2), function(x) {
    tryCatch(expr = {summary(G2_lmem2[[x]])$coefficients[,1]},
             error = function(cond){NullBindSlope})
}, mc.cores = 20) %>% do.call('rbind',.)

saveRDS(G1_slopes, 'Group1_LinearModel_Slopes.rds') 
saveRDS(G2_slopes, 'Group2_LinearModel_Slopes.rds') 
                    
saveRDS(G1_Slopes2, 'With_Zeros/Group1_LinearModel_Slopes.rds') 
saveRDS(G2_Slopes2, 'With_Zeros/Group2_LinearModel_Slopes.rds') 
                                
G1_slopes <- readRDS('Group1_LinearModel_Slopes.rds') 
G2_slopes <- readRDS('Group2_LinearModel_Slopes.rds')
                           
                                
G1_Slopes2 <- readRDS('With_Zeros/Group1_LinearModel_Slopes.rds') 
G2_Slopes2 <- readRDS('With_Zeros/Group2_LinearModel_Slopes.rds')
                                
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
                   
slopeRange <- max(abs(range(volcTiles$Slopes[volcTiles$Significant == 'Significant'])))
                                
pdf('VolcanoPlots.pdf')
         
ggplot(G1_combined, aes(x = Slopes, y = -log10(P_values))) + 
                                geom_point() + theme_minimal()+
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope') +
            + xlim(-slopeRange,slopeRange) 

ggplot(G2_combined, aes(x = Slopes, y = -log10(P_values))) + 
                                geom_point() + theme_minimal()+
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')+
            xlim(-slopeRange,slopeRange) 
                                
ggplot(volcTiles, aes(x = Slopes, y = -log10(P_values), 
                     color = Significant)) + 
            geom_point() + theme_minimal() +
            facet_wrap(~ Group, ncol = 1 ) +
            scale_color_manual(values=c('Significant' = 'red', 'Not' = 'black') ) +
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')+
            xlim(-slopeRange,slopeRange) 
            
                                
ggplot(volcTiles, aes(x = Slopes, y = -log10(P_values), 
                     color = Significant)) + 
            geom_point() + theme_minimal() +
            facet_wrap(~ Group, nrow = 1 ) +
            scale_color_manual(values=c('Significant' = 'red', 'Not' = 'black') ) +
            xlab('Slope Accessibility Change over Disease Time Course') +
            ylab('-log of P-values for Slope')+
            xlim(-slopeRange,slopeRange) 
            
 
                                
dev.off()
                                

PromoGR_G1 <-  rowRanges(STM_G1)
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
                   
saveRDS(PromoGR_G1, 'Group1_CD16Mono_TimeTiles.rds')                                
saveRDS(PromoGR_G2, 'Group2_CD16Mono_TimeTiles.rds')

PromoGR_G1 <- readRDS('Group1_CD16Mono_TimeTiles.rds')
PromoGR_G2 <- readRDS('Group2_CD16Mono_TimeTiles.rds') 
                                
write.csv(as.data.frame(plyranges::filter(PromoGR_G2, Pvalue_days_since_symptoms < 0.05)),
          'SupplementalTable_Group2_Significant_PromoterTiles.csv')
write.csv(as.data.frame(plyranges::filter(PromoGR_G1,  Pvalue_days_since_symptoms < 0.05)),
          'SupplementalTable_Group1_Significant_PromoterTiles.csv')
                                

                   
                   
library(ComplexHeatmap)
accMat2 <- getCellPopMatrix(STM_G2, 'CD16 Mono')
subAccMat2 <- accMat2[rownames(accMat2)  %in% names(PromoGR_G2)[PromoGR_G2$Pvalue_days_since_symptoms < 0.05],]    
column_ha2 = HeatmapAnnotation(Days =  anno_barplot(colData(STM_G2)$days_since_symptoms[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)]),
                        Donor = as.character(colData(STM_G2)$PTID[
    match(colnames(subAccMat2), colData(STM_G2)$Sample)]))
                                
#### Shift in promoter-accessibility per time point.
            
specGenes <- c('FPR2',#'GPR32',
               'CMKLR1', 'LTB4R', 'GPR18')
specPromo <- which(rowRanges(STM)$tileType ==  'Promoter' &
                   grepl(paste(specGenes, collapse = '|'), rowRanges(STM)$Gene))
cellMat <- getCellPopMatrix(STM[specPromo,], 'CD16 Mono') %>% as.data.frame() %>%
                    dplyr::mutate(Tile = rownames(.)) %>%
                   dplyr::mutate(Gene = rowRanges(STM)$Gene[specPromo]) %>%
                   tidyr::pivot_longer(cols = colnames(STM),
                                       names_to = 'Sample',
                                       values_to = 'Log2_Accessibility') %>%
                   dplyr::full_join(as.data.frame(colData(STM)), by = 'Sample') %>%
                   dplyr::filter(COVID_status == 'Positive') %>%
                   dplyr::mutate(DonorGroup = ifelse(PTID %in% Group1, 'Group1', 'Group2')) %>%
                   dplyr::mutate(TileGroup = paste(PTID, Tile, sep = ' ')) %>%
                   dplyr::mutate(Gene = 
                                 case_when(grepl(specGenes[1], Gene) ~ specGenes[1],
                                           grepl(specGenes[2], Gene) ~ specGenes[2],
                                           grepl(specGenes[3], Gene) ~ specGenes[3],
                                           grepl(specGenes[4], Gene) ~ specGenes[4],
                                           grepl(specGenes[5], Gene) ~ specGenes[5])) %>%
                   dplyr::mutate(GeneTile = paste(Gene, Tile, sep = "_"))  %>%
                   dplyr::mutate(NA_Accessibility = ifelse(Log2_Accessibility == 0, NA, 
                                                           Log2_Accessibility))

                   
                   
pdf('InflammationGenes_Trajectory.pdf')
                   
lapply(1:4, function(x){

    ggplot(filter(cellMat, Gene == specGenes[x]), 
           aes(x = days_since_symptoms, y = Log2_Accessibility, 
                    group = TileGroup,color = DonorGroup)) + 
    geom_line() +  facet_grid(rows = vars(DonorGroup), cols = vars(Tile)) + 
    theme_bw()  + ggtitle(specGenes[x])
                   
    
})
                   
lapply(1:4, function(x){
                   
    ggplot(filter(cellMat, Gene == specGenes[x]), 
           aes(x = days_since_symptoms, y = NA_Accessibility, 
                    group = TileGroup,color = DonorGroup)) + 
    geom_line() + facet_grid(rows = vars(DonorGroup), cols = vars(Tile)) +
    theme_bw() + ggtitle(specGenes[x])
    
})
    
                   
dev.off()
                   
pdf('InflammationGenes_Trajectory2.pdf')
                   
lapply(1:4, function(x){

    ggplot(filter(cellMat, Gene == specGenes[x]), 
           aes(x = days_since_symptoms, y = Log2_Accessibility, 
                    group = TileGroup,color = DonorGroup)) + 
    geom_line() +  facet_wrap(~ Tile) + 
    theme_bw()  + ggtitle(specGenes[x])
                   
    
})
                   
lapply(1:4, function(x){
                   
    ggplot(filter(cellMat, Gene == specGenes[x]), 
           aes(x = days_since_symptoms, y = NA_Accessibility, 
                    group = TileGroup,color = DonorGroup)) + 
    geom_line() + facet_wrap(~ Tile) + 
    theme_bw() + ggtitle(specGenes[x])
    
})
    
                   
dev.off()
           
                                       
                   
            
                   
     
#### Heatmap objects
accMat1 <- getCellPopMatrix(STM_G1, 'CD16 Mono')
subAccMat1 <- accMat1[rownames(accMat1)  %in% names(PromoGR_G1)[PromoGR_G1$Pvalue_days_since_symptoms < 0.05],]   
                                
column_ha1 = HeatmapAnnotation(Days =  anno_barplot(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)]),
                        Donor = as.character(colData(STM_G1)$PTID[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)]))
                                
col_fun = circlize::colorRamp2(c(0, 12, 20), c("blue", "white", "red"))
col_fun2 = circlize::colorRamp2(c(0, 5, 20), c("blue", "white", "red"))
        
bothSlopes <- lapply(list(PromoGR_G1,PromoGR_G2), function(x){ 
    plyranges::filter(x, Pvalue_days_since_symptoms <= 0.05) })
bothMats <- list(subAccMat1, subAccMat2)
                     
## Split into two heatmaps: Slope < 0.025
lowSlopeMat <- lapply(seq_along(bothSlopes), function(x){
        lowSlope <- plyranges::filter(bothSlopes[[x]],
                                      abs(Coefficient_days_since_symptoms) <= 0.025)
        bothMats[[x]][rownames(bothMats[[x]]) %in% names(lowSlope),]
    })
                   
highSlopeMat <- lapply(seq_along(bothSlopes), function(x){
      highSlope <- plyranges::filter(bothSlopes[[x]],
                                      abs(Coefficient_days_since_symptoms) > 0.025)     
      bothMats[[x]][rownames(bothMats[[x]]) %in% names(highSlope),]
    })
      
bothGroups <- list(STM_G1, STM_G2)
bothAnnotations <- list(column_ha1, column_ha2)
                   

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
                   
ComplexHeatmap::Heatmap(2^subAccMat1-1, top_annotation = column_ha1, 
                        column_order = order(as.numeric(colData(STM_G1)$days_since_symptoms[
    match(colnames(subAccMat1), colData(STM_G1)$Sample)])), 
                        column_title = 'Group 1 Infection Time Course',
                        #col = col_fun,
                        heatmap_legend_param = list(col_fun = col_fun),
                          show_column_names = FALSE,
                        show_row_names= FALSE)   
#Do original lambda1, not log-transformed
ComplexHeatmap::Heatmap(2^subAccMat2-1, top_annotation = column_ha2, 
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
                   
lapply(seq_along(lowSlopeMat), function(x){
ComplexHeatmap::Heatmap(lowSlopeMat[[x]], top_annotation = bothAnnotations[[x]], 
                        column_order = order(as.numeric(colData(bothGroups[[x]])$days_since_symptoms[
    match(colnames(bothMats[[x]]), colData(bothGroups[[x]])$Sample)])), 
                            column_title = 'Low Slope Infection Time Course',
                        show_column_names = FALSE,
                        show_row_names= FALSE)      
    
  })

lapply(seq_along(highSlopeMat), function(x){
        ComplexHeatmap::Heatmap(highSlopeMat[[x]], top_annotation = bothAnnotations[[x]],
                        column_order = order(as.numeric(colData(bothGroups[[x]])$days_since_symptoms[
                match(colnames(bothMats[[x]]), colData(bothGroups[[x]])$Sample)])), 
                            column_title = 'High Slope Infection Time Course',
                        show_column_names = FALSE,
                        show_row_names= FALSE)      
    
  })
                   
lapply(seq_along(highSlopeMat), function(x){
ComplexHeatmap::Heatmap(lowSlopeMat[[x]], top_annotation = bothAnnotations[[x]], 
                        column_order = order(as.numeric(colData(bothGroups[[x]])$days_since_symptoms[
    match(colnames(bothMats[[x]]), colData(bothGroups[[x]])$Sample)])), 
                            column_title = 'Low Slope Infection Time Course',
                        show_column_names = FALSE,
                        show_row_names= FALSE)      %v%                    
 ComplexHeatmap::Heatmap(highSlopeMat[[x]], top_annotation = bothAnnotations[[x]],
                        column_order = order(as.numeric(colData(bothGroups[[x]])$days_since_symptoms[
                match(colnames(bothMats[[x]]), colData(bothGroups[[x]])$Sample)])), 
                            column_title = 'High Slope Infection Time Course',
                        show_column_names = FALSE,
                        show_row_names= FALSE) 
})                
dev.off()
                                
                                
## Pathway enrichment:
                                
G1_geneset <- PromoGR_G1$Gene[PromoGR_G1$Pvalue_days_since_symptoms < 0.05]  %>% strsplit(., split = ", ") %>% unlist()   %>%
                                unique()
 
G2_geneset <- PromoGR_G2$Gene[PromoGR_G2$Pvalue_days_since_symptoms< 0.05 & !is.na(PromoGR_G2$Pvalue_days_since_symptoms)]  %>% strsplit(., split = ", ") %>% 
                   unlist()   %>% unique()
                   
write.table(G1_geneset, 'Group1_Geneset.txt', row.names= F)
write.table(G2_geneset, 'Group2_Geneset.txt', row.names= F)
                   
allGenes <- plyranges::filter(rowRanges(STM), !is.na(Gene)) %>% .$Gene %>% 
                   strsplit(., split = ", ") %>% unlist()  %>%
                                unique()
                   

txList <- suppressWarnings(GenomicFeatures::transcriptsBy(TxDb.Hsapiens.UCSC.hg38.refGene,
                                                          by = ("gene")))
  names(txList) <- suppressWarnings(AnnotationDbi::mapIds(org.Hs.eg.db, names(txList), "SYMBOL", "ENTREZID"))
                   
                   
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
                                   values_from = "FDR")  %>%
                    dplyr::mutate(Group1 = ifelse(is.na(Group1), 1, Group1),
                                 Group2 = ifelse(is.na(Group2), 1, Group2))
                                
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
           
     
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                   
                   
pathwayAnnot = unlist(lapply(rownames(allPathFDRMat), 
                             function(x) findTrunk(x, Hierarch, ReactID)))
#28, 29, 23
pathwayAnnot[grepl('external stimuli', pathwayAnnot)] = 'Cellular responses to stimuli'
                                                
pathwayAnnot2 =  unlist(lapply(allPathways$description, 
                             function(x) findTrunk(x, Hierarch, ReactID)))          
pathwayAnnot2[grepl('external stimuli', pathwayAnnot2)] = 'Cellular responses to stimuli'                                
allPathwaysFDR$pathwayAnnot <- pathwayAnnot
allPathwaysEnrich$pathwayAnnot <- pathwayAnnot  
allPathways$pathwayAnnot <- pathwayAnnot2

col_fun1 = colorRamp2(c(1, 2.5), c("white", "red"))
                                
all13Colors = c('darkslategray3','azure4','bisque2','brown1','brown4',
                'cadetblue1','chartreuse2','coral1', 'cornflowerblue',
                'darkmagenta', 'darkorange3', 'darkgreen', 'darkslategrey')
names(all13Colors) <- unique(pathwayAnnot)
                                
rowA2 <- rowAnnotation(
    Group1 = as.data.frame(allPathwaysFDR[,c('Group1')]) < 0.05,
    Group2 = as.data.frame(allPathwaysFDR[,c('Group2')]) < 0.05,
    Type = pathwayAnnot,
    col = list(Group1 = c('FALSE' = 'white', 'TRUE' = 'blue'), 
               Group2 = c('FALSE' = 'white', 'TRUE' = 'blue'),
              Type = all13Colors))
                                
                                
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
                             row_names_gp = gpar(fontsize = 6),
                               column_title= 'Enrichment')
                               
                                
sumSigPathways <- allPathwaysFDR %>% 
                    tidyr::pivot_longer(cols = c('Group1', 'Group2'), names_to = 'Group',
                                        values_to = 'PValues') %>%
                        dplyr::filter(PValues < 0.05) %>% #dplyr::mutate(PValues = PValues < 0.05) %>%
                    dplyr::group_by(Group, pathwayAnnot) %>% dplyr::summarize(PathNumb = dplyr::n()) 
     
library(RColorBrewer)
                               
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
                    theme_minimal() +
                   facet_wrap(~ GroupType)

ggplot(sumSigPathways, aes(x = Group, y = PathNumb,
                            fill = pathwayAnnot)) +
         scale_fill_manual(values = all13Colors) + 
         geom_col(position = 'dodge') + 
    theme_minimal() + xlab('Response Group') + 
    # scale_color_brewer(palette = "Spectral") +
                               ylab('Number of Enriched Pathways')    
                                
dev.off()
                               
                               
#### Create two heatmaps: For with and without zeros, along with the matching chart.
### Need to coordiante colorscheme. 
STM <- readRDS('scMACS_SampleTileMatrix.RDS')
PromoGR_G1 <- readRDS('Group1_CD16Mono_TimeTiles.rds')
PromoGR_G2 <- readRDS('Group2_CD16Mono_TimeTiles.rds') 
PromoGR_G1v2 <- readRDS('With_Zeros/Group1_CD16Mono_TimeTiles.rds')
PromoGR_G2v2 <- readRDS('With_Zeros/Group2_CD16Mono_TimeTiles.rds') 
   
allGenes <- plyranges::filter(rowRanges(STM), !is.na(Gene)) %>% .$Gene %>% 
                   strsplit(., split = ", ") %>% unlist()  %>%
                                unique()
                   
geneSets <- lapply(list(PromoGR_G1, PromoGR_G2,PromoGR_G1v2,PromoGR_G2v2), function(x){
    
        x$Gene[x$Pvalue_days_since_symptoms < 0.05]  %>% 
                strsplit(., split = ", ") %>% unlist()   %>%  unique()
                    
    })
enrichDataBaseList <- WebGestaltR::listGeneSet()    
                               
GroupsNames <- c('Group1_NoZeros', 'Group2_NoZeros', 'Group1_Zeros',
                               'Group2_Zeros')
allEnrichments <- lapply(seq_along(geneSets), function(y){
    
        WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = enrichDataBaseList$name[9],
                         interestGene = geneSets[[y]] , 
                        interestGeneType ="genesymbol",
                        referenceGene =   allGenes, 
                                    fdrThr = 1,
                            referenceGeneType= "genesymbol")   %>% as.data.frame()  %>%
            dplyr::mutate(Group = GroupsNames[y])
    
    })
                                   
combEnrich <- data.table::rbindlist(allEnrichments, use.names = TRUE) %>% as.data.frame() %>%
                    group_by(description) %>% filter(any(FDR < 0.05))           
                    
                   
    
                               
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                               
                   
pathwayAnnot_all = unlist(lapply(combEnrich$description, 
                             function(x) findTrunk(x, Hierarch, ReactID))) 
pathwayAnnot_all[grepl('external stimuli', pathwayAnnot_all)] = 
                                 'Cellular responses to stimuli'  
                                 
combEnrich$Type <- pathwayAnnot_all
                                 
allPathwaySum <- combEnrich %>% 
                    dplyr::select(description, FDR, enrichmentRatio, Group, Type) %>% 
                        dplyr::filter(FDR < 0.05) %>% 
                    dplyr::mutate(Group = gsub("_Z"," With Z",gsub("_No"," Without ",Group))) %>% 
                    dplyr::group_by(Group, Type) %>% dplyr::summarize(PathNumb = dplyr::n()) 
                             
allColors = c('darkslategray3','azure4','bisque2','brown1','brown4',
                'cadetblue1','chartreuse2','coral1', 'cornflowerblue',
                'darkmagenta', 'darkorange3', 'darkgreen', 'darkslategrey','darkred')
names(allColors) <- unique(pathwayAnnot_all) 
                                 
matrixList <- lapply(list('_Zeros','_NoZeros'), function(x){
                
            tmp <- filter(combEnrich, grepl(x,Group)) %>% 
                dplyr::mutate(Group = gsub('_.*', '', Group)) %>%
                tidyr::pivot_wider(id_cols = c(description, Type), names_from = 'Group',
                    values_from = 'enrichmentRatio')
            tmp$Group1[is.na(tmp$Group1)] = 0
            tmp$Group2[is.na(tmp$Group2)] = 0  
            tmp
                                   
})
                                                             
binaryList <- lapply(list('_Zeros','_NoZeros'), function(x){
                
            tmp <- filter(combEnrich, grepl(x,Group)) %>% 
                dplyr::mutate(Group = gsub('_.*', '', Group)) %>%
                dplyr::mutate(FDR = FDR < 0.05) %>% 
                tidyr::pivot_wider(id_cols = c(description, Type), names_from = 'Group',
                    values_from = 'FDR') 
            tmp$Group1[is.na(tmp$Group1)] = FALSE
            tmp$Group2[is.na(tmp$Group2)] = FALSE
            tmp     
})


                              
rowAList <- lapply(1:2, function(x){
    
    rowAnnotation(
        Group1 = as.data.frame(binaryList[[x]][,c('Group1')]) < 0.05,
        Group2 = as.data.frame(binaryList[[x]][,c('Group2')]) < 0.05,
        Type = binaryList[[x]]$Type,
        col = list(Group1 = c('FALSE' = 'white', 'TRUE' = 'blue'), 
               Group2 = c('FALSE' = 'white', 'TRUE' = 'blue'),
              Type = allColors))
    
})

allMat <-  combEnrich %>% 
                dplyr::mutate(Group =
                              gsub("_Z"," With Z",gsub("_No"," Without ",Group))) %>%
                tidyr::pivot_wider(id_cols = description,
                                   names_from = 'Group',
                    values_from = 'enrichmentRatio') %>% as.data.frame()
allMat[is.na(allMat)] = 0     
rownames(allMat) <- allMat$description
allBinary <- combEnrich %>% 
                dplyr::mutate(Group =
                              gsub("_Z"," With Z",gsub("_No"," Without ",Group))) %>%
                dplyr::mutate(FDR = FDR < 0.05) %>% 
                tidyr::pivot_wider(id_cols = c(description, Type),
                                   names_from = 'Group',
                    values_from = 'FDR') 
allBinary[is.na(allBinary)] = FALSE
       
GroupsNames2 = gsub("_Z"," With Z",GroupsNames) %>% 
                                 gsub("_No"," Without ", .)
allGroups <- lapply(GroupsNames2, function(x){
                allBinary[,x] < 0.05
            })
names(allGroups) <- GroupsNames2

grep('Mito',allBinary$description[allBinary$'Group1 Without Zeros'], value= TRUE)
                                 
allRows <-  rowAnnotation(
        "Group1 Without Zeros" = unlist(allBinary[,"Group1 Without Zeros"]),
        "Group1 With Zeros" = unlist(allBinary[,"Group1 With Zeros" ]),
        "Group2 Without Zeros" = unlist(allBinary[,"Group2 Without Zeros"]),
        "Group2 With Zeros" = unlist(allBinary[,"Group2 With Zeros"]),
        Type = allBinary$Type,
        col = list("Group1 Without Zeros"   = c('FALSE' = 'white', 'TRUE' = 'blue'), 
                 "Group1 With Zeros" =  c('FALSE' = 'white', 'TRUE' = 'blue'), 
            "Group2 Without Zeros" =  c('FALSE' = 'white', 'TRUE' = 'blue'), 
            "Group2 With Zeros" =  c('FALSE' = 'white', 'TRUE' = 'blue'), 
              Type = allColors))
                            
                                 
library(UpSetR)
                                 
sigPathways <- lapply(GroupsNames, function(x){
    
                dplyr::filter(combEnrich, Group == x & FDR < 0.05) %>%
                dplyr::select(description) %>% unlist()
    
    })
names(sigPathways) <- gsub("_Z"," With Z",GroupsNames) %>% 
                                 gsub("_No"," Without ", .)
names(geneSets)  <- gsub("_Z"," With Z",GroupsNames) %>% 
                                 gsub("_No"," Without ", .)  
colors12 <- c(RColorBrewer::brewer.pal(n = 6, name = 'Set2'),
              RColorBrewer::brewer.pal(n = 6, name = 'Dark2'))   
colors15 <- c(RColorBrewer::brewer.pal(n = 7, name = 'Set2'),
              RColorBrewer::brewer.pal(n = 8, name = 'Dark2'))


col_fun1 = circlize::colorRamp2(c(1, 2.5), c("white", "red"))
                                 
pdf('Fig5_Modeling_Both_Pathways.pdf')
  
lapply(1:2, function(x){
    enrichMatrix <- as.matrix(matrixList[[x]][,c(3,4)])
    rownames(enrichMatrix) <- unlist(matrixList[[x]][,1])
    ComplexHeatmap::Heatmap(enrichMatrix, name = 'Enrichment Ratio',
                       cluster_columns=  FALSE,
                        width = unit(2, "cm"),
                        #show_row_names = FALSE, 
                       #right_annotation = rowA,
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                       col = col_fun1,
                    right_annotation = rowAList[[x]],
                      heatmap_legend_param = list(direction = "horizontal"),
                     row_names_gp = gpar(fontsize = 6),
                       column_title= 'Enrichment')
    
})
                              
ComplexHeatmap::Heatmap(as.matrix(allMat[,-1]), name = 'Enrichment Ratio',
                       cluster_columns=  TRUE,
                        width = unit(2, "cm"),
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                       col = col_fun1,
                    right_annotation = allRows,
                      heatmap_legend_param = list(direction = "horizontal"),
                     row_names_gp = gpar(fontsize = 6),
                       column_title= 'Enrichment')
                                 
UpSetR::upset(fromList(sigPathways), order.by = "freq", main.bar.color = colors12 ) 
UpSetR::upset(fromList(geneSets), order.by = "freq", main.bar.color=  colors15)
                                 
ggplot(allPathwaySum, aes(x = Group, y = PathNumb,
                            fill = Type)) +
         scale_fill_manual(values = allColors) + 
         geom_col(position = 'dodge') + 
    theme_minimal() + xlab('Response Group') + 
                               ylab('Number of Enriched Pathways')                                    
                                 
dev.off()
                    
#################################################################
                   
#### 
library(chromVAR)
library(chromVARmotifs)
data(human_pwms_v2)

STM <- addMotifSet(STM, pwms = human_pwms_v2,
                      w = 7, returnObj = TRUE, motifSetName = "Motifs")
                   
saveRDS(STM, 'scMACS_SampleTileMatrix.RDS')
              
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(tidyverse)
library(plyranges)
              
assayList <- getCellPopMatrix(STM, 'CD16 Mono')
        
Obj1 <- SummarizedExperiment(
    assays = list('counts' = assayList),
    colData = colData(STM),
    rowRanges = rowRanges(STM),
    metadata = STM@metadata
)

CisbpAnno <- chromVAR::getAnnotationslinear(STM@metadata$Motifs,
                                      rowRanges = rowRanges(Obj1))

Obj1 <- chromVAR::addGCBias(Obj1, genome = BSgenome.Hsapiens.UCSC.hg38)
backPeaks <- getBackgroundPeaks(Obj1) 
                               
dev <- chromVAR::computeDeviations(object = Obj1, 
                         annotations = CisbpAnno)
saveRDS(dev, 'STM_ChromVAR.rds')
                    
Group1 = c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409)
Group2 = c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)   
              
dev_G1 <- subsetDev(dev, 
                            subsetBy =  'PTID', groupList = Group1)
dev_G2 <- subsetDev(dev,
                            subsetBy =  'PTID', groupList = Group2)
 
library(lmerTest)
G1_dev_lmem <- modelDeviations(dev_G1 ,formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID), 
                          type = 'z', numCores = 5)
G2_dev_lmem <- modelDeviations(dev_G2 ,formula = exp ~ Age + sex_at_birth + days_since_symptoms + (1|PTID), 
                          type = 'z', numCores = 5)             

saveRDS(list(G1_dev_lmem, G2_dev_lmem), 'Fig5_ChromVar_Modeling.rds')
   
allPval <- lapply(list(G1_dev_lmem, G2_dev_lmem), function(y){
                        
    pullLmerCoef(y, 5, numCores = 5)  %>% do.call('rbind',.)
        
    })
                               
allSlopes <- lapply(list(G1_dev_lmem, G2_dev_lmem), function(y){
                        
    pullLmerCoef(y, 1, numCores = 5) %>% do.call('rbind',.)
        
    })
                                
saveRDS(list(allSlopes, allPval ), 'Fig5_ChromVar_Pval_slopes.rds')                         
                               
combinedDim <- lapply(seq_along(allSlopes), function(y){
        
             pval <- allPval[[y]][,'days_since_symptoms']
                
             slope1 <- allSlopes[[y]][,'days_since_symptoms']
            
             data.frame('TF' = rownames(dev), 'Pval' = pval, 'Slope' = slope1, 
                       Group = y) %>% 
                dplyr::mutate(Label = ifelse(pval <= 0.05, 'Significant', 'Insignificant'))
            
    })
         
library(ggrepel)
                               
pdf('Fig5_ChromVar_modeling.pdf')
         
    ggplot(combinedDim[[1]], aes(x = Slope, y = -log10(Pval), color = Label)) + 
                               geom_point() +
                ggtitle('Group 1 Motifs') + 
        geom_label_repel(data =  dplyr::filter(combinedDim[[1]], Pval < 0.05), 
                         aes(label = TF))
                               
    ggplot(combinedDim[[2]], aes(x = Slope, y =  -log10(Pval), color = Label)) + 
                               geom_point() +
                ggtitle('Group 2 Motifs') + 
        geom_label_repel(data = dplyr::filter(combinedDim[[2]], Pval < 0.05), 
                         aes(label = TF))
                         
dev.off()

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

subsetDev <- function(Object,
                              subsetBy,
                              groupList){

    sampleData <- SummarizedExperiment::colData(Object)
    
      if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
    }
  }


    keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]]))]
    return(Object[, keep])
    
}
                  
modelDeviations <- function(Obj, formula, type = 'z', numCores = 1){
    

    meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))

    mat1 <- SummarizedExperiment::assays(Obj)[[type]] 

    meta <- meta1[meta1$Sample %in% colnames(mat1),]
    
    suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
        function(x) {
            df <- data.frame(exp = as.numeric(mat1[x, ]), 
                meta, stringsAsFactors = FALSE)
            lmerTest::lmer(formula = formula, data = df)
        }, cl = numCores), classes = "message")

    return(lmem_res)

    
}
              
findTrunk <- function(pathway, TreeStructure, IDMatrix){

    specID <- IDMatrix[match(pathway,IDMatrix[,2]),1]

    if(!(any(TreeStructure[,2] %in% specID))){

        return(pathway)

    }

    treeID <- TreeStructure[which(TreeStructure[,2] %in% specID)[1],]

    while(any(TreeStructure$V2 %in% treeID[,1])){

        treeID <- TreeStructure[which(TreeStructure$V2 == treeID[,1])[1],]


    }

    TrunkDescription <- IDMatrix[match(treeID[,1],IDMatrix[,1]),2]

    return(TrunkDescription)

}


###
                        
pullLmerCoef <- function(LmerList, Coeff_Num, numCores = 1){
    
    parallel::mclapply(seq_along(LmerList), function(x) {
        summary(LmerList[[x]])$coefficients[,Coeff_Num]
        }, mc.cores = numCores )
    
}