###############################
### Add just time-related peaks
    
    
varPeaks <- read.csv('variance_decomposition_peaks.csv') 

subVarPeaks <- varPeaks[,c('PTID','Time','Age','Sex','Residual')]
maxComp <- colnames(subVarPeaks)[apply(subVarPeaks, 1, which.max)]
varPeaks$MaxVar <- maxComp
                                
Time_Regs <- filter(varPeaks, MaxVar == 'Time')
TimeRegs <- Time_Regs$Gene %>% sub("\\.",":",.) %>%
            sub("\\.","-", .) %>% StringsToGRanges()    

CD16s <- subsetArchRProject(MonoDCE, cells = getCellNames(MonoDCE)[MonoDCE$predictedGroup_Co2 == 'CD16 Mono'],
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
                   
CD16s <- addSlingShotTrajectories(CD16s, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  reducedDims = 'TimeLSI',
                                  embedding = 'Time_UMAP', force  = TRUE)
                   
#sudo apt-get install gdal-bin libgdal-dev
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
                   
CD16s <- addTrajectory(CD16s, name = 'StageTrajectory', 
                      trajectory = c('Early Infection', 'Late Infection', 
                                                     'Recovery Phase', 'Uninfected'),
                           groupBy = 'InfectionStages',
                                 reducedDims = 'TimeLSI')   
                   
CD16s <- addTrajectory(CD16s, name = 'KTrajectory', 
                      trajectory = c('1','3','2'),
                           groupBy = 'kClusters',
                                 reducedDims = 'TimeLSI')                   
                   
                   
CD16s <- addImputeWeights(CD16s)    
                   
saveArchRProject(CD16s)
                   
trajCells <- getCellColData(CD16s) %>% as.data.frame() %>% mutate(Cells = rownames(.)) %>%
                   dplyr::select(Cells, Trajectory, InfectionStages) %>% arrange(Trajectory)
ggplot(trajCells, aes(x = Trajectory, y =1, fill = InfectionStages)) + geom_tile()
                   
trajGS <- lapply(c("SlingShot.Curve1", "Trajectory","StageTrajectory","KTrajectory"), function(x) {
    
            getTrajectory(ArchRProj = CD16s, name = x, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
    
    })
trajDev <- lapply(c("SlingShot.Curve1", "Trajectory","StageTrajectory","KTrajectory"), function(x) {
    
            getTrajectory(ArchRProj = CD16s, name = x, useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
    
    })
    
    

pdf('CD16s_Trajectories.pdf')
plotTrajectory(CD16s, name = "SlingShot.Curve1", trajectory = "SlingShot.Curve1", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS[[1]],  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev[[1]],  pal = paletteContinuous(set = "horizonExtra"))
plotTrajectory(CD16s, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS[[2]],  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev[[2]],  pal = paletteContinuous(set = "horizonExtra"))
plotTrajectory(CD16s, name = "StageTrajectory", trajectory = "StageTrajectory", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS[[3]],  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev[[3]],  pal = paletteContinuous(set = "horizonExtra"))
plotTrajectory(CD16s, name = "KTrajectory", trajectory = "KTrajectory", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS[[4]],  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev[[4]],  pal = paletteContinuous(set = "horizonExtra"))
dev.off()
            
pdf('CD16s_Trajectories_Monocle.pdf')
plotTrajectory(CD16s, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]
trajectoryHeatmap(trajGS[[2]],  pal = paletteContinuous(set = "horizonExtra"))
trajectoryHeatmap(trajDev[[2]],  pal = paletteContinuous(set = "horizonExtra"))
dev.off()
                   
                   
### Correlate genescore to motifs
                   
#corMat <- correlateTrajectories(trajGS[[2]],trajDev[[2]], fix1 = 'start', fix2 = 'start', log2Norm2 = FALSE)
#Only correlates between named-matched Genescores and TFs                 
                   
GenScore <- assays(trajGS[[2]])$smoothMat
Dev <- assays(trajDev[[2]])$smoothMat[grepl("z:",rownames(assays(trajDev[[2]])$smoothMat)),] 

pdf('StandardDeviation_Motifs&Genescores.pdf')
hist(apply(GenScore, 1, sd), breaks = 1000)
hist(apply(Dev, 1, sd), breaks = 1000)
dev.off()

GenScoreF <- GenScore[apply(GenScore, 1, sd) > 0.025, ]
DevF <- Dev[apply(Dev, 1, sd) > 50, ]
                
                   
library(CCA)                   
### Let's clean up each matrix
                   
### calculate rank of matrix 
qr(GenScore)$rank
qr(Dev)$rank                   
### return only full rank columns
GenScoreF <- robustbase::fullRank(GenScoreF) 
DevF <- robustbase::fullRank(DevF)
                   
model2 <- CCA::rcc(t(GenScoreF), t(DevF), lambda1 = 0.1, lambda2 = 0.1)
model2$xcoef
CCA::plt.cc(model2, var.label = TRUE)
dev.off()
                   
plt.cc(model2, var.label = FALSE)
plot(x = model2$xcoef[,1], y = model2$ycoef[,1])

model2 <- readRDS('CCA_Object.RDS')
                   
chosenCor <- model2$cor > .90

head(model2$corr.X.xscores)
head(model2$scores$corr.Y.xscores)
##### Let's pull out all the GeneScores that correlate to the Comp1 > 0.9
##### Repeat for Motif Deviations - Correlations > 0.9


findCorrGroup <- function(CCA_model, ComponentNumber, threshold, DevMat,GenMat, plot = FALSE, numCores = 1){

    Comp1Genes <- names(which(abs(CCA_model$scores$corr.X.xscores[,ComponentNumber]) > 0.75))
    Comp1Dev <- names(which(abs(CCA_model$scores$corr.Y.yscores[,ComponentNumber]) > 0.75))

    correlationList <- mclapply(Comp1Genes, function(y){

           tmp <- lapply(Comp1Dev, function(x){
        
                cor(DevMat[which(rownames(DevMat) %in% x),], 
                     GenMat[which(rownames(GenMat) %in% y),], method = 'spearman')
           })
           
           names(tmp) <- Comp1Dev
           
           unlist(tmp)
           
    }, mc.cores = numCores) 
    
    names(correlationList) <- Comp1Genes

    corrTable <- do.call('cbind', correlationList) 

    Comp1Corr <- corrTable %>% 
          as.data.frame() %>%
          rownames_to_column("Dev") %>%
          pivot_longer(-c(Dev), names_to = "Genes", values_to = "cor")
          
    if(plot){
    
       p1 <-  ggplot(Comp1Corr, aes(x=samples, y=Dev, fill=cor)) + 
                  geom_raster() +
                  scale_fill_viridis_c()
      
      return(p1)
    
    }else{
    
        return(Comp1Corr)
        
    }
          
}

allComp <- lapply(1:2, function(x){
                print(x) 
                findCorrGroup(model2, x, 0.9, DevF, GenScoreF, numCores = 40)
                
               })


allComp2 <- lapply(allComp, function(x){

                x %>% mutate(Motifs = gsub('z:','', Dev), GeneScore = gsub('.*:','', Genes)) %>%
                  mutate(Motifs = gsub('_.*','',Motifs))

})

plotList <- lapply(allComp2, function(x){

                  ggplot(x, aes(x=Genes, y=Dev, fill=cor)) + 
                  geom_raster() +
                  scale_fill_viridis_c()
                  
             })

plotList2 <- lapply(allComp2, function(x){

                 x %>% pivot_wider(., id_cols = 'Motifs', 
                         names_from = 'GeneScore', values_from = 'cor') %>% 
                         as.data.table() %>% as.matrix(., rownames = 'Motifs' ) %>%
                         pheatmap::pheatmap()
             })


pdf('CCA_Components.pdf')

plotList
ggplot()
plotList2[[1]]
ggplot()
plotList2[[2]]
  
dev.off()



write.csv(allComp2[[1]], 'CCA_Component1_Correlation.csv')
write.csv(allComp2[[2]], 'CCA_Component2_Correlation.csv')

enrichDataBaseList <- WebGestaltR::listGeneSet()
allGenes <- gsub(".*:",'',rownames(GenScore))
                             
Comp1ORA <- lapply(enrichDataBaseList$name[c(1:9,58:65)], function(x)
                      simplifiedORA(x,
                               allComp2[[1]]$GeneScore, allGenes))
                               
Comp2ORA <- lapply(enrichDataBaseList$name[c(1:9,58:65)], function(x)
                      simplifiedORA(x,
                               allComp2[[2]]$GeneScore, allGenes))
any(lapply(Comp2ORA, length) > 0)
Comp1ORA[[(unlist(lapply(Comp1ORA, length)) > 0]]








## Just females to match modeling?
CD16s <- DietArchRProject(CD16s, keep = c('TileMatrix', "CD16LongCISBPMatrix", "GeneIntegrationMatrix", "GeneScoreMatrix",
                                          "PeakMatrix"), verbose = TRUE) 
CD16f <- subsetArchRProject(CD16s, cells = getCellNames(CD16s)[
                CD16s$predictedGroup_Co2 == 'CD16 Mono' & CD16s$COVID_status == 'Positive'],
                            outputDirectory = 'CD16Infected', dropCells = TRUE, force = TRUE)
saveArchRProject(CD16f)
    
CD16f <- loadArchRProject('CD16Infected')

#CD16f <- addPeakSet(CD16f, TimeRegs, force = TRUE)
#CD16f <- addPeakMatrix(CD16f, force = TRUE)
addArchRThreads(50)
CD16f <- addIterativeLSI(CD16f,  useMatrix = "PeakMatrix", name = "TimeLSI",   iterations = 5, force = TRUE)
CD16f <- addUMAP(CD16f,  reducedDims ="TimeLSI",
  name = "Time_UMAP", force = TRUE)
saveArchRProject(CD16f)
                   

##### Let's try a different type of density plot
EmbedDF <- getEmbedding(CD16s, 'TimeU_UMAP')
   
allDen <- plotDensity(ArchRProj= CD16s, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
    
allDenPASC <- plotDensity(ArchRProj= CD16s, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "PASC_InfectionStages", returnObj = TRUE)
                   
allDen2 <- plotDensity(ArchRProj= CD16f, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
    
allDenPASC2 <- plotDensity(ArchRProj= CD16f, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "PASC_InfectionStages", returnObj = TRUE)

library('ggpubr')
    
plotAll <- ggpubr::ggarrange(plotlist = allDen, nrow = 2, ncol = 2)

plotAllP <- ggpubr::ggarrange(plotlist = allDenPASC , nrow = 4, ncol = 2)
plotAll2 <- ggpubr::ggarrange(plotlist = allDen2, nrow = 2, ncol = 2)

plotAllP2 <- ggpubr::ggarrange(plotlist = allDenPASC2 , nrow = 4, ncol = 2)



    
pdf('CD16_LongitudinalUMAP_Density.pdf')
annotate_figure(plotAll, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))
    

annotate_figure(plotAllP, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))
    
annotate_figure(plotAll2, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))
    

annotate_figure(plotAllP2, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))
    
dev.off()
    
DonorDensity <- lapply(unique(CD16s$PTID), function(x){
    
        plotDensity(ArchRProj= CD16s[CD16s$PTID == x], embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
    
    })
    
DonorDensity2 <- lapply(seq_along(unique(CD16s$PTID)), function(x){
                    ggpubr::ggarrange(plotlist = DonorDensity[[x]], nrow = 2, ncol = 2)
                })
        
pascS <- getCellColData(CD16s) %>% as.data.frame() %>% dplyr::select(PTID, PASC_status) %>% distinct() 
   

                   
pdf('CD16_LongitudinalUMAP_Density_Donor.pdf')
                
                   
for(i in 1:length(DonorDensity2)){
    
    p1 <- annotate_figure(DonorDensity2[[i]], top = text_grob(paste("Donor ",unique(CD16s$PTID)[i], ": ",
                                                                    pascS$PASC_status[pascS$PTID == unique(CD16s$PTID)[i]],sep=""), 
              face = "bold", size = 14))
    print(p1)
    
}
                          
dev.off()
                   
DonorDensity3 <- lapply(unique(CD16s$InfectionStages), function(x){
    
        plotDensity(ArchRProj= CD16s[CD16s$InfectionStages == x], embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "PTID", returnObj = TRUE)
    
    })
                   

pdf('CD16s_LongitudinalUMAP_By_Stage.pdf')
                   
for(i in 1:length(DonorDensity3)){
    
    p1 <- annotate_figure(ggarrange(plotlist = DonorDensity3[[i]]), 
                          top = text_grob(unique(CD16s$InfectionStages)[i], 
              face = "bold", size = 14))
    print(p1)
    
}                          
dev.off() 
                   
############## plot by Days PSO per donor
                   
tmp <- data.frame(COVID = CD16s$COVID_status, PTID = CD16s$PTID, DaysPSO = CD16s$DaysPSO) %>% distinct()
                
COVIDPos <- unique(CD16s$PTID[CD16s$COVID_status == 'Positive'])
DonorDayDensity <- lapply(COVIDPos, function(x){
    
        tmp <- plotDensity(ArchRProj= CD16s[CD16s$PTID == x], embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, plotAxis= FALSE, plotTitle = TRUE, 
                           titleAdjust = c(hjust = 0, vjust= 0), 
                      axisTextSize = 10, addText = 'Day ',
                      identity = "DaysPSO", returnObj = TRUE)
        nEmpty = 5-length(tmp)
         print(nEmpty)
        if(nEmpty  > 0){
         for(i in 1:nEmpty){
            
                tmp <- append(tmp, list(ggplot()))
            
         }F
        }
    
        tmp
    
    })
                  
library(ggpubr)          
library(patchwork)   
pdf('CD16s_LongitudinalUMAP_By_DaysPSO.pdf' )
for(i in 1:length(DonorDayDensity)){
    
    p1 <- annotate_figure(ggarrange(plotlist = DonorDayDensity[[i]]), 
                          top = text_grob(COVIDPos [i], face = "bold", size = 14))
    print(p1)
    
}    
dev.off()
                   

donorPlots <- function(ArchRProj, embeddingName,
                         DonorCol = 'PTID', LongitudinalID = 'days_since_onset',
                         LongitudinalTitleAdd = 'Day', rowNumber = 5){

 metadata <- as.data.frame(getCellColData(ArchRProj)) %>% dplyr::mutate(Cells = rownames(.)) 
  
 embed_df <- getEmbedding(ArchRProj, embedding = embeddingName) %>% as.data.frame() %>% 
    dplyr::mutate(Cells = rownames(.)) %>%
    dplyr::rename(UMAP_1 = grep("_1", colnames(.)), UMAP_2 = grep("_2", colnames(.))) %>%
                dplyr::inner_join(metadata, by = "Cells") 
    
 df_list <- embed_df %>%
         dplyr::group_by(!!sym(DonorCol)) %>% 
        dplyr::mutate(First = ifelse(visit == 1,
                                     TRUE, FALSE)) %>% 
        dplyr::mutate(First = ifelse(visit == 2 & !!sym(DonorCol) == 32245,
                                     TRUE, First)) %>%
                                     
        dplyr::ungroup() %>%
         dplyr::group_split(!!sym(DonorCol), .keep = TRUE) 

    
 endList <- lapply(df_list, function(x){
    
    tmp <- x %>% dplyr::arrange(!!sym(LongitudinalID)) %>% dplyr::group_split(!!sym(LongitudinalID), .keep = TRUE) %>%
             purrr::map( ~ ggplot(., aes(UMAP_1, UMAP_2)) + 
         geom_density_2d_filled(aes(alpha = (..level..))) +  theme_void() + 
        scale_alpha_discrete(range = c(0, 1)) +
        scale_fill_viridis_d(na.value = "transparent", aes(fill = ..density..*10^4)) + 
        theme(legend.position = 'none',  
          plot.margin = margin(t = 0.5, r = 0.5, b =0.5, l = 0.5, unit = "pt"),
         panel.spacing=unit(-1.5, "lines")) +
        ggtitle(ifelse(any(.$First), 
                       paste(unique(unlist(.[,DonorCol])), 
                             LongitudinalTitleAdd, unique(.[,LongitudinalID]), sep =' '),
                       paste(LongitudinalTitleAdd, unique(unlist(.[,LongitudinalID])), sep = ' '))) 
        )
    
     
    if(length(tmp) < rowNumber){
        nEmpty = rowNumber - length(tmp)
        p1 <- ggplot() + theme_void()
         tmp <- append(tmp, as.list(rep(list(p1),nEmpty)))
    }
    tmp
    
    }) 
    
    unlist(endList, recursive = FALSE)

}

                   
finalList <- donorPlots(ArchRProj= CD16s[CD16s$COVID_status == 'Positive'], 
                        embeddingName= 'Time_UMAP',
                         DonorCol = 'PTID', LongitudinalID = 'DaysPSO',
                         LongitudinalTitleAdd = 'Day', rowNumber = 5)

pdf('CD16s_LongitudinalUMAP_By_DaysPSO_All.pdf', width = 7, height = 21 )
   
  ggpubr::ggarrange(plotlist = finalList, ncol = 5, nrow = 18)    
    
dev.off()

                   
pdf('CD16s_LongitudinalUMAP_By_DaysPSO_All.pdf', width = 7, height = 21 )                  
    patchwork::wrap_plots(tmp1, 
                        nrow = 18, ncol = 5)
dev.off()    
                   

                     
CD16s2 <- saveArchRProject(CD16s, outputDirectory = 'CD16s_FullPeakset')
addArchRThreads(50)
CD16s2 <- addPeakSet(CD16s2, StringsToGRanges(results$SampleTileMatrix$tileID), force = TRUE)
CD16s2 <- addPeakMatrix(CD16s2, force = TRUE)
CD16s2 <- addIterativeLSI(CD16s2,  useMatrix = "PeakMatrix", name = "LongLSI", force = TRUE,
                        varFeatures = 237407, iterations = 1)
CD16s2 <- addUMAP(CD16s2,  reducedDims ="LongLSI",
  name = "Long_UMAP", force = TRUE)
saveArchRProject(CD16s2)
                     
LongDen <- plotDensity(ArchRProj= CD16s2, embeddingName= "Long_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)                    
plotLong <- ggpubr::ggarrange(plotlist = LongDen, nrow = 2, ncol = 2)                
                     
pdf('CD16_LongPeakSet_UMAP.pdf')
annotate_figure(plotLong, top = text_grob("CD16s in All Donors & All Peaks", 
              face = "bold", size = 14))
dev.off()
    
################################################################

##Subset donors by end-point

Group1 = c(31207, 31945, 32054, 32140, 32245,32255, 32416, 42409, 32131)  #Low
Group2 = c(31874, 31924, 32038, 32124, 32196, 32209, 32251, 32415) #High

tmp2 <- c(31924, 32196, 32415, 32251, 32209, 32038, 32124, 31874)
tmp1 <- c(32416,  42409,  32255,  31207,  32245,  32209,  32054,  32131,  32140,  31945) 

CD16f1 <- subsetArchRProject(CD16f,
                cells = getCellNames(CD16f)[CD16f$PTID %in% Group1],
                            outputDirectory = 'CD16s_Group1', dropCells = TRUE)
saveArchRProject(CD16f1)

CD16f2 <- subsetArchRProject(CD16f,
                cells = getCellNames(CD16f)[CD16f$PTID %in% Group2],
                            outputDirectory = 'CD16s_Group2', dropCells = TRUE)
saveArchRProject(CD16f2)


MonocleTrajf1 <- getMonocleTrajectories(
                  CD16f1,
                  name = "Trajectory",
                useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  embedding = 'Time_UMAP')
                                  
MonocleTrajf2 <- getMonocleTrajectories(
                  CD16f2,
                  name = "Trajectory",
                useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  embedding = 'Time_UMAP')

CD16f1 <- addMonocleTrajectory(CD16f1, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase'),
                                 monocleCDS = MonocleTrajf1,
                                  groupBy = 'InfectionStages', force = TRUE)
CD16f2 <- addMonocleTrajectory(CD16f2, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase'),
                                 monocleCDS = MonocleTrajf2,
                                  groupBy = 'InfectionStages', force = TRUE)
saveArchRProject(CD16f1)          
saveArchRProject(CD16f2)

CD16G1 <- plotDensity(ArchRProj= CD16f1, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
CD16G2 <- plotDensity(ArchRProj= CD16f2, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)
                                  
pdf('CD16s_Trajectories_Subset_Monocle.pdf')
cowplot::plot_grid(plotlist = CD16G1, ncol = 3)
plotTrajectory(CD16f1, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]
cowplot::plot_grid(plotlist = CD16G2, ncol = 3)
plotTrajectory(CD16f2, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]  
dev.off()



