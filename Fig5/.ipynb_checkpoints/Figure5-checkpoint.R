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


tileResults <- callOpenTiles(
    MonoDCE,
    cellPopLabel = "predictedGroup_Co2",
    cellPopulations= 'CD16 Mono',
    TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
    Org = org.Hs.eg.db,
    studySignal = studySignal,
    outDir = NULL,
    numCores = 45
)

saveRDS(tileResults, 'CD16_tileResults.rds')

tileResults <- readRDS('CD16_tileResults.rds')

SampleTileObj <-  MOCHA::getSampleTileMatrix( 
    tileResults,
    groupColumn= 'InfectionStages',
    threshold = 0.2,
    numCores = 35,
)
SampleTileObj <- annotateTiles(SampleTileObj)

saveRDS(SampleTileObj, 'CD16_Full_SampleTileMatrix.RDS')

###############################

SampleTileObj <- readRDS('CD16_Full_SampleTileMatrix.RDS')
CD16s <- subsetArchRProject(MonoDCE, 
                            cells = getCellNames(MonoDCE)[
                                MonoDCE$predictedGroup_Co2 == 'CD16 Mono'],
                            outputDirectory = 'CD16sLong', dropCells = TRUE)
saveArchRProject(CD16s)
CD16s <- loadArchRProject('CD16sLong')

CD16s <- addPeakSet(CD16s, rowRanges(SampleTileObj), force = TRUE)
CD16s <- addPeakMatrix(CD16s, force = TRUE)
addArchRThreads(50)
## Run LSI on all tiles called by MOCHA to look for inherent time-related bias
CD16s <- addIterativeLSI(CD16s,  useMatrix = "PeakMatrix", name = "TimeLSI",   iterations = 5, force = TRUE)
CD16s <- addUMAP(CD16s,  reducedDims ="TimeLSI",
  name = "Time_UMAP", force = TRUE)
CD16s <- addClusters(CD16s,  reducedDims ="TimeLSI",
  name = "Clusters3", force = TRUE)
saveArchRProject(CD16s)

metadata1 <- getCellColData(CD16s) %>% as.data.frame() %>%
                dplyr::left_join(as.data.frame(colData(SampleTileObj))[,c('Sample','DaysPSO')], by = 'Sample')


CD16s <- addCellColData(CD16s, data = metadata1$DaysPSO, name = "DaysPSO", 
                        cells = getCellNames(CD16s), force= TRUE)
saveArchRProject(CD16s)

##### Let's try a different type of density plot

   
allDen <- plotDensity(ArchRProj= CD16s, embeddingName= "Time_UMAP", 
                       plotLegend = FALSE, 
                      axisTextSize = 10, 
                      identity = "InfectionStages", returnObj = TRUE)

library('ggpubr')
    
plotAll <- ggpubr::ggarrange(plotlist = allDen, nrow = 2, ncol = 2)
    
pdf('CD16_LongitudinalUMAPs.pdf')

plotEmbedding(CD16s, embedding = "Time_UMAP", name = "DaysPSO")
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "InfectionStages")
plotEmbedding(CD16s, embedding = "Time_UMAP", name = "Clusters3")      

annotate_figure(plotAll, top = text_grob("CD16s in All Donors", 
              face = "bold", size = 14))

dev.off()
    
                   
MonocleTraj <- getMonocleTrajectories(
                  CD16s,
                  name = "Trajectory",
                useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                  principalGroup = 'Early Infection', 
                                  groupBy = 'InfectionStages',
                                  embedding = 'Time_UMAP')
saveRDS(MonocleTraj, 'MonocleTrajectoryObject.rds')
CD16s <- addMonocleTrajectory(CD16s, useGroups = c('Early Infection', 'Late Infection', 
                                                       'Recovery Phase', 'Uninfected'),
                                 monocleCDS = MonocleTraj,
                                  groupBy = 'InfectionStages', force =TRUE)
                   
                  
saveArchRProject(CD16s)

CD16s <- addMotifAnnotations(CD16s, name = "cisbp", force = TRUE)
CD16s <- addDeviationsMatrix(CD16s, peakAnnotation = 'cisbp', matrixName = 'CD16LongCISPMatrix')
                   
trajCells <- getCellColData(CD16s) %>% as.data.frame() %>% mutate(Cells = rownames(.)) %>%
                   dplyr::select(Cells, Trajectory, InfectionStages) %>% arrange(Trajectory)
ggplot(trajCells, aes(x = Trajectory, y =1, fill = InfectionStages)) + geom_tile()
                   
trajGS <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory', 
                        useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajDev <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
trajPeak <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "PeakMatrix", log2Norm = FALSE)    
    

pdf('CD16s_Fig5_Trajectories.pdf')

plotTrajectory(CD16s, name = "Trajectory", trajectory = "Trajectory", 
                embedding = 'Time_UMAP')[[1]]

plotTrajectoryHeatmap(trajGS,  pal = paletteContinuous(set = names(ArchRPalettes)[19]))

plotTrajectoryHeatmap(trajDev,  pal = paletteContinuous(set = names(ArchRPalettes)[17]))

trajectoryHeatmap(trajPeak,  pal = paletteContinuous(set = "horizonExtra"))

dev.off()

            
## generate a bottom annotation for the heatmap , that marks DayPSO across the pseudotime trajectory
trajectoryAnnotation = data.frame(DaysPSO = CD16s$days_since_symptoms, Trajectory = CD16s$Trajectory2) %>%
                   dplyr::filter(!is.na(DaysPSO)) %>% dplyr::mutate(y = 1)
                   
pdf('CD16s_Fig5_DaysPSO_Over_Trajectory.pdf')
ggplot(trajectoryAnnotation, aes(x = Trajectory, y = y, fill = DaysPSO)) + geom_tile() + theme_minimal()
dev.off()
       

#### Now look for pathway enrichment along pseudotime for GeneScore and PeakMatrix
CD16s <- loadArchRProject('CD16sLong')

allGenes <- getGenes(CD16s)

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

#Get the trajectory object
trajGS <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory', 
                        useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajDev <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
trajPeak <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "PeakMatrix", log2Norm = FALSE)    
    
## Extract the matrix, so that we can pull the genes. 
trajGSMat <- plotTrajectoryHeatmap(trajGS, returnMatrix = TRUE)
trajPeakMat <- plotTrajectoryHeatmap(trajPeak, returnMatrix = TRUE)
                                
write.csv(trajGSMat, 'Pseudotime_GeneScoreMatrix.csv')
write.csv(trajPeakMat, 'Pseudotime_PeakMatrix.csv')

trajTiles <- annotateTiles(MOCHA::StringsToGRanges(gsub("_","-",rownames(trajPeakMat))),
                          TxDb =TxDb.Hsapiens.UCSC.hg38.refGene, Org = org.Hs.eg.db)
write.csv(trajTiles, 'Pseudotime_MOCHA_Peaks.csv')
trajPromoterTiles <- unique(unlist(
                                lapply(plyranges::filter(trajTiles, tileType == 'Promoter')$Gene,
                                       function(x){strsplit(x, split= ', ')})))      

GE_Trajectory <- simplifiedORA(enrichDataBaseList$name[9],
                               unique(gsub('.*:','', rownames(trajGSMat))),  
                        unique(as.character(allGenes$symbol)))
GE_TileTrajectory <- simplifiedORA(enrichDataBaseList$name[9],
                               unique(trajPromoterTiles),  
                        unique(as.character(allGenes$symbol)))   

trajPathways <- rbind(data.frame(GE_Trajectory, Trajectory = 'GeneScore'), 
      data.frame(GE_TileTrajectory, Trajectory = 'PromoterTiles'))

write.csv(trajPathways, 'Fig5_PseudotimePathways.csv')
                
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                   
                   
pathwayAnnot_Traj = unlist(lapply(trajPathways$description, 
                             function(x) findTrunk(x, Hierarch, ReactID)))
### Manually correct the few unannotate pathways
pathwayAnnot_Traj[grepl('mutants|Diseases|Influenza|Signaling by BRAF and RAF fusions', pathwayAnnot_Traj)] = 'Disease'
pathwayAnnot_Traj[grepl('external stimuli|ATF4 activates genes|HSP90 chaperone cycle for steroid hormone receptors|Regulation of Hypoxia-inducible Factor', 
                        pathwayAnnot_Traj)] = 'Cellular responses to stimuli'      
pathwayAnnot_Traj[grepl('Mitotic G1-G1/S phases|G1|Nuclear Envelope|Ubiquitin-dependent degradation of Cyclin D1', 
                        pathwayAnnot_Traj)] = 'Cell Cycle'                           
pathwayAnnot_Traj[grepl('Clathrin derived vesicle budding', pathwayAnnot_Traj)] = 'Developmental Biology'         
pathwayAnnot_Traj[grepl('Intrinsic Pathway for Apoptosis', pathwayAnnot_Traj)] = 'Programmed Cell Death'       
pathwayAnnot_Traj[grepl('CDT1 association with the CDC6:ORC:origin complex|Cleavage of Growing Transcript', pathwayAnnot_Traj)] = 'DNA Replication'  
pathwayAnnot_Traj[grepl('Mitophagy', pathwayAnnot_Traj)] = 'Autophagy'
pathwayAnnot_Traj[grepl('Peroxisome proliferator-activated receptor alpha', 
                        pathwayAnnot_Traj)] = 'Metabolism'
pathwayAnnot_Traj[grepl('Rho GTPase cycle|Signaling by TGF-beta family members|mTOR signalling', pathwayAnnot_Traj)] = 'Signal Transduction'   
pathwayAnnot_Traj[grepl('ROS, RNS production in phagocytes|Innate Immune', pathwayAnnot_Traj)] = 'Immune System'   
trajPathways$TopLevel <- pathwayAnnot_Traj   
                                  
write.csv(trajPathways, 'Fig5_PseudotimePathways.csv')
                                  
trajPathways <- dplyr::group_by(trajPathways, TopLevel, Trajectory) %>% 
                                  dplyr::mutate(PathwayN = dplyr::n())
                                  
library(ggVennDiagram)
                                  
immuneSubset <- trajPathways %>% dplyr::filter(TopLevel == 'Immune System') 
immuneSubset$SecondLevel = immuneSubset$description
immuneSubset$SecondLevel[grepl(
    'Toll|FCGR|Dectin|NLR|TLR|MyD88|FCERI|Neutrophil|phagocyte|TAK1|IKK|Advanced glycosylation endproduct receptor|CLRs|DC-SIGN|pathogen-associated DNA|mediated induction of type I IFNs|phagocytic|MAP kinase', 
    immuneSubset$SecondLevel)] = 'Innate Immune System'
immuneSubset$SecondLevel[grepl('IFNG signaling|Interleukin|IL|Interferon|AP-1|IFN-stimulated|ISG15 antiviral|signaling by CBL|noncanonical NF-kB signaling', immuneSubset$SecondLevel)] = 'Cytokine Signaling in the Immune System'
immuneSubset$SecondLevel[grepl('NFAT|B cells|BCR|Adaptive|Antigen|CTLA4|CD28|TCR|ER-Phagosome|Immunoregulatory interactions|MHC|antigen', immuneSubset$SecondLevel)] = 'Adaptive Immune System'         
                                  
TrajpieChartSummary <- dplyr::filter(immuneSubset,Trajectory == 'PromoterTiles') %>%
                             group_by(SecondLevel) %>% summarize(PathwayType = dplyr::n()) %>%
                             mutate(prop = PathwayType / sum(.$PathwayType)) %>%
                             mutate(ypos = cumsum(prop) - 0.5*prop)               TrajpieChartSummary2 <- dplyr::filter(immuneSubset,Trajectory == 'GeneScore') %>%
                             group_by(SecondLevel) %>% summarize(PathwayType = dplyr::n()) %>%
                             mutate(prop = PathwayType / sum(.$PathwayType)) %>%
                             mutate(ypos = cumsum(prop) - 0.5*prop)                    

                                  
## now plot them. 
                                  
pdf('Fig5_TrajectoryPathwayHits.pdf')
         
ggplot(trajPathways, 
       aes(x = reorder(TopLevel, - PathwayN), fill = TopLevel)) + 
                                  geom_bar(stat = "count") + 
                             theme_bw() + 
    theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0.25)) +
    facet_wrap(~Trajectory, scales = "free_x") + theme(legend.position = "none")
                             
ggplot(TrajpieChartSummary, aes(x= '', y = PathwayType, fill = SecondLevel)) +
      geom_bar(stat="identity", width=1) +
      coord_polar('y', start = 0) +
      geom_text(aes(y = ypos, label = SecondLevel), color = "white", size=6) +
      theme_void()
                           
ggplot(TrajpieChartSummary2, aes(x= '', y = PathwayType, fill = SecondLevel)) +
      geom_bar(stat="identity", width=1) +
      coord_polar('y', start = 0) +
      geom_text(aes(y = ypos, label = SecondLevel), color = "white", size=6) +
      theme_void()
                                                  
dev.off()

                

################################################################################


#### Now for pseudobulk analysis

#################################################################################

### Linear Modeling of Promoter Accessibility
                   
                   
##Insert Samir's code here
                   
pathwayHits <- read.csv('covid_pathway_hits.csv')
                
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                   
                   
pathwayAnnot = unlist(lapply(pathwayHits$description, 
                             function(x) findTrunk(x, Hierarch, ReactID)))
### Manually correct the few unannotate pathways

pathwayAnnot[grepl('TGF-beta family members', pathwayAnnot)] = 'Signal Transduction' 
pathwayAnnot[grepl('Innate Immune', pathwayAnnot)] = 'Immune System' 
pathwayAnnot[grepl('Diseases', pathwayAnnot)] = 'Disease'                             
pathwayHits$TopLevel <- pathwayAnnot               
                     
Level2 <- pathwayAnnot
Level2[grepl('TAK1|ZBP1|IRAK1|NOD1|TLR|MyD88|TRAF6|Toll-like', pathwayHits$description) & 
       grepl('Immune System', pathwayAnnot)] = 'Innate Immune System'
Level2[grepl('Interleukin|MAP kinase', pathwayHits$description) & 
       grepl('Immune System', pathwayAnnot)] = 'Signaling by Interleukins'
pathwayHits$Level2 <-Level2
write.csv(pathwayHits, 'covid_pathway_hits_annotated.csv')
                       
pieChartSummary <- dplyr::filter(pathwayHits, TopLevel == 'Immune System') %>%
                             group_by(Level2) %>% summarize(PathwayType = dplyr::n()) %>%
                             mutate(prop = PathwayType / sum(.$PathwayType)) %>%
                             mutate(ypos = cumsum(prop) - 0.5*prop)
library(treemap)
                             
pathwayHits <- dplyr::group_by(pathwayHits, TopLevel) %>% mutate(count_name = dplyr::n())
                             
pdf('Fig5_PathwayHits.pdf')
         
ggplot(pathwayHits, aes(x = reorder(TopLevel, -count_name), fill = TopLevel)) + geom_bar(stat = 'count') + 
                             theme_bw() + theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0.25))
                             
ggplot(pieChartSummary, aes(x= '', y = PathwayType, fill = Level2)) +
      geom_bar(stat="identity", width=1) +
      coord_polar('y', start = 0) +
      geom_text(aes(y = ypos, label = Level2), color = "white", size=6) +
      theme_minimal()
                             
treemap(pieChartSummary, index = "Level2", vSize = "PathwayType", type = "index", palette = "Dark2",
         inflate.labels = T)
                             
dev.off()

                                
#########################################################
### lets do upset plots to compare
                             
pathwayHits <- read.csv('covid_pathway_hits_annotated.csv')

                                
library(UpSetR)
                        
upsetInput <- list('MOCHA: Zero-inflated Modeling' = pathwayHits$description,
                'Monocle3: GeneScore Pseudotime' = dplyr::filter(trajPathways, 
                                                Trajectory == 'GeneScore')$description,
                'Monocle3: Promoter Tile Pseudotime' = dplyr::filter(trajPathways, 
                                                Trajectory == 'PromoterTiles')$description )
             
##Calculate the percent overlap with any other method
    
pathwayRep <- data.frame(Approach = names(upsetInput),
                         Reproducibility = c(
                         sum(upsetInput[[1]] %in% unlist(upsetInput[-1]))/length(upsetInput[[1]]),
                         sum(upsetInput[[2]] %in% unlist(upsetInput[-2]))/length(upsetInput[[2]]),
                         sum(upsetInput[[3]] %in% unlist(upsetInput[-3]))/length(upsetInput[[3]])))                   

                             
pdf('Figure5_TrajectoryComparison.pdf')
                                
UpSetR::upset(fromList(upsetInput), order.by = "freq", sets.bar.color=c("blue","red", "green"),
            main.bar.color = c("blue","purple","red", "green", 'black', "yellow3", "cyan") )

ggplot(pathwayRep, aes(y =  Reproducibility, x = Approach, fill = Approach)) + geom_col()    +
            theme_bw() + ylim(0,1) 
dev.off()                  
                                
                                
#### Let's do a quick comparison to test whether the overlap is more correlated with true time than the trajectory unique. 
                             
peakMatrix <- getMatrixFromProject(CD16s, 'PeakMatrix')
peakMat <- assays(peakMatrix)[[1]]
rownames(peakMat) <- MOCHA::GRangesToString(getPeakSet(CD16s))
                  
modelingGenes <- unlist(stringr::str_split(paste(pathwayHits$userId, sep = ";"), ";")) %>% unique() 
promoGenes <- unlist(stringr::str_split(paste(dplyr::filter(trajPathways, 
                        Trajectory == 'PromoterTiles')$userId, sep = ";"), ";")) %>% unique() 
overlapGenes <- promoGenes[promoGenes %in% modelingGenes]
                             
trajTilesOverlap <- trajTiles[grepl(paste(overlapGenes, collapse = "|"), trajTiles$Gene)]
                             

NonoverlapGenes <- promoGenes[!promoGenes %in% modelingGenes]
                      
trajTilesNonOverlap <- trajTiles[grepl(paste(NonoverlapGenes, collapse = "|"), trajTiles$Gene) &
                                !grepl(paste(overlapGenes, collapse = "|"), trajTiles$Gene)]
trajTilesNonOverlap <- trajTiles[(grepl(paste(NonoverlapGenes[c(1:1600)], collapse="|"), trajTiles$Gene) |
                                  grepl(paste(NonoverlapGenes[c(1601:3197)], collapse="|"), trajTiles$Gene)) &
                                 !grepl(overlapGenes, trajTiles$Gene]
                                        
nonoverlapMat <- peakMat[rownames(peakMat) %in% MOCHA::GRangesToString(trajTilesNonOverlap),]
overlapMat <- peakMat[rownames(peakMat) %in% MOCHA::GRangesToString(trajTilesOverlap),]


## Test relationship to time
timeDf <- data.frame(Cells = colnames(nonoverlapMat),
                     'Time' = getCellColData(CD16s, 'DaysPSO')[colnames(nonoverlapMat),1], 
                     'Trajectory' = getCellColData(CD16s, 'Trajectory')[colnames(nonoverlapMat),1],
                    'PTID' = getCellColData(CD16s, 'PTID')[colnames(nonoverlapMat),1],
                    'Stages' = getCellColData(CD16s, 'InfectionStages')[colnames(nonoverlapMat),1])
all(timeDf$Cells == colnames(nonoverlapMat)) #TRUE
                             
timeDf$Time[timeDf$Time < 0] = 200
# Test for percent variance explained by time.
cl <- parallel::makeCluster(35)
formula1 = as.formula(exp ~ c(1|Stages) + c(1|Trajectory)+c(1|PTID))
parallel::clusterExport(cl=cl, varlist=c("formula1","nonoverlapMat", "overlapMat", "timeDf",'extractVarDecomp'), envir=environment())
                             
CorrVal <- pbapply::pblapply(c(1:dim(nonoverlapMat)[1]),
        function(x) {
           #print(x)
          data.frame(Time = MOCHA:::weightedZISpearman(nonoverlapMat[x,], timeDf$Time), 
                     Trajectory = MOCHA:::weightedZISpearman(nonoverlapMat[x,], timeDf$Trajectory))

}, cl = cl) %>% do.call('rbind',.)
                             
sum(abs(CorrVal$Time) > abs(CorrVal$Trajectory))
sum(abs(CorrVal$Time) < abs(CorrVal$Trajectory))                          
CorrVal2 <- pbapply::pblapply(c(1:dim(overlapMat)[1]),
        function(x) {
           #print(x)
          data.frame(Time = MOCHA:::weightedZISpearman(overlapMat[x,], timeDf$Time), 
                     Trajectory = MOCHA:::weightedZISpearman(overlapMat[x,], timeDf$Trajectory))

}, cl = cl) %>% do.call('rbind',.)
                             
sum(abs(CorrVal2$Time) > abs(CorrVal2$Trajectory))
sum(abs(CorrVal2$Time) < abs(CorrVal2$Trajectory))
                             
data.frame(Corr = abs(c(CorrVal$Time, CorrVal2$Time, CorrVal$Trajectory, CorrVal2$Trajectory)), 
           Type = c(rep('Time', length(CorrVal$Time)+length(CorrVal2$Time)),
                    rep('Trajectory', length(CorrVal$Time)+length(CorrVal2$Time))),
            Group = c(rep('TrajectoryOnly', length(CorrVal$Time)),
                        rep('Common', length(CorrVal2$Time)),
                    rep('TrajectoryOnly', length(CorrVal$Time)),
                        rep('Common', length(CorrVal2$Time)))) %>%                                                                    ggplot(., aes(x = Corr, fill = Group)) + 
                      geom_histogram(alpha = 0.3) + theme_bw() +
                      facet_wrap(~ Type)
                             
data.frame(Corr = abs(c(CorrVal$Time, CorrVal2$Time, CorrVal$Trajectory, CorrVal2$Trajectory)), 
           Type = c(rep('Time', length(CorrVal$Time)+length(CorrVal2$Time)),
                    rep('Trajectory', length(CorrVal$Time)+length(CorrVal2$Time))),
            Group = c(rep('TrajectoryOnly', length(CorrVal$Time)),
                        rep('Common', length(CorrVal2$Time)),
                    rep('TrajectoryOnly', length(CorrVal$Time)),
                        rep('Common', length(CorrVal2$Time)))) %>%                                                                    ggplot(., aes(x = Corr, fill = Group)) + 
                      geom_density(alpha = 0.3) + theme_bw() +
                      facet_wrap(~ Type, scales = 'free')     
                                        
cor(timeDf$Time[timeDf$Time != 200], timeDf$Trajectory[timeDf$Time != 200], method = 'spearman')
    
exp ~ c(1 | Stages) + c(1 | Trajectory) + c(1 | PTID)
lmerTest::lmer(formula = as.formula(exp ~ c(1|Stages) + c(1 | PTID)), data = df1)
                             
suppressMessages(varDecompRes <- pbapply::pblapply(c(1:dim(nonoverlapMat)[1]),
        function(x) {
           
           tryCatch({
                df1 <-  data.frame(exp = as.numeric(nonoverlapMat[x,]), 
                timeDf, stringsAsFactors = TRUE)
                lmem1 <- lmerTest::lmer(formula = formula1, data = df1)
                tmp_df <- extractVarDecomp(lmem1, c('Time', 'PTID', 'Trajectory'))
                rm(lmem1, df1)
                tmp_df
           }, error = function(e){
                tmp_df <- rep(NA, (2 + 1)) 
                names(tmp_df) <- c('Time','Trajectory', 'Residual')
                tmp_df
           })

}, cl = cl), classes = "message")
varDecomp_df <- do.call('rbind',varDecompRes)


extractVarDecomp <- function(lmem1, variableList){

    lmem_re <- as.data.frame(lme4::VarCorr(lmem1))
    row.names(lmem_re) <- c(lmem_re$grp)

    lmem_re <- lmem_re[c(variableList,'Residual'), ]
    fix_effect <- lme4::fixef(lmem1)  #get fixed effect
    lmem_re$CV <- lmem_re$sdcor/fix_effect

    normVar <- (lmem_re$vcov)/sum(lmem_re$vcov)
    names(normVar) <- row.names(lmem_re) 

    return(normVar)

}
                    
#################################################################
                   
#### 
library(chromVAR)
library(chromVARmotifs)
data(human_pwms_v2)
                                  
STM_All <- readRDS('../COVID scATAC Manuscript/SampleTileObject.rds')

STM_All <- addMotifSet(STM_All, pwms = human_pwms_v2,
                      w = 7, returnObj = TRUE, motifSetName = "Motifs")
                   
saveRDS(STM_All, 'MOCHA_All_SampleTileMatrix.RDS')
STM_All <- readRDS('MOCHA_All_SampleTileMatrix.RDS')
                                  
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(tidyverse)
library(plyranges)
              
object1 <- MOCHA::combineSampleTileMatrix(STM)


                                                       
##### Compare with Monocle3 ChromVar trajectory
                                                       
devModels <- read.csv('FullCovid_TF_results.csv')
    
CD16s <- loadArchRProject('CD16sLong')
trajDev <- getTrajectory(ArchRProj = CD16s, name = 'Trajectory',
                         useMatrix = "CD16LongCISBPMatrix", log2Norm = FALSE)
                             
trajDev_res <- plotTrajectoryHeatmap(trajDev, returnMatrix = TRUE)                   

write.csv(trajDev_res, 'Pseudotime_ChromVARMatrix.csv')                      
                                
library(UpSetR)
                        
upsetInput2 <- list('MOCHA: ChromVAR Modeling' = devModels$TF[devModels$pval < 0.05],
                'Monocle3: ChromVAR Pseudotime' = gsub("z:|_.*","",rownames(trajDev_res)))
             
##Calculate the percent overlap with any other method
    
chromVarRep <- data.frame(Approach = names(upsetInput2),
                         Reproducibility = c(
                         sum(upsetInput2[[1]] %in% unlist(upsetInput2[-1]))/length(upsetInput2[[1]]),
                         sum(upsetInput2[[2]] %in% unlist(upsetInput2[-2]))/length(upsetInput2[[2]])))
                             
library(UpSetR)
                             
library(venneuler)

vennList <- c(sum(fromList(upsetInput2)[1] & !fromList(upsetInput2)[2]),
              sum(!fromList(upsetInput2)[1] & fromList(upsetInput2)[2]),
              sum(fromList(upsetInput2)[1] & fromList(upsetInput2)[2]))
                             
names(vennList) <- c("M", "P", "M&P")
                             
pdf('Figure5_ChromVAR_TrajectoryComparison.pdf')
                                
UpSetR::upset(fromList(upsetInput2), order.by = "freq", sets.bar.color=c("blue","red"),
            main.bar.color = c("blue","purple","red") )
                             
ggplot(chromVarRep, aes(y =  Reproducibility, x = Approach, fill = Approach)) + geom_col()    +
            theme_bw() + ylim(0,1) 
                             
plot(venneuler(vennList))
                             
                          
dev.off()                  
               

###### Functions used. 
                        
simplifiedORA <- function(database, foreground, background){
    
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                            referenceGeneType= "genesymbol")
    
    }

              
findTrunk <- function(pathway, TreeStructure, IDMatrix, exportTree = FALSE){

    specID <- IDMatrix[match(pathway,IDMatrix[,2]),1]

    if(!(any(TreeStructure[,2] %in% specID))){

        return(pathway)

    }

    treeID <- TreeStructure[which(TreeStructure[,2] %in% specID)[1],]
    descriptionTree = pathway

    while(any(TreeStructure$V2 %in% treeID$V1)){

        treeID <- TreeStructure[which(TreeStructure$V2 %in% treeID$V1)[1],]
        descriptionTree <- paste(descriptionTree,  
                                 IDMatrix[match(treeID[,1],IDMatrix[,1])[1],2], sep = ", ")

    }

    TrunkDescription <- IDMatrix[match(treeID[,1],IDMatrix[,1]),2]


    if(exportTree){
    
        return(descriptionTree)
        
    }
    return(TrunkDescription)

}


###
                        
pullLmerCoef <- function(LmerList, Coeff_Num, numCores = 1){
    
    parallel::mclapply(seq_along(LmerList), function(x) {
        summary(LmerList[[x]])$coefficients[,Coeff_Num]
        }, mc.cores = numCores )
    
}