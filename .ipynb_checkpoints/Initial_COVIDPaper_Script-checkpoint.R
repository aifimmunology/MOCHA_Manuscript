library(SummarizedExperiment)
library(lmerTest)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)

packrat::on()
STM <- readRDS('SampleTileObject.rds')
library(chromVARmotifs)
data(human_pwms_v2)


library(BiocParallel)
set.seed(2017)

register(MulticoreParam(55, progressbar = TRUE))
STM <- correctGenome(STM, BSgenome.Hsapiens.UCSC.hg38)

STM <- addMotifSet(STM, pwms= human_pwms_v2, w = 7, returnObj = TRUE, motifSetName = "Motifs")

saveRDS(STM, 'SampleTileObject.rds')

cellTypes <- names(assays(STM))
cellTypes <- cellTypes[cellTypes %in% colnames(rowData(STM))]

subSTM <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList= 'Positive')
subSTM <- subsetMOCHAObject(subSTM, subsetBy = 'celltypes', groupList= cellTypes)

rawDataUMAPs <- parallel::mclapply(cellTypes, function(x){
    accMat <- getCellPopMatrix(subSTM, cellPopulation = x)
    accUMAP <- as.data.frame(uwot::umap(t(2^accMat + 1)))
    colnames(accUMAP) <- c('UMAP1', 'UMAP2')
    accUMAP$Sample = rownames(accUMAP)
    accUMAP <- dplyr::left_join(accUMAP, as.data.frame(colData(STM)), by = 'Sample') %>%
                   group_by(PTID) %>%
                        dplyr::mutate(Last =
                                dplyr::case_when(days_since_symptoms == max(days_since_symptoms) ~ 'Last',
                                     days_since_symptoms == min(days_since_symptoms) ~ "First",
                                     TRUE ~ 'Other')) %>%
                dplyr::arrange(days_since_symptoms)

    accUMAP$KMeans <- kmeans(accUMAP[,c('UMAP1', 'UMAP2')], centers =3, nstart = 200)$cluster
    accUMAP
}, mc.cores =30)
saveRDS(rawDataUMAPs, 'UMAPs_AllCellTypes_AllTiles.rds')

##Find optimal Kmeans cluster values
pdf('Identify_OptimalKMeans.pdf')
    lapply(rawDataUMAPs, function(x){

        fviz_nbclust(x[c(1:2)], kmeans, method = "gap_stat")

        })
dev.off()

optimalKMeans <- unlist(lapply(rawDataUMAPs, function(x){
    gapStat <- cluster::clusGap(x[c(1:2)],  kmeans, K.max = 10)
    cluster::maxSE(gapStat$Tab[,'gap'], gapStat$Tab[,'SE.sim'])

}))
names(optimalKMeans) <- CellTypes
## For which cell types does the gap statistic fail?
which(optimalKMeans == 1)
## cDCs, Treg, B Effector, Transitional B cells, and HSPCs


## Re-run KMeans
rawDataUMAPs <- lapply(seq_along(rawDataUMAPs), function(x){
    if(optimalKMeans[x] > 1){
    rawDataUMAPs[[x]]$KMeans <- kmeans(rawDataUMAPs[[x]][,c('UMAP1', 'UMAP2')],
                                       centers = optimalKMeans[x], nstart = 200)$cluster
    }else{rawDataUMAPs[[x]]$KMeans <- 1}
    rawDataUMAPs[[x]]
})


saveRDS(rawDataUMAPs, 'UMAPs_AllCellTypes_AllTiles.rds')

pdf('AllCellTypes_Longitudinal_UMAPs_AllTiles.pdf')

lapply(seq_along(rawDataUMAPs), function(x){
    p1 <- ggplot(rawDataUMAPs[[x]], aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(x = UMAP1, y= UMAP2, shape = Stage, color = Last),
                   size = 2.5) +
               #scale_color_manual(values = c('First' = '#2CB04A', 'Last' = '#E9282A', 'Other' = ) +
        geom_path(aes(x = UMAP1, y = UMAP2, group = PTID), alpha = 0.40,
                 arrow = arrow(angle = 30,
                               length = unit(0.25, "inches"),
          ends = "last", type = "open")) + theme_bw() + ggtitle(cellTypes[x])
    p2 <- ggplot(rawDataUMAPs[[x]], aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(x = UMAP1, y= UMAP2, shape = Stage, color = as.character(KMeans)),
                   size = 2.5) +
        geom_path(aes(x = UMAP1, y = UMAP2, group = PTID), alpha = 0.40,
                 arrow = arrow(angle = 30,
                               length = unit(0.25, "inches"),
          ends = "last", type = "open")) + theme_bw() + ggtitle(cellTypes[x])
    list(p1, p2)
})
dev.off()

### Is it the same group of 10-11 donors whose last time point fall in the same cluster?

largestGroup <- lapply(rawDataUMAPs, function(x){

    domKMeans <- filter(x, Last == 'Last') %>% dplyr::select(Sample, Last, KMeans) %>%
        group_by(KMeans) %>% summarize(Sum1 = dplyr::n()) %>% arrange(Sum1)
    filter(x, KMeans == domKMeans$KMeans[which.max(domKMeans$Sum1)] & Last == 'Last') %>%
        ungroup() %>%
        select(PTID) %>% unlist()

})
names(largestGroup) <- cellTypes

largestGroup[which(optimalKMeans != 1)]
## Ignore cDCs, Treg, B Effector, Transitional B cells, and HSPCs since there was no optimal>
## We can cluster 11 out of 16 cell types. 
## Same group of donors for 9/11 cell types: CD4 Effector, CD14 Mono, CD8 Effector, CD4 Naiv>
##### NK, CD8 Naive, Transitional. 
## Almost identical for:
# CD16 Mono: Missing PTID 32415 
# B memory: Includes additional PTID 32251 - only two KMeans cluster
# NK CD56 Bright: Includes PTID 32220, & 32038


## Core set of eleven PTIDs: c(32416, 42409, 32415, 32255, 31207, 32245, 32209, 32054, 32131>
## This is nearly identical to what we found before - just 32415
## The rest: Group 2 --> c(32220, 32038, 31874, 32124, 31924, 32251, 32196)

############################################################################################>

### Run Variance Decomposition for each cell type longitudinally for all tiles

## Find tiles with at least 80% non-zeros 

rowList <- parallel::mclapply(cellTypes, function(x){
    accMat <- getCellPopMatrix(subSTM, cellPopulation = x)
    colData1 <- colData(subSTM) %>% as.data.frame() %>% filter(Sample %in% colnames(accMat))
    which(rowMeans(accMat[,colData1$Sample[colData1$Stage == 'Early Infection']] != 0) >= 0.>
          rowMeans(accMat[,colData1$Sample[colData1$Stage == 'Late Infection']] != 0) >= 0.8>
          rowMeans(accMat[,colData1$Sample[colData1$Stage == 'Recovery']] != 0) >= 0.8)
}, mc.cores = 35)
names(rowList) <- cellTypes

rm(STM)
gc()



varForm <- paste0(unlist(lapply(c('Sex','Age','days_since_symptoms'),
                                function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
formula1 <- as.formula(paste('exp ~ ',varForm, sep = ''))

## subset down to one cell type and all the rows that passed threshold. 

newsubSTM <- subsetMOCHAObject(subSTM, subsetBy = 'CellTypes', groupList = cellTypes[12])[ro>
newsubSTM$Freq = metadata(newsubSTM)$CellCounts[colnames(newsubSTM), cellTypes[12]]
linRes <- linearModeling(newsubSTM[c(1:50),], formula = formula1, CellType =cellTypes[12],
                         NAtoZero = TRUE, numCores = 1)

varDecomp1 <- calculateVarDecomp(newsubSTM, cellTypes[1], c('Sex','Age','days_since_symptoms'),
                               NAtoZero = TRUE, numCores = 1)
write.csv(varDecomp1, 'VarDecomp_CD4_Effector.csv')
rm(varDecomp1)

saveRDS(linRes, 'ModelForVarDecomp_CD4TEM.rds')

TregRes <-calculateVarDecomp(newsubSTM[c(1:200),], cellTypes[12], c('Sex','Age','days_since_symptoms'),
                               NAtoZero = TRUE, numCores = 50) #11m 35s

TregRes2 <-calculateVarDecomp(newsubSTM[c(1:200),], cellTypes[12], c('Sex','Age','days_since_symptoms'),
                               NAtoZero = TRUE, numCores = 50)

for(y in seq_along(cellTypes)){
    print(cellTypes[y])
    newsubSTM <- subsetMOCHAObject(subSTM, subsetBy = 'CellTypes',
                                   groupList = cellTypes[y])[rowList[[y]],]             
    varDecomp1 <- calculateVarDecomp(newsubSTM, cellTypes[y],
                                      c('Sex','Age','days_since_symptoms'),
                                   NAtoZero = TRUE, numCores = 50)                      
    write.csv(varDecomp1, paste('VarDecomp_', cellTypes[y], '.csv', sep = ''))
    rm(varDecomp1)
    rm(newsubSTM)
    gc()

}

#names(allModelResults) <- cellTypes
#saveRDS(allModelResults, 'VarDecomp_AllCellTypes.rds')

############################################################################################>
### Let's run chromVAR longitudinally for all cell types(within each cell type, over all con>
devList <- runChromVAR(STM, motifSetName = 'Motifs', cellTypeSpecific = TRUE)
saveRDS(devList, 'Dev_AllCellTypes_List.rds')
devList <- readRDS('Dev_AllCellTypes_List.rds')

subsDevlist <- lapply(devList, function(x){
    tmp <- subsetDev(x, subsetBy = 'COVID_status', groupList = 'Positive')
    tmp$time <- tmp$days_since_symptoms - median(tmp$days_since_symptoms)
    tmp$time_sqrd <- tmp$time^2
    tmp
    })
saveRDS(subsDevlist,'Dev_AllCellTypes_SubList.rds')
Group1 <- c(32416, 42409, 32415, 32255, 31207, 32245, 32209, 32054, 32131, 32140, 31945)
Group2 <- c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)

set.seed(2017)
subsDevlist <- readRDS('Dev_AllCellTypes_SubList.rds')

## Run UMAP on TF scores
TFumaps <- lapply(subsDevlist, function(x){
    accUMAP <- as.data.frame(uwot::umap(t(assays(x)[['z']])))
    colnames(accUMAP) <- c('UMAP1', 'UMAP2')
    accUMAP$Sample = rownames(accUMAP)
    accUMAP <- dplyr::left_join(accUMAP, as.data.frame(colData(x)), by = 'Sample') %>%
                   group_by(PTID) %>%
                        dplyr::mutate(Last =
                                dplyr::case_when(days_since_symptoms == max(days_since_symptoms) ~ 'Last',
                                     days_since_symptoms == min(days_since_symptoms) ~ 'First',
                                     TRUE ~ 'Other')) %>%
                dplyr::arrange(days_since_symptoms)
    accUMAP
})
saveRDS(TFumaps, 'UMAPs_AllCellTypes_AllChromVARDeviations.rds')

pdf('AllCellTypes_Longitudinal_UMAPs_chromVARDeviations.pdf')

lapply(seq_along(TFumaps), function(x){
    ggplot(TFumaps[[x]], aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(x = UMAP1, y= UMAP2, shape = Stage, color = Last),
                   size = 2.5) +
               #scale_color_manual(values = c('First' = '#2CB04A', 'Last' = '#E9282A', 'Other>
        geom_path(aes(x = UMAP1, y = UMAP2, group = PTID), alpha = 0.40,
                 arrow = arrow(angle = 30,
                               length = unit(0.25, "inches"),
          ends = "last", type = "open")) + theme_bw() + ggtitle(names(subsDevlist)[x])

})
dev.off()

## Now run modeling

innateCells <- c("CD14 Mono", "NK", "cDC", 'CD16 Mono')

for(y in innateCells){


        Group1Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group1)

        Group1List = modelDeviations(Group1Dev, formula = exp ~ Age + Sex + time + time_sqrd,
                                     numCores = 15)

        saveRDS(Group1List, paste("ChromVAR_Modeling_",y,'Group1.rds',sep = ''))

        rm(Group1Dev)
        rm(Group1List)

        gc()

        Group2Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group2)

        Group2List = modelDeviations(Group2Dev, formula = exp ~ Age + Sex + time + time_sqrd,
                                     numCores = 5)

        saveRDS(Group2List, paste("ChromVAR_Modeling_",y,'Group2.rds',sep = ''))


        rm(Group2Dev)
        rm(Group2List)

        gc()

    }


adaptCells <- c( "CD4 Effector",  "CD8 Effector" ,
                 "CD4 Naive", "B Naive",  "MAIT",
                "NK",  "CD8 Naive",  "NK CD56Bright", "Treg",
                  "B Memory", "B Effector", "Transitional" , "HSPC")

for(y in adaptCells){

for(y in innateCells){


        Group1Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group1)

        Group1List = modelDeviations(Group1Dev, formula = exp ~ Age + Sex + time + time_sqrd>
                                     numCores = 15)

        saveRDS(Group1List, paste("ChromVAR_Modeling_",y,'Group1.rds',sep = ''))

        rm(Group1Dev)
        rm(Group1List)

        gc()

        Group2Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group2)

        Group2List = modelDeviations(Group2Dev, formula = exp ~ Age + Sex + time + time_sqrd,
                                     numCores = 5)

        saveRDS(Group2List, paste("ChromVAR_Modeling_",y,'Group2.rds',sep = ''))


        rm(Group2Dev)
        rm(Group2List)

        gc()

    }


adaptCells <- c( "CD4 Effector",  "CD8 Effector" ,
                 "CD4 Naive", "B Naive",  "MAIT",
                "NK",  "CD8 Naive",  "NK CD56Bright", "Treg",
                  "B Memory", "B Effector", "Transitional" , "HSPC")
        Group1Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group1)

   Group1List = modelDeviations(Group1Dev, formula = exp ~ Age + Sex + time + time_sqrd>
                                     numCores = 35)

        saveRDS(Group1List, paste("ChromVAR_Modeling_",y,'Group1.rds',sep = ''))

        rm(Group1Dev)
        rm(Group1List)

        gc()

        Group2Dev <- subsetDev(subsDevlist[[y]], subsetBy = 'PTID', groupList = Group2)

        Group2List = modelDeviations(Group2Dev, formula = exp ~ Age + Sex + time + time_sqrd,
                                     numCores = 35)

        saveRDS(Group2List, paste("ChromVAR_Modeling_",y,'Group2.rds',sep = ''))


        rm(Group2Dev)
        rm(Group2List)

        gc()

    }

#############################################################################
###### Check whether Macrophage-associated receptors follow similar trends in CD14 monocytes>

tfList <- rownames(readRDS('Dev_AllCellTypes_List.rds')[[1]])
allG1 <- parallel::mclapply(names(assays(subSTM)), function(y){
      modRes <- readRDS(paste("Longitudinal Modeling Results/ChromVAR_Modeling_",y,'Group1.rds', sep = ''))
        names(modRes) <- tf
        modRes
    }, mc.cores = 60)
names(allG1) <- names(assays(subSTM))

gc()

allG2 <- parallel::mclapply(names(assays(subSTM)), function(y){

        modRes <- readRDS(paste("Longitudinal Modeling Results/ChromVAR_Modeling_",y,'Group2.rds'), sep = '')
        names(modRes) <- names(assays(subSTM))
        modRes
    }, mc.cores = 60)
names(allG2) <- names(assays(subSTM))

gc()

All_Coeff <- lapply(names(assays(subSTM)), function(y){

    G1_res <- readRDS(paste("Longitudinal Modeling Results/ChromVAR_Modeling_",
                            y,'Group1.rds', sep=''))
    names(G1_res) <- tfList

    G2_res <- readRDS(paste("Longitudinal Modeling Results/ChromVAR_Modeling_",
                            y,'Group2.rds', sep=''))
    names(G2_res) <- tfList

    tmp <- parallel::mclapply(tfList, function(x){

                G1 <- summary(G1_res[[x]])$coefficients
                G2 <- summary(G2_res[[x]])$coefficients

                #rbind(data.frame(Intercept = rep(G1[1,1],5),PVals = G1[,5], EffectSize = G1[>
                #             Variable = rownames(G1), Group = rep('Group1',5),
                #             TF = rep(x,5), CellType = rep(y, 5)),
                #   # data.frame(Intercept = rep(G2[1,1],5), PVals = G2[,5], EffectSize = G2[,>
                #           Variable = rownames(G2), Group = rep('Group2',5),
                #           TF = rep(x,5), CellType = rep(y, 5)))

  }, mc.cores = 60) %>% do.call('rbind',.)

    gc()

    tmp

    })

names(All_Coeff) <- names(assays(subSTM))

saveRDS(All_Coeff, 'ChromVAR_Modeling_AllCellTypes.rds')

#################################################
###### Let's look specifically at Rexinoids
allRXR <- lapply(All_Coeff, function(x){

        dplyr::filter(x,
                      grepl('RXR',TF) & grepl('time',Variable) & PVals < 0.05)

    })
names(allRXR) <- names(All_Coeff)

devList <- readRDS('Dev_AllCellTypes_List.rds')

RXR_Slopes <- lapply(seq_along(All_Coeff), function(y){
    tmp <- dplyr::filter(All_Coeff[[y]], 
                      grepl('RXR',TF) & grepl('time',Variable) & PVals < 0.05)
    if(dim(All_Coeff[[y]])[1] != 0 & dim(tmp)[1] != 0){
        bothGroups <- dplyr::filter(All_Coeff[[y]],
                      grepl('RXR',TF) & grepl('time',Variable)) %>%
                    tidyr::pivot_wider(id_cols = c('TF','Group', 'Intercept'),
                                       names_from = 'Variable', 
                       values_from = c('PVals', 'EffectSize')) %>%
                    dplyr::group_split(Group)
       lapply(1:2, function(z){
        data.frame(Group = rep(unlist(unique(bothGroups[[z]][,'Group'])),
                               61*dim(bothGroups[[z]])[1]),
                   TF = unlist(lapply(unlist(bothGroups[[z]][,'TF']),
                                      function(x) rep(x,61))),
                   x = rep(seq(-timeList[[y]], -timeList[[y]]+60,1),dim(bothGroups[[z]])[1]),
                   time = rep(seq(0,60,1),dim(bothGroups[[z]])[1]),
                   y = unlist(lapply(seq_along(unlist(bothGroups[[z]][,'Intercept'])),
                        function(x) {
                            unlist(bothGroups[[z]][x,'Intercept']) +
                            unlist(bothGroups[[z]][x,'EffectSize_time'])*
                                seq(-timeList[[y]], -timeList[[y]]+60,1) +
                            unlist(bothGroups[[z]][x,'EffectSize_time_sqrd'])*
                                seq(-timeList[[y]], -timeList[[y]]+60,1)^2
                        }))
                    )
        }) %>% do.call('rbind',.)


    }else { NULL}

})
saveRDS(RXR_Slopes,'RXR_SlopesAcrossCellTypes.rds')

pdf('RXR_TFs_byGroup_AcrossCellTypes.pdf')
lapply(seq_along(RXR_Slopes), function(z){

    if(!is.null(RXR_Slopes[[z]])){
        ggplot(RXR_Slopes[[z]], aes(x = time, y = y, Group= TF, linetype = Group, color = TF\
            ylab('ChromVAR Z-score') + xlab('Days Since Symptoms') + theme_bw() + ggtitle(na>

    }

    })

dev.off()

#p;-\                    PVals_Group1_time < 0.05 & PVals_Group2_time >= 0.05 &
                                    
                        EffectSize_Group1_time > 0 ~ "Group1 Up", 
                     PVals_Group1_time < 0.05 & PVals_Group2_time >= 0.05 &
                        EffectSize_Group1_time < 0 ~ "Group1 Down",
                     PVals_Group1_time >= 0.05 & PVals_Group2_time < 0.05 &
                        EffectSize_Group2_time > 0 ~ "Group2 Up",
                     PVals_Group1_time >= 0.05 & PVals_Group2_time < 0.05 &
                        EffectSize_Group2_time < 0 ~ "Group2 Down", 
                     PVals_Group1_time >= 0.05 & PVals_Group2_time >= 0.05 ~ "No Linear Chan>
                QuadType =
                 dplyr::case_when(PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd < 0>
                        EffectSize_Group1_time_sqrd < 0 & EffectSize_Group2_time_sqrd < 0 ~ >
                     PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd < 0.05 &
                        EffectSize_Group1_time_sqrd > 0 & EffectSize_Group2_time_sqrd > 0 ~ >
                     PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd < 0.05 &
                        EffectSize_Group1_time_sqrd > 0 & EffectSize_Group2_time_sqrd < 0 ~ >
                     PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd < 0.05 &
                        EffectSize_Group1_time_sqrd < 0 & EffectSize_Group2_time_sqrd > 0 ~ >
                     PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd >= 0.05 &
                        EffectSize_Group1_time_sqrd > 0 ~ "Group1 Up",
                     PVals_Group1_time_sqrd < 0.05 & PVals_Group2_time_sqrd >= 0.05 &
                        EffectSize_Group1_time_sqrd < 0 ~ "Group1 Down",
                     PVals_Group1_time_sqrd >= 0.05 & PVals_Group2_time_sqrd < 0.05 &
                        EffectSize_Group2_time_sqrd > 0 ~ "Group2 Up",
                     PVals_Group1_time_sqrd >= 0.05 & PVals_Group2_time_sqrd < 0.05 &
                        EffectSize_Group2_time_sqrd < 0 ~ "Group2 Down",
                     PVals_Group1_time_sqrd >= 0.05 & PVals_Group2_time_sqrd >= 0.05 ~ "No Q>
    write.csv(.,'ChromVAR_Modeling_CD14Mono_forCytoscape.csv')

cDC_Coeff <- parallel::mclapply(rownames(devList[[1]]), function(x){

            G1 <- summary(allG1[['cDC']][[x]])$coefficients
            G2 <- summary(allG2[['cDC']][[x]])$coefficients

     rbind(data.frame(PVals = G1[,5], EffectSize = G1[,1], 
                             Variable = rownames(G1), Group = rep('Group1',5),
                             TF = rep(x,5)),
                data.frame(PVals = G2[,5], EffectSize = G2[,1],
                           Variable = rownames(G2), Group = rep('Group2',5),
                           TF = rep(x,5)))

    }, mc.cores = 55) %>% do.call('rbind', .)

dplyr::filter(cDC_Coeff, TF %in% MNR & grepl('time',Variable) & PVals < 0.05)
## Replicates for NR1H3, RXRG, RXRA, RXRB, NR1H2
dplyr::filter(cDC_Coeff, grepl('SMAD|SP|ESR|ARID|FOXO|BCL|NR2F1|POU2F1',TF) & grepl('time',V>

dplyr::filter(cDC_Coeff, TF %in% c('AHR','ARNT') & 
              grepl('time',Variable) & PVals < 0.05)

########################################################################
###### repeat for promoter accessibility modeling (linear?)

STM <- readRDS('SampleTileObject.rds')

Group1 <- c(32416, 42409, 32415, 32255, 31207, 32245, 32209, 32054, 32131, 32140, 31945)
Group2 <- c(32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)

innateCells <- c("CD14 Mono", "NK", "cDC", 'CD16 Mono')

subSTM <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList = 'Positive')
Group1STM <- subsetMOCHAObject(subSTM[rowRanges(STM)$tileType == 'Promoter',], subsetBy = 'P>
Group2STM <- subsetMOCHAObject(subSTM[rowRanges(STM)$tileType == 'Promoter',], subsetBy = 'P>

for(y in innateCells){
    
    

        Group1List = linearModeling(Group1STM, CellType = y,
                                     formula = exp ~ Age + Sex + days_since_symptoms + (1|PT>
                                    NAtoZero = TRUE,
                                     numCores = 45)
        gc()

        Group1List2 = linearModeling(Group1STM, CellType = y,
                                      formula = exp ~ Age + Sex + days_since_symptoms + (1|P>
                                    NAtoZero = FALSE,
                                     numCores = 45)

        Group1Both <- list(Group1List, Group1List2)
        names(Group1Both) <- c('WithZeros', 'WithoutZeros')
        saveRDS(Group1Both, paste("Promoter_Modeling_",y,'Group1.rds',sep = ''))

        rm(Group1Both)
        rm(Group1List)
        rm(Group1List2)

        gc()

        Group2List = linearModeling(Group2STM,
                                    formula = exp ~ Age + Sex + days_since_symptoms + (1|PTI>
                                    NAtoZero = TRUE,
                                     numCores = 45)
        gc()
         Group2List2 = linearModeling(Group2STM,
                                    formula = exp ~ Age + Sex + days_since_symptoms + (1|PTI>

                                    NAtoZero = FALSE,
                                     numCores = 45)

        Group2Both <- list(Group2List, Group2List2)
        names(Group2Both) <- c('WithZeros', 'WithoutZeros')
               saveRDS(Group2Both, paste("Promoter_Modeling_",y,'Group2.rds',sep = ''))

        rm(Group2Both)
        rm(Group2List)
        rm(Group2List2)

        gc()

    }


adaptCells <- c( "CD4 Effector",  "CD8 Effector" ,
                 "CD4 Naive", "B Naive",  "MAIT",
                "NK",  "CD8 Naive",  "NK CD56Bright", "Treg",
                  "B Memory", "B Effector", "Transitional" , "HSPC")

for(y in adaptCells){

        Group1List = linearModeling(Group1STM, CellType = y,
                                     formula = exp ~ Age + Sex + days_since_symptoms + (1|PT>
                                    NAtoZero = TRUE,
                                     numCores = 45)
        gc()
        Group1List2 = linearModeling(Group1STM, CellType = y,
                                      formula = exp ~ Age + Sex + days_since_symptoms + (1|P>
                                    NAtoZero = FALSE,
                                     numCores = 45)

        Group1Both <- list(Group1List, Group1List2)
        names(Group1Both) <- c('WithZeros', 'WithoutZeros')
        saveRDS(Group1Both, paste("Promoter_Modeling_",y,'Group1.rds',sep = ''))

        rm(Group1Dev)
        rm(Group1List)

        gc()

        Group2STM <- subsetMOCHAObject(cellType1, subsetBy = 'PTID', groupList = Group2)

        Group2List = linearModeling(Group2STM,
                                    formula = exp ~ Age + Sex + days_since_symptoms + (1|PTI>
                                    NAtoZero = TRUE,
                                     numCores = 45)
        gc()
                                    
           Group1List2 = linearModeling(Group1STM,
                                    formula = exp ~ Age + Sex + days_since_symptoms + (1|PTI>

                                    NAtoZero = FALSE,
                                     numCores = 45)

        Group2Both <- list(Group2List, Group2List2)
        names(Group2Both) <- c('WithZeros', 'WithoutZeros')
        saveRDS(Group2Botht, paste("Promoter_Modeling_",y,'Group2.rds',sep = ''))


        rm(Group2Dev)
        rm(Group2List)

        gc()

    }

  ############################################################################################>

######## 
STM <- readRDS('SampleTileObject.rds')

unique(colData(STM)$PTID)

COVIDDonors <- subsetMOCHAObject(STM, subsetBy = 'COVID_status', groupList = 'Positive') %>%>

COVIDDonors %in% unique(Other$patientID)

### Load scRNA
load('CombinedRData.rds')
patientIDs <- unique(Other$patientID)
merge(subset(x = BCells,
       subset = patientID == patientIDs[1]),
      y = c(subset(x = scRNA[[2]],
         subset = patientID == patientIDs[1]),
            subset(x = scRNA[[3]],
         subset = patientID == patientIDs[1])),
      add.cell.ids = c("BCells", "TCells", "MonoDC")

############## Older scRNA analysis
all(MonoDC$pbmc_sample_id %in% sampleList) # TRUE
all(sampleList %in% MonoDC$pbmc_sample_id) # TRUE
## All samples present           

scRNA <- lapply(list('BCell.RDS','TCell.RDS','MonoDC.RDS'), readRDS)
names(scRNA) <- list('BCell.RDS','TCell.RDS','MonoDC.RDS')                                      
                                        
 future::plan("multicore", workers = 60)
options(future.globals.maxSize = 400 * 1024 ^ 3)

#cl <- parallel::makeCluster(3)
#parallel::clusterExport(cl, "Seurat", envir = environment())
scRNA <- lapply(scRNA, function(x){
     print('Normalizing...')
     tmp <- Seurat::NormalizeData(x, verbose = FALSE)
     tmp <- Seurat::FindVariableFeatures(tmp)
     print('Scaling...')
     tmp <- Seurat::ScaleData(tmp, vars.to.regress = "percent.mt")
     print('PCAing...')
     Seurat::RunPCA(tmp, verbose = FALSE)
})

scRNA <- lapply(scRNA, function(x){
     RunUMAP(x, dims = 1:30, return.model = TRUE,verbose = FALSE)
})

scRNA <- lapply(scRNA, function(x){
     tmp <- FindNeighbors(x, dims = 1:30)
     FindClusters(tmp, verbose = FALSE)
})
#stopCluster(cl)

scRNA <- lapply(scRNA, function(x){
     tmp <- FindNeighbors(x, dims = 1:30)
     FindClusters(tmp, verbose = FALSE)
})



lapply(seq_along(scRNA), function(x) { saveRDS(scRNA[[x]], names(scRNA)[[x]])})


patientIDs <- unique(scRNA[[1]]$patientID)
      
      
  merge(subset(x = scRNA[[1]],
       subset = patientID == patientIDs[1]),
      y = c(subset(x = scRNA[[2]],
         subset = patientID == patientIDs[1]),
            subset(x = scRNA[[3]],
         subset = patientID == patientIDs[1])),
      add.cell.ids = c("BCells", "TCells", "MonoDC")



 <- merge(scRNA[[1]],
                   y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBM>

## Plot cell type labels and clusters    
pdf('scRNA_CellType_UMAPs.pdf')
lapply(c('predicted.celltype.l2.5','predicted.Phase','seurat_clusters'), function(y){
    lapply(scRNA, function(x){
     DimPlot(x, group.by = y)
    })
})
dev.off()


## Generate new labels for clusters based on cell types

CellType <- MonoDC$predicted.celltype.l2.5
CellType <- gsub("ASDC|cDC|pDC","cDC",CellType)




## Compress into pseudobulk in SummarizedExperiment format                                    
