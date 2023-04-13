

library(ArchR)
library(tidyverse)


setwd('scMACS_Analysis')



#### Let's run the ArchR DAP N-1 tests on CD16s. 
MonoDCE <- loadArchRProject('MonoDC_Edits')

load('template_CD16_summarized_exp.RData')

mySe <- DAPsToSE(results[[3]])

MonoDCE <- addPeakSet(MonoDCE,sort(StringsToGRanges(results[[3]]$Peak)), force = TRUE)
MonoDCE <- addMotifAnnotations(MonoDCE,  motifSet = "cisbp", name = "cisbpMotif", force = TRUE)


# #Now let's downsample by 1 and do it 5 times

metadf <- as.data.frame(getCellColData(MonoDCE))


MonoDCE <- addPeakMatrix(MonoDCE)
saveArchRProject(MonoDCE)

##This is the order of sample list to remove
sampleList <- c("B011-AP0C1W3", "B011-AP0C1W8", "B011-AP0C2W1", 
               "B011-AP0C2W2", "B011-AP0C2W7", "B021-AP0C1W6", 
               "B021-AP0C2W8", "B024-AP0C1W3", "B024-AP0C1W7",
               "B024-AP0C2W7", "B025_FSQCAZ0BP8K-01", "B025_FSQEAZ0BTRD-01",
               "B025_FSQHAZ0BXYJ-01", "B025_FSQJAZ0BH54-01", "B025_FSQJAZ0BKGD-01", 
               "B025_FSQJAZ0BS71-01", "B025_FSQJAZ0BSS6-01", "B033_FSQEAZ0BCPQ-03",
               "B033_FSQHAZ0BH94-03", "B033_FSQHAZ0BVVL-04", "B033_FSQKAZ0BWRX-03",
               "B037_FSQAAZ0BZHX-01", "B037_FSQDAZ0C3CM-01", "B037_FSQFAZ0BFPC-01",
               "B037_FSQFAZ0BNXT-03", "B037_FSQHAZ0C2DL-01", "B038_FSQCAZ0BFDX-03",
               "B038_FSQDAZ0CB41-01", "B038_FSQEAZ0CBL1-01", "B038_FSQFAZ0CBQN-01",
               "B038_FSQKAZ0BG6W-01", "B038_FSQKAZ0C234-01", "FSQAAZ0C1RW-02",
               "FSQBAZ0BSRC-01", "FSQBAZ0C08Q-02", "FSQDAZ0BDS2-03",
               "FSQDAZ0BXW1-02", "FSQEAZ0C2D3-02", "FSQFAZ0BZJQ-02")

CD16List <- lapply(seq_along(sampleList), function(x){
    
    metadf1 <- metadf %>% filter(Sample != sampleList[x]) %>%
            filter(grepl('CD16', EarlyComparison)) 
     subsetCells(MonoDCE, cellNames = rownames(metadf1))
})


CD16DAPs <- lapply(CD16List, function(x){
    
            getMarkerFeatures(x, 
                             groupBy = "EarlyComparison", 
                              useGroups = 'CD16.Mono.Positive',
                              bgdGroups = 'CD16.Mono.Negative',
                              useMatrix = 'PeakMatrix',
                              maxCells = 24744)
    })
CD16DAPsf <- lapply(CD16DAPs, function(x){
    
              getMarkers(x, cutOff = 'FDR < 0.05')[[1]]  
    
    })
orig <- getMarkerFeatures(MonoDCE, 
                             groupBy = "EarlyComparison", 
                              useGroups = 'CD16.Mono.Positive',
                              bgdGroups = 'CD16.Mono.Negative',
                              useMatrix = 'PeakMatrix',
                              maxCells = 24744)
saveRDS(orig, 'CD16_ArchR_DAPs_All.rds')
origM <- getMarkers(orig, cutOff = 'FDR < 0.05')[[1]]  

##
CumulOverlap = list()
for(i in 1:length(CD16DAPsf)){
    
    if(i > 1){
    
        tmp <- plyranges::filter_by_overlaps(
                    makeGRangesFromDataFrame(CD16DAPsf[i], keep.extra.columns = TRUE), 
                                 CumulOverlap[[i-1]])
        
    }else{
        
           tmp <- plyranges::filter_by_overlaps(
                       makeGRangesFromDataFrame(CD16DAPsf[i], keep.extra.columns = TRUE), 
                                  makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE))
        
    }
    
    CumulOverlap = append(CumulOverlap, list(tmp))
}

CumulOverlapCount <- append(dim(origM)[1],lapply(CumulOverlap, length))

sampleSigArch <- lapply(CD16DAPsf, nrow)

CumulNonOverlap = list()
for(i in 1:length(CD16DAPsf)){
    
    if(i > 1){
    
        tmp <- c(plyranges::filter_by_non_overlaps(
                    makeGRangesFromDataFrame(CD16DAPsf[i], keep.extra.columns = TRUE), 
                                 CumulNonOverlap[[i-1]]),
                     CumulNonOverlap[[i-1]])
        
    }else{
        
           tmp <- plyranges::filter_by_non_overlaps(
                       makeGRangesFromDataFrame(CD16DAPsf[i], keep.extra.columns = TRUE), 
                                  makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE))
        
    }
    
    CumulNonOverlap  = append(CumulNonOverlap, list(tmp))
}

CumulNonOverlapCount <- append(0,lapply(CumulNonOverlap, length))

SampleRecall = lapply(CD16DAPsf, function(x){
    
            length(plyranges::filter_by_overlaps(
                makeGRangesFromDataFrame(x, keep.extra.columns = TRUE), 
                        makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE)))/
                length(makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE))
    
    })


SamplePrec = lapply(CD16DAPsf, function(x){
    
            length(plyranges::filter_by_overlaps(
                makeGRangesFromDataFrame(x, keep.extra.columns = TRUE), 
                        makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE)))/
               nrow(x)
    
    })

SampleUnique = lapply(CD16DAPsf, function(x){
    
            length(plyranges::filter_by_non_overlaps(
                makeGRangesFromDataFrame(x, keep.extra.columns = TRUE), 
                        makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE)))
    
    })

SampleConserved = lapply(CD16DAPsf, function(x){
    
            length(plyranges::filter_by_overlaps(
                makeGRangesFromDataFrame(x, keep.extra.columns = TRUE), 
                        makeGRangesFromDataFrame(origM, keep.extra.columns = TRUE)))
    
    })


saveRDS(CD16DAPsf, 'CD16_ArchR_NMinus_AllIterations.RDS')

sumNMinus <- data.frame(Iteration = c(1:(length(sampleList)+1)), 
                SampleRemoved = c('None',sampleList), 
               Recall = unlist(CumulOverlapCount),
               Unique = unlist(CumulNonOverlapCount),
               PeakNumber = c(dim(origM)[1], unlist(sampleSigArch)),
               SampleUnique = c(0, unlist(SampleUnique)),
                SampleConserved = c(dim(origM)[1],unlist(SampleConserved))) %>%
            mutate(UniqueAccum = Unique/dim(origM)[1],
                   Recall = Recall/dim(origM)[1],
                  SampleSpecificPrec = SampleConserved/PeakNumber,
                   SampleSpecificRecall = SampleConserved/dim(origM)[1],
                  PercentUnique = SampleUnique/PeakNumber)

sumNMinus2 <- sumNMinus %>% select(Iteration, SampleRemoved, PercentRecall, PercentUnique) %>%
                 pivot_longer(cols = c(PercentRecall, PercentUnique), names_to = 'Metric', 
                              values_to = 'Percentage')
ggplot(sumNMinus, aes(x = Iteration, y = PercentUnique)) + geom_line()
ggplot(sumNMinus, aes(x = Iteration, y = PercentRecall)) + geom_line()
ggplot(sumNMinus, aes(x = Iteration, y = SampleRecall1)) + geom_point(size =2)+ xlim(2,40) + theme_bw()
ggplot(sumNMinus, aes(x = Iteration, y = SamplePrec1)) + geom_point(size =2 ) + xlim(2,40) + theme_bw()

ggplot(sumNMinus, aes(x = Iteration, y = PercentUnique)) + geom_line()

write.csv(sumNMinus, 'ArchR_DAPs_N_Minus_One_CD16.csv')

pdf('ArchR_NMinus1_Perturbations_CD16.pdf')
ggplot(sumNMinus2, aes(x = Iteration, y = Percentage, group = Metric, color = Metric)) + 
    geom_line() + theme_bw() + ylab('Percentage of Original DAPs') + xlab('N-1 Iteration') + 
    scale_color_manual(labels = c("Conserved Peaks","New Peaks"),
                     values = c("red", "blue")) +  theme(legend.position = c(0.75,0.5)) +
    ggtitle('ArchR DAPs are highly sensitive to N-1 perturbations.')
dev.off()

write.csv(origM, 'CD16_DAPs_ArchR_on_scMACSTiles.csv')

##########################################################################
###### Now for the same for Signac
###########################################################################
##
library(Seurat)
library(Signac)

chromAssay <- CreateChromatinAssay(counts = assays(peakMatf)[[1]], ranges = rowRanges(peakMatf))

seurObj <- CreateSeuratObject(counts = chromAssay, assay = 'ATAC', meta.data= as.data.frame(colData(peakMatf)))

Idents(seurObj) <- as.data.frame(colData(peakMatf))$EarlyComparison

fripScores <-  as.data.frame(getCellColData(MonoDCE, 'FRIP')) %>% mutate(Cells = rownames(.))

newTmp <- as.data.frame(colData(peakMatf)) %>% mutate(Cells = rownames(.)) %>% 
            left_join(.,fripScores, by ='Cells')

saveRDS(seurObj, 'CD16_Mono_SignacObject.RDS')


seurObj@meta.data$FRIP <- newTmp$FRIP

SeurDaps <- FindMarkers(
  object = seurObj,
  ident.1 = "CD16.Mono.Positive",
  ident.2 = "CD16.Mono.Negative",
  min.pct = 0.001,
  logfc.threshold = 0.05,
  test.use = 'LR',
  latent.vars = 'nFrags'
)
##Setting logfc.threshold = 0.025 means a runtime of 36 minutes. 
#### Would take 24 hours and use all the cores to test N-1
##logfc.threshold = 0.05 means a runtime of 8 minutes. 

write.csv(SeurDaps, 'CD16_DAPs_Seurat.csv')

SeurDaps_NMinus <- lapply(seq_along(sampleList), function(x){
   tmp2 <- as.data.frame(colData(peakMatf)) %>% as.data.frame() %>% 
    mutate(newLabel = ifelse(!Sample %in% sampleList[x],EarlyComparison,'Other')) 

   Idents(seurObj) = tmp2$newLabel
    
    FindMarkers(
          object = seurObj,
          ident.1 = "CD16.Mono.Positive",
          ident.2 = "CD16.Mono.Negative",
          min.pct = 0.001,
          logfc.threshold = 0.05,
          test.use = 'LR',
      latent.vars = 'nFrags'
    )


})
saveRDS(SeurDaps_NMinus, 'CD16_Seurat_NMinus.rds')
SeurDaps_NMinus <- readRDS('CD16_Seurat_NMinus.rds')

SeurDapsf <- SeurDaps %>% filter(p_val_adj < 0.05)


#### Let's find features in common
SeurSampleCommon = lapply(SeurDaps_NMinus, function(x){
            Seur_tmp <- x %>%   filter(p_val_adj < 0.05)
            rownames(Seur_tmp)[rownames(Seur_tmp) %in% rownames(SeurDapsf)]
    })

SeurSampleUnique = lapply(SeurDaps_NMinus, function(x){
             Seur_tmp <- x %>% filter(p_val_adj < 0.05)
            rownames(Seur_tmp)[! rownames(Seur_tmp) %in% rownames(SeurDapsf)]
    })


sampleSigSeur <- lapply(SeurDaps_NMinus, function(x){
    
    x %>% filter(p_val_adj < 0.05) %>% nrow()
    
    })

SeurCumulOverlap = SeurDaps_NMinus
for(i in 1:length(SeurDaps_NMinus)){

    if(i == 1){
        SeurCumulOverlap[[i]]  = SeurSampleCommon[[i]]
    }else{
        
       SeurCumulOverlap[[i]] = SeurSampleCommon[[i]][SeurSampleCommon[[i]] %in% 
                                              SeurCumulOverlap[[i-1]]]
    }    

}
SeurCumulOverlapCount <- append(dim(SeurDapsf)[1],lapply(SeurCumulOverlap, length))


SeurCumulNonOverlap = SeurDaps_NMinus
for(i in 1:length(SeurDaps_NMinus)){

    if(i == 1){
        SeurCumulNonOverlap[[i]]  = SeurSampleUnique[[i]][!SeurSampleUnique[[i]] %in% 
                                             rownames(SeurDapsf) ]
    }else{
        
       SeurCumulNonOverlap[[i]] = c(SeurSampleUnique[[i]][!SeurSampleUnique[[i]] %in% 
                                              SeurCumulNonOverlap[[i-1]]],
                                   SeurCumulNonOverlap[[i-1]])
    }    

}
SeurCumulNonOverlapCount <- append(0,lapply(SeurCumulNonOverlap, length))

SeursumNMinus <- data.frame(Iteration = c(1:(length(sampleList)+1)), 
                SampleRemoved = c('None',sampleList), 
               Conserved = unlist(SeurCumulOverlapCount),
               Unique = unlist(SeurCumulNonOverlapCount),
               PeakNumber = c(dim(SeurDapsf)[1], unlist(sampleSigSeur)),
               SampleUnique = c(0, unlist(lapply(SeurSampleUnique, length))),
                SampleConserved = c(dim(SeurDapsf)[1],
                                    unlist(lapply(SeurSampleCommon, length)))) %>%
            mutate(UniqueAccum = Unique/dim(SeurDapsf)[1],
                   Recall = Conserved/dim(SeurDapsf)[1],
                  SampleSpecificPrec = SampleConserved/PeakNumber,
                   SampleSpecificRecall = SampleConserved/dim(SeurDapsf)[1],
                  PercentUnique = SampleUnique/PeakNumber)

pdf('Seur_NMinus1_Perturbations_CD16.pdf')

ggplot(SeursumNMinus, aes(x = Iteration, y = Recall)) + geom_line() + ylim(0,1)
ggplot(SeursumNMinus, aes(x = Iteration, y = UniqueAccum)) + geom_line() 
ggplot(SeursumNMinus, aes(x = Iteration, y = SampleSpecificRecall)) + geom_point(size =2)+ 
    xlim(2,40) + theme_bw() + ylim(0,1)
ggplot(SeursumNMinus, aes(x = Iteration, y = SampleSpecificPrec)) + geom_point(size =2 ) + 
    xlim(2,40) + theme_bw() + ylim(0,1)

dev.off()


# ##########################################################


# ## ArchR Run time tests:

# # We'll test 

peakToTest <- c(204103, 190779, 162531, 110009, 55251, 20966, 4210, 229, 3)

totalPeakSet <- getPeakSet(MonoDCE)
totalPeakSet[sample(c(1:215649), 204103)]

for(i in 1:length(peakToTest)){
    
        MonoDCE <- addFeatureMatrix(
                  input = MonoDCE,
                  features = totalPeakSet[sample(c(1:215649),peakToTest[i])],
                    matrixName = paste('SubMatrix_',i,sep=''))
}


runTimeTest <- lapply(c(1:7), function(x){
    
    print(x)
    startTime = Sys.time()
    
    tmp1 <- getMarkerFeatures(MonoDCE,
                             groupBy = "EarlyComparison", 
                              useGroups = 'CD16.Mono.Positive',
                              bgdGroups = 'CD16.Mono.Negative',
                              useMatrix = paste('SubMatrix_',x,sep=''),
                              maxCells = 24744)
    tmp2 <- getMarkers(tmp1, cutOff = "FDR < 0.05")[[1]]
    runTime <- Sys.time() - startTime 
    tmp2 %>% as.data.frame() %>% 
        mutate(RunTime = runTime, PeakNumber = peakToTest[x])

})

runSummary <- runTimeTest %>% do.call('rbind', .) %>% as.data.frame() %>%
                    group_by(PeakNumber) %>% summarize(DAPNumber = dplyr::n(), RunTime = unique(RunTime)) 


write.csv(runSummary, 'ArchR_RunTime.csv')


# ###########################################################
# ### Functions:


DAPsToSE <- function(daps){
    
    #Generate a GRanges from the character strings found in the peaks column
    tmp1 <- StringsToGRanges(daps$Peak)
    #Generate a data table with the rest of the info.
    tmp2 <- daps[,!'Peak']
    mcols(tmp1) <- tmp2
    tmp3 <- sort(tmp1)
    #Generate the assay list to put into the summarized experiment
    assaysList = lapply(colnames(tmp2), function(x){
                    as.matrix(mcols(tmp3)[,x])
              })
    #Make sure it's named appropriately.
    names(assaysList) = colnames(tmp2)

    rowData = as.data.frame(granges(tmp3))
    rowData$idx = seq(1:length(tmp3))
    #Now load into a summarized experiment format. 
    SummarizedExperiment(
         assays = assaysList,
          #rowRanges = granges(tmp3),
          rowData = rowData,
          metadata = list(Params = list(useMatrix = "PeakMatrix"))
    )
}



getCellPeakMat <- function(ArchR, cellLabel = NULL, column = NULL, 
                           metaColNames  = NULL, numCores =1){
    
    if(any(unlist(lapply(list(cellLabel, column, metaColNames), is.NULL))){
        stop('Error: Inputs missing.')
    }
    
    metadf <- getCellColData(ArchR)
    arrowList <- getArrowFiles(ArchR)
    
    metadf_f <- metadf %>% as.data.frame(.) %>% 
           dplyr::filter(Sample %in% names(arrowList) &
                         grepl(cellLabel,{{ column }}))
    arrowList2 <- arrowList[grepl(paste(unique(metadf_f$Sample),collapse="|"), arrowList)]
       
    peakMatList <- mclapply(seq_along(arrowList2), function(x){
    
        arr1 <- metadf %>% as.data.frame(.) %>% 
           dplyr::filter(Sample %in% names(arrowList2)[x] &
                         grepl(cellLabel,{{ column }}))

        if(dim(arr1)[1] > 0){
            getMatrixFromArrow(arrowList2[x], useMatrix = 'PeakMatrix', 
                                     cellNames = rownames(arr1))
            }else{ NULL}
    }, mc.cores = numCores)
       
    peakMatList <- peakMatList[!unlist(lapply(peakMatList, is.null))]
   
    tmp1 <-peakMatList[[1]]
    colData(tmp1) <- colData(tmp1)[,metaColNames]
    tmp2 <-peakMatList[[2]]
    colData(tmp2) <- colData(tmp2)[,metaColNames]
    peakMat <- cbind(tmp1, tmp2)
    
    for(i in 3:length(peakMatList)){
        tmp <-peakMatList[[i]]
        colData(tmp) <- colData(tmp)[,colnames(colData(tmp)) %in% 
                                     colnames(colData(peakMat))]
        peakMat <- cbind(peakMat, tmp)
    }
    return(peakMat)

}
