############################################################
#      Simulating Single Cell ATAC-seq Data
############################################################

library(MOCHA)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

cellNumber = c(10,25, 50, 75, 100, 150, 200, 250, 500, 1000, 1500, 2000, 3000, 5000)
iteration1 = c(1:10)

grid <- expand.grid(cellNumber, iteration1)
## didn't work - parallelization failed

cl = parallel::makeCluster(15)

parallel::clusterEvalQ(cl, {
        library('BSgenome.Hsapiens.UCSC.hg38')
    library('GenomicRanges')
      })

## Get all peaksets
locationList <- mclapply(1:10, function(XX){
    getSimulatedPeakSet(peakNumber = 200000, peakSize = 500, peakPeriodicity = 10, 
                                Genome = BSgenome.Hsapiens.UCSC.hg38, allLocations = TRUE)
}, mc.cores = 10)

names(locationList) <- c(1:10)
saveRDS(locationList, 'TruePeakSet.rds')

plyranges::write_bed(locationList, 'GroundTruth_SimulatedPeakset.bed')

#duplicate to match grid
locationList2 <- locationList[grid[,2]]
#Make into a list and join with truePeakList
iterList1 <- split(grid,seq(nrow(grid)))
iterList2 <- mapply(append, iterList1, locationList2 , SIMPLIFY=FALSE)

pbapply::pblapply(cl = cl, X = iterList2, simulateAndWriteFragments)

simulateAndWriteFragments <- function(YY, force = FALSE){
    cellNum2 = as.numeric(YY[[1]])
    y = YY[[2]]
    locPeaks <- YY[[3]]

    fileName = paste('SimulatedFragments2/SimulatedData_CellNumber',  cellNum2, 'Round',y, sep = '_')
    #existingFiles <- list.files(pattern = paste("^",gsub(".*/","", fileName), ".csv", sep =''), recursive = TRUE)
    fragTmp <- MOCHA::simulateFragments(nCells = cellNum2, meanFragsPerCell = 4000, fragThreshold = 1000, 
                                            allLocationsGR = locPeaks, FRIP = 0.95, meanLengths = c(75, 200), 
                                            lengthProbability = c(0.9, 0.1))
            
        
    write.csv(as.data.frame(fragTmp), paste(fileName,'.csv', sep = ''))
    rm(fragTmp)
  
}

stopCluster(cl)
gc()

############################################################
#      Coverage Tracks and Peak Calling
############################################################

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db


blackList <- makeGRangesFromDataFrame(data.frame(seqnames = 'chr1', start = 20, end = 250, strand = '+'))


source("comparisonsHelperFunctions_v2.R")
library(ArchR)
library(MOCHA)
library(ggplot2)
library(data.table)
library(pbapply)
library(RPushbullet)
library(scales)
library(parallel)

studySignal <-  4000
study_prefactor <- 3668 / studySignal


csvDir <- "/home/jupyter/scMACS_Analysis/simulatedData/SimulatedFragments2"

globalNCores <- 1

# Generate coverage bedgraphs from simulated data up front
# for MACS2 and HOMER

covFileBedGraphList <- mclapply(mc.cores =10, X = split(grid,seq(nrow(grid))), function(XX){
  cellNumber <- as.numeric(XX[1])
  frip <- XX[2]

  popFrags <- getSimulatedFrags(cellNumber, frip, csvDir)
  covFileBedGraph <- newExportFragsList(popFrags, cellNumber, frip, TxDb)
  print(covFileBedGraph)
    gc()
  covFileBedGraph
})



globalNCores <- 35
############################################################
#      HOMER


Sys.setenv(
  PATH=paste(
    Sys.getenv("PATH"), 
    ":/home/jupyter/scMACS_Analysis/HOMER/bin/", 
    sep="")
)

homerOutDir <- "./homer_simulations_results"
dir.create(homerOutDir)

homer_run_simulated(covFileBedGraphList, rep = 1, homerOutDir)

############################################################
#      MACS2

globalNCores <- 35
gc()
macs2OutDir <- "./macs2__simulations_results"
dir.create(macs2OutDir)

macs2_run_simulated(covFileBedGraphList, rep = 1, macs2OutDir)


############################################################
#      MOCHA
mochaOutFile <- "mocha_simulations_results"
dir.create(mochaOutFile)
MOCHA_run_simulated(cellNumber, iteration1, csvDir = csvDir, outDir = mochaOutFile)


############################################################
#      Analysis of Peak calls
############################################################

homerFiles <- grep("^homer_simulations_results/", list.files(pattern = 'peaks.txt', recursive = TRUE), value = TRUE)
hPeaks <- pbapply::pblapply(cl = 10, X = homerFiles, function(XX){

    tmp1 <- read.table(XX)
    tmp1 <- tmp1[,-1]
    colnames(tmp1) = c('seqnames', 'start', 'end', 
                'strand', 'Num1', 'Num2', 'Num3')
    tilePeaks(tmp1, blackList)
    #makeGRangesFromDataFrame(tmp1)
    })
names(hPeaks) <- gsub(".txt|homer_simulations_results/|/peaks", "", sub(".*_results/", "", homerFiles))
saveRDS(hPeaks, 'ProcessedPeaks_HOMER_v2.rds')

macsFiles <- grep("^macs2__simulations_results/", list.files(pattern = '.broadPeak', recursive = TRUE), value = TRUE)
mPeaks <- pbapply::pblapply(cl = 10, X = macsFiles, function(XX){

    tmp1 <- read.table(XX)
    colnames(tmp1) = c('seqnames', 'start', 'end', 
                'File', 'Num1', 'Num2', 'Num3', 'Num4', 'Num5')
    tilePeaks(tmp1, blackList)

    })
names(mPeaks) <- gsub(".broadPeak|^macs2__simulations_results/", "", sub(".*_results/", "", macsFiles))
saveRDS(mPeaks, 'ProcessedPeaks_MACS2_v2.rds')

mochaPeaks <- readRDS("mocha_simulations_results/tileGrangesList.rds")
mocPeaks <- lapply(mochaPeaks, function(XX){
            tmp <- dplyr::filter(as.data.frame(XX), !is.na(seqnames))
            plyranges::filter(makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE), peak)
    })
names(mocPeaks) <- paste('SimulatedData', names(mochaPeaks), 'peaks', sep = '_')



## Get the true peak set
fullLoc <- readRDS('TruePeakSet.rds')
truePeaks <- lapply(fullLoc, function(XX) MOCHA::GRangesToString(plyranges::filter(XX, isPeak)))
names(truePeaks) <- paste('Round', c(1:10), sep = "_")

library(plyranges)
## Let's calculate the accuracy rate
testResults <- function(newPeaks, truePeaks, type = 'truePos'){

    subNames <- paste(gsub("_peaks", "", gsub(".*Round", "Round", names(newPeaks))), "$", sep = '')
    allOut <- unlist(mclapply(mc.cores = 10, seq_along(newPeaks), function(x){
        
        whichTrue <- truePeaks[[which(grepl(subNames[x], names(truePeaks)))]]
        if(type == 'truePos'){
             sum(MOCHA::GRangesToString(newPeaks[[x]]) %in% whichTrue)
        }else if(type == 'falsePos'){
            sum(!MOCHA::GRangesToString(newPeaks[[x]]) %in% whichTrue)

        }else if(type == 'falseNeg'){
            sum(!whichTrue %in% MOCHA::GRangesToString(newPeaks[[x]]))
        }
        }))
        gc()
    return(allOut)
}
                            
## Function for calculating the number of peaks called

sumDF <- data.frame(TruePositives = c(testResults(hPeaks,truePeaks, type = 'truePos'),
                            testResults( mPeaks,truePeaks, type = 'truePos'),
                             testResults(mocPeaks, truePeaks, type = 'truePos')
                            ),
           FalseNegatives =  c(testResults(hPeaks,truePeaks, type =  'falseNeg'),
                             testResults( mPeaks,truePeaks, type = 'falseNeg'),
                              testResults(mocPeaks, truePeaks, type = 'falseNeg')
                              ),
           FalsePositives =  c(testResults(hPeaks, truePeaks, type = 'falsePos'),
                              testResults( mPeaks, truePeaks, type = 'falsePos'),
                              testResults(mocPeaks, truePeaks, type = 'falsePos')
                              ),
            totalOpenTiles = c(lengths(hPeaks), lengths(mPeaks), lengths(mocPeaks)),
            Sample = c(names(hPeaks), names(mPeaks),names(mocPeaks)),
            Method = c(rep('HOMER', length(hPeaks)), rep('MACS2', length(mPeaks)), rep('MOCHA', length(mocPeaks)))
    )
sumDF <- dplyr::mutate(sumDF, CellNumber = as.numeric(gsub("_Round.*|.*CellNumber_","", Sample)),
                       Iteration = as.numeric(gsub(".*Round_|_peaks","", Sample)))
write.csv(sumDF, 'PeakCalls_AbsoluteCountSummary.csv')
peakNumber = unique(lengths(truePeaks))
                              
## Let's plot it! -- in log10 base? CellNumber?
pdf('SimulatedResults_AbsoluteCounts_v2.pdf')

ggplot(sumDF, 
       aes(x = CellNumber, y = totalOpenTiles, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + ggtitle('Total Open Tiles') +
        geom_hline(yintercept = peakNumber , color = 'red') +
        xlab('Number of Simulated Cells') + ylab('Number of Open Tiles Identified')

ggplot(sumDF, 
       aes(x = CellNumber, y = TruePositives, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + ggtitle('True Positives') +
        geom_hline(yintercept = peakNumber , color = 'red') +
       ylab('Number of True Open Tiles Identified')

ggplot(sumDF, 
       aes(x = CellNumber, y = FalsePositives, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw()+ scale_x_log10() + scale_y_log10() + ggtitle('False Positives')+
        ylab('Number of Open Tiles Incorrectly Identified')
ggplot(sumDF,
       aes(x = CellNumber, y = FalseNegatives, color = Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + ggtitle('False Negatives')+
        ylab('Number of Open Tiles Incorrectly Missed')

tidyr::pivot_longer(sumDF, cols = c('totalOpenTiles', 'TruePositives', 'FalsePositives', "FalseNegatives"), 
                    names_to = 'Metric', values_to = 'PeakNumber') %>%
    ggplot(.,
       aes(x = CellNumber, y = PeakNumber, color = Method, linetype = Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + 
        facet_wrap(~ Metric, nrow = 2, ncol = 2, scales = 'free_y') 
dev.off()

## Let's do this as a percentage....

sumDF2 <- dplyr::mutate(sumDF, 
                        Precision = TruePositives/totalOpenTiles, 
                        Recall = 1 - FalseNegatives/(peakNumber),
                        FalsePositives = FalsePositives/totalOpenTiles)
sumDF2 <- dplyr::mutate(sumDF2, F1 = 2*(Precision * Recall)/(Precision + Recall))

write.csv(sumDF2, 'PeakCalls_Summary_v2.csv')
                                       

pdf('SimulatedResults_Percentages_v2.pdf')
ggplot(sumDF2, 
       aes(x = CellNumber, y = Precision, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() +  geom_smooth()  + theme_bw() + scale_x_log10() + 
        ylab(' Precision (True Positives / All Tiles)') + ggtitle('Peak calling Precision') 

ggplot(sumDF2,
       aes(x = CellNumber, y = Recall, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() +  geom_smooth()  + theme_bw() + scale_x_log10() + 
        ylab('Recall (Called Peaks / True Peaks)')+ ggtitle('True Peak Recall') 

ggplot(sumDF2,
       aes(x = CellNumber, y = FalsePositives, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + 
        ylab('False Positives /All Called Tiles')+ ggtitle('False Positive Rate') 

ggplot(sumDF2, 
       aes(x = CellNumber, y = F1, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() +  geom_smooth() + theme_bw()+ scale_x_log10() +
        ylab('F1 Score')+ ggtitle('Open Chromatin Model Performance Measured by F-Score') 

tidyr::pivot_longer(sumDF2, cols = c('Precision', 'Recall', 'FalsePositives', "F1"), 
                    names_to = 'Metric', values_to = 'PeakNumber') %>%
    ggplot(.,
       aes(x = CellNumber, y = PeakNumber, color = Method, linetype = Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + 
        facet_wrap(~ Metric, nrow = 2, ncol = 2) 
dev.off()
                              
pdf('SimulatedResults_FigurePlots.pdf')

ggplot(dplyr::filter(sumDF, CellNumber > 50), 
       aes(x = CellNumber, y = totalOpenTiles, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + 
        scale_y_log10() + ggtitle('Total Open Tiles') +
        geom_hline(yintercept = peakNumber , color = 'red') +
        xlab('Number of Simulated Cells') + ylab('Number of Open Tiles Identified')
                              
ggplot(dplyr::filter(sumDF2, CellNumber > 50), 
       aes(x = CellNumber, y = Recall, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + 
        ggtitle('Recall') +
        xlab('Number of Simulated Cells') + ylab('Number of True Open Tiles Identified')
                              
ggplot(dplyr::filter(sumDF2, CellNumber > 50), 
       aes(x = CellNumber, y = FalsePositives, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() +
                              ggtitle('Incorrect Peak Calls')  +
        xlab('Number of Simulated Cells') + ylab('Number of Incorrect Open Tiles Identified')

ggplot(dplyr::filter(sumDF2, CellNumber > 50), 
       aes(x = CellNumber, y = F1, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() +  geom_smooth() + theme_bw()+ scale_x_log10() +
        ylab('F1 Score')+ ggtitle('Open Chromatin Model Performance Measured by F-Score') 
dev.off()


############################################################
#      MACS2 at higher numbers
############################################################

library(MOCHA)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

cellNumber = c(100, 200, 500, 1000, 2000, 3000, 5000, 10000, 50000)
iteration1 = c(12:14)

grid <- expand.grid(cellNumber, iteration1)
                              

cl = parallel::makeCluster(35)

parallel::clusterEvalQ(cl, {
        library('BSgenome.Hsapiens.UCSC.hg38')
    library('GenomicRanges')
    library(MOCHA)
      })


truePeaks1 <- readRDS('TruePeakSet.rds')[[1]]
                        
#Make into a list and join with truePeakList
iterList1 <- split(grid,seq(nrow(grid)))
iterList2 <- mapply(append, iterList1, list(truePeaks1) , SIMPLIFY=FALSE)

pbapply::pblapply(cl = cl, X = iterList2, simulateAndWriteFragments)

simulateAndWriteFragments <- function(YY, force = FALSE){
    cellNum2 = as.numeric(YY[[1]])
    y = YY[[2]]
    locPeaks <- YY[[3]]

    fileName = paste('SimulatedFragments2/SimulatedData_CellNumber',  cellNum2, 'Round',y, sep = '_')
    #existingFiles <- list.files(pattern = paste("^",gsub(".*/","", fileName), ".csv", sep =''), recursive = TRUE)
    fragTmp <- MOCHA::simulateFragments(nCells = cellNum2, meanFragsPerCell = 4000, fragThreshold = 1000, 
                                            allLocationsGR = locPeaks, FRIP = 0.95, meanLengths = c(75, 200), 
                                            lengthProbability = c(0.9, 0.1))
            
        
    write.csv(as.data.frame(fragTmp), paste(fileName,'.csv', sep = ''))
    rm(fragTmp)
  
}

stopCluster(cl)
gc()

studySignal <-  4000
study_prefactor <- 3668 / studySignal

csvDir <- "/home/jupyter/scMACS_Analysis/simulatedData/SimulatedFragments2"

globalNCores <- 1

# Generate coverage bedgraphs from simulated data up front
# for MACS2 and HOMER

covFileBedGraphList2 <- mclapply(mc.cores =10, X = split(grid,seq(nrow(grid))), function(XX){
  cellNumber <- as.numeric(XX[1])
  frip <- XX[2]

  popFrags <- getSimulatedFrags(cellNumber, frip, csvDir)
  covFileBedGraph <- newExportFragsList(popFrags, cellNumber, frip, TxDb)
  print(covFileBedGraph)
    gc()
  covFileBedGraph
})
                              
macs2OutDir <- "./macs2__simulations_results"
macs2_run_simulated(covFileBedGraphList2, rep = 1, macs2OutDir)

globalNCores <- 35
                              
macsFiles2 <- grep("^macs2__simulations_results/", list.files(pattern = '.broadPeak', recursive = TRUE), value = TRUE)
macsFiles2 <- grep('12|13|14', macsFiles2, value = TRUE)
mPeaks2 <- pbapply::pblapply(cl = 10, X = macsFiles2, function(XX){

    tmp1 <- read.table(XX)
    colnames(tmp1) = c('seqnames', 'start', 'end', 
                'File', 'Num1', 'Num2', 'Num3', 'Num4', 'Num5')
    tilePeaks(tmp1, blackList)

    })
names(mPeaks2) <- gsub(".broadPeak|^macs2__simulations_results/", "", sub(".*_results/", "", macsFiles2))
saveRDS(mPeaks2, 'ProcessedPeaks_MACS2_v3.rds')

truePeaks2 <- MOCHA::GRangesToString(plyranges::filter(truePeaks1, isPeak))

## Let's calculate the accuracy rate
testResults2 <- function(newPeaks, truePeaks, type = 'truePos'){

    subNames <- paste(gsub("_peaks", "", gsub(".*Round", "Round", names(newPeaks))), "$", sep = '')
    allOut <- unlist(mclapply(mc.cores = 10, seq_along(newPeaks), function(x){
        
        whichTrue <- truePeaks2
        if(type == 'truePos'){
             sum(MOCHA::GRangesToString(newPeaks[[x]]) %in% whichTrue)
        }else if(type == 'falsePos'){
            sum(!MOCHA::GRangesToString(newPeaks[[x]]) %in% whichTrue)

        }else if(type == 'falseNeg'){
            sum(!whichTrue %in% MOCHA::GRangesToString(newPeaks[[x]]))
        }
        }))
        gc()
    return(allOut)
}
                              
sumDF3 <- data.frame(TruePositives = testResults2(mPeaks2, truePeaks2, type = 'truePos'),
           FalseNegatives =  testResults2(mPeaks2, truePeaks2, type = 'falseNeg'),
           FalsePositives = testResults2(mPeaks2, truePeaks2, type = 'falsePos'),
            totalOpenTiles = lengths(mPeaks2),
            Sample =  names(mPeaks2),
            Method = rep('MACS2', length(mPeaks2))
    ) %>%
    dplyr::mutate(CellNumber = as.numeric(gsub("_Round.*|.*CellNumber_","", Sample)),
                       Iteration = as.numeric(gsub(".*Round_|_peaks","", Sample))) %>%
     dplyr::mutate(Precision = TruePositives/totalOpenTiles, 
                        Recall = 1 - FalseNegatives/(length(truePeaks2)),
                        FalsePositivesRate = FalsePositives/totalOpenTiles) %>%
        dplyr::mutate(F1 = 2*(Precision * Recall)/(Precision + Recall)) %>%
        tidyr::pivot_longer(cols = c('totalOpenTiles', 'TruePositives', 'FalsePositives', "FalseNegatives",
                                      'Precision', 'Recall', 'FalsePositivesRate', 'F1'), 
                    names_to = 'Metric', values_to = 'PeakNumber')
pdf('SimulatedResults_MACS2_BulkCellNumbers.pdf')
                              
ggplot(dplyr::filter(sumDF3, Metric %in% c('F1', "Recall", 'FalsePositiveRate', 'totalOpenTiles')), 
       aes(x = CellNumber, y = PeakNumber, color = Method, linetype = Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + 
        facet_wrap(~ Metric, nrow = 2, ncol = 2, scales = 'free_y') 
                              
dev.off()

