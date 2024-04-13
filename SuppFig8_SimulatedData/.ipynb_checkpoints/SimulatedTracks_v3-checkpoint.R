############################################################
#      Simulating Single Cell ATAC-seq Data
############################################################

library(MOCHA)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

devtools::load_all('../../MOCHA')

peakNumber = c(150000, 250000, 350000)
fragsPerCell = c(1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)

grid <- expand.grid(peakNumber, fragsPerCell)

## Get all peaksets
locationList <- pbapply::pblapply(peakNumber, function(XX){
    getSimulatedPeakSet(peakNumber = XX, peakSize = 500, peakPeriodicity = 10, 
                                Genome = BSgenome.Hsapiens.UCSC.hg38, allLocations = TRUE)
}, cl =2)


names(locationList) <- peakNumber
saveRDS(locationList, 'TruePeakSet_VariablePeakNumber.rds')

#duplicate to match grid
locationList2 <- locationList[as.character(grid[,1])]
#Make into a list and join with truePeakList
iterList1 <- split(grid,seq(nrow(grid)))
iterList2 <- mapply(append, iterList1, locationList2 , SIMPLIFY=FALSE)

Rounds = c(1:10)

## Create function for simulating and writing fragments to disk
simulateAndWriteFragments <- function(YY, force = FALSE, round = round){
    peakNumber = as.numeric(YY[[1]])
    meanFrags = YY[[2]]
    locPeaks <- YY[[3]]

    fileName = paste('SimulatedFragments/SimulatedData_PeakNumber',  peakNumber, 
                     'MeanFrags',meanFrags, 'Round', round, sep = '_')
    #existingFiles <- list.files(pattern = paste("^",gsub(".*/","", fileName), ".csv", sep =''), recursive = TRUE)
    fragTmp <- MOCHA::simulateFragments(nCells = 250,
                                        meanFragsPerCell = meanFrags, 
                                        fragThreshold = 500, 
                                        allLocationsGR = locPeaks, 
                                        FRIP = 0.95, 
                                        meanLengths = c(75, 200), 
                                        lengthProbability = c(0.9, 0.1))

    write.csv(as.data.frame(fragTmp), paste(fileName,'.csv', sep = ''))
    rm(fragTmp)
  
}

## Now parallelize and run the function. 
cl = parallel::makeCluster(25)

parallel::clusterEvalQ(cl, {
        library('BSgenome.Hsapiens.UCSC.hg38')
    devtools::load_all('../../MOCHA')
    library('GenomicRanges')
      })

for(i in 1:10){
pbapply::pblapply(cl = cl, X = iterList2,
                  simulateAndWriteFragments, round =i)
    gc()
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


source("comparisonsHelperFunctions_v3.R")
library(ArchR)
library(MOCHA)
library(ggplot2)
library(data.table)
library(pbapply)
library(scales)
library(parallel)

csvDir <- "/home/jupyter/scMACS_Analysis/simulatedData2/SimulatedFragments"

globalNCores <- 1

# Generate coverage bedgraphs from simulated data up front
# for MACS2 and HOMER
covFileBedGraphList= list()
for(i in 1:10){
    
    covFileBedGraphList_tmp <- mclapply(mc.cores =20, X = split(grid,seq(nrow(grid))), function(XX){
      peakNumber <- as.numeric(XX[1])
      meanFrags <- XX[2]

      popFrags <- getSimulatedFrags(peakNumber, meanFrags, round = i, csvDir)
      covFileBedGraph <- newExportFragsList(popFrags, peakNumber, meanFrags, round = i, TxDb)
      print(covFileBedGraph)
        gc()
      covFileBedGraph
    })
    gc()
    covFileBedGraphList = append(covFileBedGraphList, covFileBedGraphList_tmp)
}


saveRDS(covFileBedGraphList, 'BedGraphList.rds')

covFileBedGraphList <- readRDS('BedGraphList.rds')


globalNCores <- 10

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

globalNCores <- 10
gc()
macs2OutDir <- "./macs2__simulations_results"
dir.create(macs2OutDir)

macs2_run_simulated(covFileBedGraphList, rep = 1, macs2OutDir)


############################################################
#      MOCHA
mochaOutFile <- "mocha_simulations_results"
dir.create(mochaOutFile)

peakNumber = c(150000, 250000, 350000, 450000)
fragsPerCell = c(1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)

globalNCores <- 15

for(i in 1:10){
    MOCHA_run_simulated(peakNumber, fragsPerCell, round = i, csvDir = csvDir, outDir = mochaOutFile)
}
gc()

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
saveRDS(hPeaks, 'ProcessedPeaks_HOMER_v3.rds')

macsFiles <- grep("^macs2__simulations_results/", list.files(pattern = '.broadPeak', recursive = TRUE), value = TRUE)
mPeaks <- pbapply::pblapply(cl = 10, X = macsFiles, function(XX){

    tmp1 <- read.table(XX)
    colnames(tmp1) = c('seqnames', 'start', 'end', 
                'File', 'Num1', 'Num2', 'Num3', 'Num4', 'Num5')
    tilePeaks(tmp1, blackList)

    })
names(mPeaks) <- gsub(".broadPeak|^macs2__simulations_results/", "", sub(".*_results/", "", macsFiles))
saveRDS(mPeaks, 'ProcessedPeaks_MACS2_v3.rds')

mocPeaks <- list()
for(i in c(1:10)){
    
    mochaPeaks_tmp <- readRDS(paste("mocha_simulations_results/tileGrangesList_Round_",
                                    i,".rds",sep=""))
    mochaPeaks_tmp2 <- lapply(mochaPeaks_tmp, function(XX){
            tmp <- dplyr::filter(as.data.frame(XX), !is.na(seqnames))
            plyranges::filter(makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE), peak)
    })
    names(mochaPeaks_tmp2) <- paste('SimulatedData', names(mochaPeaks_tmp2), 'peaks', sep = '_')
    mocPeaks <- append(mocPeaks, mochaPeaks_tmp2)
}
names(mocPeaks) <- paste('SimulatedData', names(mochaPeaks), 'peaks', sep = '_')

## Get the true peak set
fullLoc <- readRDS('TruePeakSet_VariablePeakNumber.rds')
truePeaks <- lapply(fullLoc, function(XX) MOCHA::GRangesToString(plyranges::filter(XX, isPeak)))
names(truePeaks) <- paste('PeakNumber_', names(truePeaks), sep ='')
                    
                    
library(plyranges)
## Let's calculate the accuracy rate
testResults <- function(newPeaks, truePeaks, type = 'truePos'){

    subNames <- paste(gsub("_MeanFrags.*", "", 
                           gsub(".*PeakNumber", "PeakNumber", names(newPeaks))), "$", sep = '')
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
                              
sumDF <- dplyr::mutate(sumDF, Sample = gsub("_peaks$","", Sample)) %>%
        dplyr::mutate(
        PeakNumber = as.numeric(gsub("_MeanFrags.*|.*PeakNumber_","", Sample)),
        MeanFrags = as.numeric(gsub(".*MeanFrags_|_peaks|_Round.*","", Sample)),
        Round = as.numeric(gsub(".*_Round_","", Sample))) 
                              
write.csv(sumDF, 'PeakCalls_AbsoluteCountSummary.csv')
                        
## Let's do this as a percentage....

sumDF2 <- dplyr::mutate(sumDF, 
                        Precision = TruePositives/totalOpenTiles, 
                        Recall = 1 - FalseNegatives/(PeakNumber),
                        FalsePositives = FalsePositives/totalOpenTiles)
sumDF2 <- dplyr::mutate(sumDF2, F1 = 2*(Precision * Recall)/(Precision + Recall))

write.csv(sumDF2, 'PeakCalls_Summary_v3.csv')
                              
pdf('SimulatedResults_FigurePlots_CellTypes.pdf')

lapply(unique(sumDF$PeakNumber), function(i){
   ggplot(dplyr::filter(sumDF, PeakNumber ==i), 
       aes(x = MeanFrags, y = totalOpenTiles, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + 
        scale_y_log10() + ggtitle(paste('Total Open Tiles for Peakset with ', i, ' Open Regions', sep='')) +
        geom_hline(yintercept = i, color = 'black') +
        geom_vline(xintercept = 3000*0.95, color = 'red') +
        geom_vline(xintercept = 1000, color = 'red') +   
        geom_text(aes(3000*0.95+1250, .4*10^5,label = 'Signac & SnapATAC Thresholds'), color = 'black') + 
        geom_text(aes(1250, 0.4*10^5, label = 'ArchR Threshold'), color = 'black') +
        coord_cartesian(clip = "off") +
        xlab('Sequencing Depth Per Cell') + 
        ylab('Number of Open Tiles Identified') 
})
                              
lapply(unique(sumDF2$PeakNumber), function(i){
      ggplot(dplyr::filter(sumDF2, PeakNumber ==i), 
       aes(x = MeanFrags, y = Recall, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() + scale_y_log10() + 
        ggtitle(paste('Recall for Peakset with ', i, ' Open Regions', sep='')) +
          geom_vline(xintercept = 3000*0.95, color = 'red') +
        geom_vline(xintercept = 1000, color = 'red') +   
        geom_text(aes(3000*0.95+1250, .4,label = 'Signac & SnapATAC Thresholds'), color = 'black') + 
        geom_text(aes(1250, 0.4, label = 'ArchR Threshold'), color = 'black') +
        coord_cartesian(clip = "off") +
        xlab('Sequencing Depth Per Cell') + 
        ylab('Recall')
})
                              
lapply(unique(sumDF2$PeakNumber), function(i){
      ggplot(dplyr::filter(sumDF2, PeakNumber ==i), 
       aes(x = MeanFrags, y = FalsePositives, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() + geom_smooth() + theme_bw() + scale_x_log10() +
      geom_vline(xintercept = 3000*0.95, color = 'red') +
        geom_vline(xintercept = 1000, color = 'red') +   
        geom_text(aes(3000*0.95+1250, .4,label = 'Signac & SnapATAC Thresholds'), color = 'black') + 
        geom_text(aes(1250, 0.4, label = 'ArchR Threshold'), color = 'black') +
        coord_cartesian(clip = "off") +
        ggtitle(paste('Incorrect Peak Calls for Peakset with ', i, ' Open Regions', sep='')) +
        xlab('Sequencing Depth Per Cell') + 
        ylab('False Positive Rate')
})                           
lapply(unique(sumDF2$PeakNumber), function(i){
      ggplot(dplyr::filter(sumDF2, PeakNumber ==i), 
       aes(x = MeanFrags, y = F1, color =  Method,
           shape = Method, group = Method)) + 
        geom_point() +  geom_smooth() + theme_bw()+ scale_x_log10() +
      geom_vline(xintercept = 3000*0.95, color = 'red') +
        geom_vline(xintercept = 1000, color = 'red') +   
        geom_text(aes(3000*0.95+1250, .4,label = 'Signac & SnapATAC Thresholds'), color = 'black') + 
        geom_text(aes(1250, 0.4, label = 'ArchR Threshold'), color = 'black') +
        coord_cartesian(clip = "off") +
        ylab('F1 Score')+ 
        ggtitle(paste('F-Scores for Peakset with ', 
                i, ' Open Regions', sep=''))
})
    
dev.off()

