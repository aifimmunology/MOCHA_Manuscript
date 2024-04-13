# Get the simulated fragments as a GRanges from a given cellnumber and frip score
getSimulatedFrags <- function(cellNumber, frip, csvDir = "/home/jupyter/cache/8115ccee-b1d2-4345-b14b-38c394bcc23c/simulatedData/"){
  #cache/8115ccee-b1d2-4345-b14b-38c394bcc23c/simulatedData/SimulatedData_CellNumber_1500_FRIPScore_0.6.csv
  df <- read.csv(paste0(
    csvDir, "/",
    "SimulatedData",
    "_CellNumber_",
    cellNumber, 
    "_Round_",
    frip,
    ".csv"
  ))
  fragGR <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  rm(df)
  list(fragGR)
}

# # Get the path to a coverage bigwig from a given cellnumber and frip score
# getCoverageBWPath <- function(cellNumber, frip, covDir = "/home/jupyter/cache/8115ccee-b1d2-4345-b14b-38c394bcc23c/simulatedData/CoverageFiles"){
#   
#   #cache/8115ccee-b1d2-4345-b14b-38c394bcc23c/simulatedData/CoverageFiles/SimulatedData_CellNumber_10_Round_2_Coverage.bw
#   
#   bwpath <- file.path(covDir, paste(
#     "SimulatedData",
#     "CellNumber",
#     cellNumber, 
#     "Round",
#     frip,
#     "Coverage.bw", sep="_"
#   ))
#   
#   if(file.exists(bwpath)){
#     return(bwpath)
#   } else {
#     stop(bwpath, "does not exist")
#   }
#   
# }

# This function creates one bedgraph per sample, and returns 
# a list of filenames from a given cellnumber and frip score
newExportFragsList <- function(popFrags, cellNumber, frip, TxDb){
  
  dir.create("./coverage_files")
  
  # calculate norm factor for each group
  # getPopFrags with normMethod=nCells returns (1000/nCells) as the norm factor
  normFactors <- list(1000/cellNumber)
  
  # source("/home/jupyter/MOCHA/R/getCoverage.R")
  coverage_gr <- MOCHA::getCoverage(
   popFrags, 
   normFactors,
   TxDb,
   cl = globalNCores,
   filterEmpty = FALSE, verbose = FALSE
  )
  
  # export to bedGraph with plyranges
  # Export fragments to bedgraph
  # for use by Homer and MACS2
  covGR <- coverage_gr[[1]]
  
  covFileBedGraph <- paste0(
    "SimulatedData_", "CellNumber_", cellNumber, "_Round_", frip, ".bedGraph"
  )
  covFileBedGraph <- file.path("./coverage_files", covFileBedGraph)
  
  plyranges::write_bed_graph(
    covGR,
    covFileBedGraph, 
    index = FALSE
  )
  
  covFileBedGraph
}

############################################################
# MACS2 
############################################################
macs2_run_simulated <- function(covFileBedGraphList, rep = 1, macs2OutDir) {
  
  #~~~~~~~~~~~~~~~~ Parallelize over 10 cores
  outFileList <- mclapply(covFileBedGraphList, function(covFileBedGraph){
    
    outFile <- file.path(macs2OutDir, tools::file_path_sans_ext(basename(covFileBedGraph)))
    
    cmd_sample <- paste(
      "macs2 callpeak -t",
      covFileBedGraph,
      "-g hs -f BED --nolambda --shift -75 --extsize 150 --broad",
      "--nomodel -n",
      outFile
    )
    
    system(cmd_sample)
    outFile
  }, mc.cores = globalNCores)
  
  #~~~~~~~~~~~~~~~~
  gc()
  outFileList
}

############################################################
# HOMER
############################################################
homer_run_simulated <- function(covFileBedGraphList, rep = 1, homerOutDir) {
  
  #~~~~~~~~~~~~~~~~ Parallelize over 10 cores
  outFileList <- mclapply(covFileBedGraphList, function(covFileBedGraph){
    
    tagDir <- file.path(
      homerOutDir, 
      tools::file_path_sans_ext(basename(covFileBedGraph)), 
      sep=''
    )

    tagdir_cmd <- paste(
      'makeTagDirectory',
      tagDir,
      covFileBedGraph
    )
    system(tagdir_cmd)

    callPeaks_cmd <- paste(
      'findPeaks',
      tagDir,
      '>',
      file.path(tagDir, 'peaks.txt'),
      '-style histone',
      '-size 150',
      sep=' '
    )
    system(callPeaks_cmd)
    
    file.path(tagDir, 'peaks.txt')
  }, mc.cores = globalNCores)
  #~~~~~~~~~~~~~~~~
  
  gc()
  outFileList
}


############################################################
# MOCHA
############################################################
MOCHA_run_simulated <- function(cellNumbers, FRIPScores, csvDir, outDir = mochaOutFile) {

  # Get a matrix with all combinations of cell# and FRIP
  grid <- expand.grid(cellNumbers, FRIPScores)
  
  #~~~~~~~~~~~~~~~~ Parallelize
  # Lapply over rows of the matrix
  tilesGRangesList <- mclapply(seq_len(dim(grid)[1]), function(x){

    # Get sample-specific fragments for MOCHA
    cellNumber <- grid[x, 1]
    frip <- grid[x, 2]
    popFrags <- getSimulatedFrags(cellNumber, frip, csvDir)[[1]]
    totalFrags <- length(popFrags)

    popFrags$cell <- popFrags$CellID
    
    # Call tiles
    tilesGRanges <- MOCHA:::callTilesBySample(
      blackList = blackList,
      returnAllTiles = TRUE,
      totalFrags = totalFrags,
      fragsList = popFrags,
      cellCol = "CellID", # This col exists in simulated data if converted to GRanges with keep.extra.columns=TRUE
      verbose = FALSE,
      StudypreFactor = study_prefactor
    )
    rm(popFrags)
    
    tilesGRanges <- list(tilesGRanges)
    names(tilesGRanges) <- c(paste0("CellNumber_", cellNumber, "_Round_", frip))

    tilesGRanges
      
  }, mc.cores = globalNCores)
  #~~~~~~~~~~~~~~~~
  
  tilesGRangesList <- unlist(tilesGRangesList, recursive = FALSE)
  saveRDS(tilesGRangesList, file.path(mochaOutFile, "tileGrangesList.rds"))
  
  gc()
  file.path(mochaOutFile, "tileGrangesList.rds")
}


################### Helper functions

## Function for tiling peaks

tilePeaks <- function(mat, blackList){
    require(data.table)
    require(plyranges)
    diff_start <- mat$start %% 500
    diff_end <- mat$end %% 500

    mat$new_start <- mat$start
    idx_start <- which(diff_start <= 75)
    mat$new_start[idx_start] <- mat$start[idx_start] + 75
    
    mat$new_end <- mat$end 
    idx_end <- which(diff_end <= 75)
    mat$new_end[idx_end] <- mat$new_end[idx_end] - 75
        
    colnames(  mat )[2:3] <-c('old_start','old_end')
    tmp <- makeGRangesFromDataFrame(mat, start.field='new_start',
                                 end.field='new_end', keep.extra.columns=TRUE)
    binSize=500
    TotalRangesFilt<- plyranges::setdiff_ranges(tmp, blackList)

    RangeBins <- plyranges::stretch(plyranges::anchor_end(TotalRangesFilt),
                                  extend = GenomicRanges::start(TotalRangesFilt)%% binSize)

    FinalBins <- plyranges::stretch(plyranges::anchor_start(RangeBins),
                                  extend = (binSize - end(RangeBins)%% binSize)) %>%
    plyranges::reduce_ranges() %>% plyranges::slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

    FinalBins$tileID <- paste(FinalBins$seqnames,':',
                           FinalBins$start,'-',
                           FinalBins$end,
                           sep='')

    return(FinalBins)
    
}



