ThemeMain <- theme(
  plot.margin = unit(c(0, 0, 0, 0), "npc"),
  plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
  axis.title.x = element_text(size = 20, color = "black"),
  axis.title.y = element_text(size = 20, color = "black"),
  axis.text.x = element_text(size = 16, color = "black", angle = 90),
  axis.text.y = element_text(size = 14, color = "black"),
  strip.text.x = element_text(size = 20, color = "black"),
  strip.text.y = element_text(size = 20, color = "black"),
  strip.text = element_text(size = 20, color = "black")
)

# Define a helper to subset the ArchR project to our cellQuants
# and extract fragments for this subsetted project
getSubsetArchRFrags <- function(covidArchR, numSamples) {
  
  # Subset cells to a random selection of samples, size numSamples
  set.seed(2022)
  selectedSamples <- sample(unique(covidArchR$Sample), numSamples)
  idxSample <- BiocGenerics::which(covidArchR$Sample %in% selectedSamples)
  subsetArchR <- covidArchR[covidArchR$cellNames[idxSample], ]
  subsetArchR

  # Extract fragments
  popFrags <- MOCHA::getPopFrags(
    ArchRProj = subsetArchR,
    cellPopLabel = "CellSubsets",
    cellSubsets = "ALL",
    numCores = globalNCores,
    verbose = FALSE
  )
  rm(subsetArchR)
  
  # return fragments
  popFrags
}

# This function creates one bedgraph per sample, and returns 
# a list of filenames.
newExportFragsList <- function(popFrags, rep, numSamples, TxDb){
  
  dir.create("./coverage_files")
  
  # Extract norm factor for each group
  # getPopFrags with normMethod=nCells
  # returns (1000/nCells) as the norm factor
  groups <- gsub("__.*", "", names(popFrags))
  cellNorm <- as.numeric(gsub(".*_", "", names(popFrags)))
  names(cellNorm) <- groups
  
  # We only use one cell population
  # for this benchmark
  normFactors <- lapply(seq_along(cellNorm), function(x) {
    normFactor <- cellNorm[groups[x]]
  })
  
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
  covFileBedGraphList <- mclapply(seq_along(coverage_gr), function(i){
    
    covGR <- coverage_gr[[i]]
    sampleName <- gsub(".*#", "",gsub("__.*", "", names(coverage_gr)[[i]]))
    
    covFileBedGraph <- paste(
      "rep", rep, "nsamples", numSamples, sampleName,
      "CD14.bedGraph", 
    sep = "_")
    covFileBedGraph <- file.path("./coverage_files", covFileBedGraph)
    
    plyranges::write_bed_graph(
      covGR,
      covFileBedGraph, 
      index = FALSE
    )
    
    covFileBedGraph
  })
  
  covFileBedGraphList
}
  


# Helper to run repeated measures
# These should not be parallelized to avoid
# influence of resource allocations
runBenchmark <- function(peakCaller, outDir, numRepeats, name){
  dir.create(outDir)
  message("Now benchmarking: ", basename(outDir))
  runtimes <- lapply(
    unique(sampleQuants),
    function(nsamples) {
      gc()
      message("Number of Samples: ", nsamples)
      lapply(
        1:numRepeats,
        function(x){
          peakCaller(numSamples = nsamples, rep = x, outDir)
        }
      )
    }
  )

  runtimes.df <- as.data.frame(rbindlist(lapply(
    runtimes[0:numRepeats+1],
    function(x) {
      rbindlist(x)
    }
  )))

  fwrite(runtimes.df, paste(basename(outDir), "results_", name, ".csv", sep=''))
  runtimes.df
}

############################################################
# MACS2 
############################################################
macs2_measure_runtimes <- function(numSamples, rep = 1, macs2OutDir) {
  
  popFrags <- getSubsetArchRFrags(covidArchR, numSamples)
  tic()
  #~~~~~~~~~~~~~~~~ Parallelize over 10 cores
  covFileBedGraphList <- newExportFragsList(popFrags, rep, numSamples, TxDb)
  rm(popFrags)
  
  mclapply(covFileBedGraphList, function(fname){
    cmd_sample <- paste(
      "macs2 callpeak -t",
      fname,
      "-g hs -f BED --nolambda --shift -75 --extsize 150 --broad",
      "--nomodel -n",
      file.path(
        macs2OutDir, tools::file_path_sans_ext(basename(fname))
      )
    )
    system(cmd_sample)
    system(paste("rm -R", file.path(
      macs2OutDir, tools::file_path_sans_ext(basename(fname))
    )))
  }, mc.cores = globalNCores)
  
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic

  res.df <- data.frame(
    Time = runtime_secs,
    Nsamples = numSamples
  )
  write.csv(res.df, file.path(
    macs2OutDir, paste(
      "rep", rep, "numSamples", numSamples, 'results.csv', sep="_"
  )))
  gc()
  
  lapply(covFileBedGraphList, file.remove)
  return(res.df)
}

############################################################
# HOMER
############################################################
homer_measure_runtimes <- function(numSamples, rep = 1, homerOutDir) {
  
  popFrags <- getSubsetArchRFrags(covidArchR, numSamples)
  
  tic()
  #~~~~~~~~~~~~~~~~ Parallelize over 10 cores
  covFileBedGraphList <- newExportFragsList(popFrags, rep, numSamples, TxDb)
  rm(popFrags)
  
  mclapply(covFileBedGraphList, function(covFileBedGraph){
    
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
      sep=' '
    )
    system(callPeaks_cmd)
    system(paste("rm -R", tagDir))
    
  }, mc.cores = globalNCores)
  
  #system(paste("rm -R", tagDir))
  
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic
  
  res.df <- data.frame(
    Time = runtime_secs,
    Nsamples = numSamples
  )
  write.csv(res.df, file.path(
    homerOutDir, paste(
      "rep", rep, "numSamples", numSamples, 'results.csv', sep="_"
    )
  ))
  gc()
  lapply(covFileBedGraphList, file.remove)
  return(res.df)
}


############################################################
# MOCHA
############################################################
MOCHA_measure_runtime <- function(numSamples, rep = NULL, outDir = NULL) {

  # Get sample-specific fragments for MOCHA
  popFrags <- getSubsetArchRFrags(covidArchR, numSamples)
  
  tic()
  #~~~~~~~~~~~~~~~~ Parallelize over 10 cores
  totalFrags <- as.integer(sapply(popFrags, length))
  tileList <- mclapply(popFrags, function(sampleFrags){
    tilesGRanges <- MOCHA:::callTilesBySample(
      blackList = blackList,
      returnAllTiles = TRUE,
      totalFrags = totalFrags,
      fragsList = sampleFrags,
      verbose = FALSE,
      StudypreFactor = study_prefactor
    )
  }, mc.cores = globalNCores)
  rm(popFrags)
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic
  
  res.df <- data.frame(
    Time = runtime_secs,
    Nsamples = numSamples
  )
  write.csv(res.df, file.path(
    mochaOutFile, paste(
      "rep", rep, "numSamples", numSamples, 'results.csv', sep="_"
  )))
  gc()
  return(res.df)
}
