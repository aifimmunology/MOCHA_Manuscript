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
getSubsetArchRFrags <- function(covidArchR, numCells, sampleSpecific = FALSE) {
  # Subset cells to a random selection of size numCells
  # This will get fragments across multiple samples
  set.seed(2022)
  idxCells <- sample(1:TotalNCells, size = numCells)
  subsetArchR <- covidArchR[idxCells, ]

  # Extract fragments
  popFrags <- MOCHA::getPopFrags(
    ArchRProj = subsetArchR,
    metaColumn = "CellSubsets",
    cellSubsets = "ALL",
    NormMethod = "nCells",
    numCores = 10,
    sampleSpecific = sampleSpecific,
    verbose = FALSE
  )
  rm(subsetArchR)
  
  # return fragments
  popFrags
}

newExportFrags <- function(popFrags, rep, numCells, TxDb){
  
  dir.create("./coverage_files")
  # Export fragments to bedgraph
  # for use by Homer and MACS2
  covFileBedGraph <- paste(
    "tmp", rep, numCells, 
    "CD14_perThousandCells.bedGraph", 
  sep = "_")
  covFileBedGraph <- file.path("./coverage_files", covFileBedGraph)
  
  # Extract norm factor for each group
  # getPopFrags with normMethod=nCells
  # returns (1000/nCells) as the norm factor
  groups <- gsub("_.*", "", names(popFrags))
  cellNorm <- as.numeric(gsub(".*_", "", names(popFrags)))
  names(cellNorm) <- groups
  
  # Have to make this a list to generalize to multiple cell
  # populations, though we only use one cell population
  # for this benchmark
  normFactors <- lapply(seq_along(cellNorm), function(x) {
    normFactor <- cellNorm[groups[x]]
  })
  
  coverage_gr <- getCoverage(
   popFrags, 
   normFactors, 
   TxDb, 
   filterEmpty = FALSE, numCores = 1, verbose = FALSE
  )
  
  # export to bedGraph with plyranges
  plyranges::write_bed_graph(
    IRanges::stack(as(coverage_gr, "GRangesList")), 
    covFileBedGraph, 
    index = FALSE
  )
  return(covFileBedGraph)
}


# Helper to run repeated measures
runBenchmark <- function(peakCaller, outDir, numRepeats, name){
  dir.create(outDir)
  message("Now benchmarking: ", basename(outDir))
  runtimes <- lapply(
    unique(cellQuants),
    function(TIME) {
      gc()
      message("Number of Cells: ", TIME)
      mclapply(
        1:numRepeats,
        function(x){
          peakCaller(numCells = TIME, rep = x, outDir)
        },
        mc.cores = 2
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
macs2_measure_runtimes <- function(numCells, rep = 1, macs2OutDir) {
  
  popFrags <- getSubsetArchRFrags(
    covidArchR, numCells, sampleSpecific = FALSE
  )
  tic()
  #~~~~~~~~~~~~~~~~
  covFileBedGraph <- newExportFrags(popFrags, rep, numCells, TxDb)
  rm(popFrags)

  cmd_sample <- paste(
    "macs2 callpeak -t",
    covFileBedGraph,
    "-g hs -f BED --nolambda --shift -75 --extsize 150 --broad",
    "--nomodel -n",
    file.path(macs2OutDir, tools::file_path_sans_ext(basename(covFileBedGraph)))
  )
  system(cmd_sample)
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic

  res.df <- data.frame(
    Time = runtime_secs,
    Ncells = numCells
  )
  gc()
  file.remove(covFileBedGraph)
  return(res.df)
}

############################################################
# HOMER
############################################################
homer_measure_runtimes <- function(numCells, rep = 1, homerOutDir) {
  
  popFrags <- getSubsetArchRFrags(
    covidArchR, numCells, sampleSpecific = FALSE
  )
  
  tic()
  #~~~~~~~~~~~~~~~~
  covFileBedGraph <- newExportFrags(popFrags, rep, numCells, TxDb)
  rm(popFrags)
  
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
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic
  
  res.df <- data.frame(
    Time = runtime_secs,
    Ncells = numCells
  )
  gc()
  file.remove(covFileBedGraph)
  return(res.df)
}


############################################################
# MOCHA
############################################################
MOCHA_measure_runtime <- function(numCells, rep = NULL, outDir = NULL) {

  # Get sample-specific fragments for MOCHA
  frags <- getSubsetArchRFrags(
    covidArchR, numCells, sampleSpecific = FALSE
  )
  
  tic()
  #~~~~~~~~~~~~~~~~
  totalFrags <- as.integer(sapply(frags, length))
  
  tilesGRanges <- MOCHA:::callTilesBySample(
    blackList = blackList,
    returnAllTiles = TRUE,
    totalFrags = totalFrags,
    fragsList = frags[[1]],
    verbose = FALSE,
    StudypreFactor = study_prefactor
  )
  rm(frags)
  #~~~~~~~~~~~~~~~~
  runtime <- toc(log = T)
  runtime_secs <- runtime$toc - runtime$tic
  
  res.df <- data.frame(
    Time = runtime_secs,
    Ncells = numCells
  )
  gc()
  return(res.df)
}
