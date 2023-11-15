setwd("MOCHA_revisions/")
install.packages("ggbio")
install.packages("tictoc")
# devtools::install("../MOCHA")

# Run in terminal:
# pip install numpy --upgrade

Sys.setenv(
  PATH=paste(
    Sys.getenv("PATH"), 
    ":/home/jupyter/MOCHA_runtime_comparison/HOMER/bin/", 
    sep="")
)
install.packages('pbapply')
install.packages('RPushbullet')

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

# source("/home/jupyter/MOCHA_fig2/ArchR_Export_ATAC_Data.R")
source("./runtimeHelperFunctions.R")
library(ArchR)
library(MOCHA)
library(ggplot2)
library(tictoc)
library(data.table)
library(pbapply)
library(RPushbullet)
library(scales)
library(parallel)

RPushbullet::pbSetup()

fullCovidArchR <- loadArchRProject("/home/jupyter/FullCovid")
# Extract these from the full covidArchR
blackList <- ArchR::getBlacklist(fullCovidArchR)
cellColData <- ArchR::getCellColData(fullCovidArchR)
studySignal <- stats::median(cellColData$nFrags)
study_prefactor <- 3668 / studySignal # Training median

# Choose a large cell population that we can downsample
cellPop <- "CD14 Mono"

# Subset full ArchR Project to this cell population
idxSample <- BiocGenerics::which(fullCovidArchR$CellSubsets %in% cellPop)
cellsSample <- fullCovidArchR$cellNames[idxSample]
covidArchR <- fullCovidArchR[cellsSample, ]

# Remove full ArchR Project to save RAM
rm(fullCovidArchR)

TotalNCells <- length(cellsSample)
# Define our quantities for downsampling
# cellQuants <- c(
#     min(50000, TotalNCells),
#     min(100000, TotalNCells)
# )
# cellQuants <- c(10, 50, 100, 250, 500, 1000)
# cellQuants <- c(
#   min(2000, TotalNCells),
#   min(3000, TotalNCells),
#   min(5000, TotalNCells),
#   min(10000, TotalNCells),
#   min(15000, TotalNCells)
# )
# cellQuants <- c(
#   min(50000, TotalNCells),
#   min(100000, TotalNCells)
# )

# Remove repeats (if cell pop) does not have enough cells
# to go beyond its sample size
# idxMax <- which.max(cellQuants)
# cellQuants <- cellQuants[1:idxMax]

sampleQuants <- seq(10, 90, by=10)
sampleQuants <- c(70,80,90)
globalNCores <- 10
name <- "parallelruntime_by10samples_70-90"

############################################################
#      MOCHA
############################################################
mochaOutFile <- "mocha_runtimePeaks_byNSamples"

MOCHA_runtimes.df <- runBenchmark(
  peakCaller = MOCHA_measure_runtime,
  outDir = mochaOutFile,
  numRepeats = 5,
  name = name
)
pbPost("note", "MOCHA complete, check output mocha_runtimes.df")

############################################################
#      MACS2
############################################################
macs2OutDir <- "./macs2_runtimePeaks_byNSamples"
dir.create(macs2OutDir)

Macs2_runtimes.df <- runBenchmark(
  peakCaller = macs2_measure_runtimes,
  outDir = macs2OutDir,
  numRepeats = 5,
  name = name
)
pbPost("note", "MACS2 complete, check output Macs2_runtimes.df")

############################################################
#      HOMER
############################################################
homerOutDir <- "./homer_runtimePeaks_byNSamples"
dir.create(homerOutDir)

homer_runtimes.df <- runBenchmark(
  peakCaller = homer_measure_runtimes,
  outDir = homerOutDir,
  numRepeats = 5,
  name = name
)
pbPost("note", "Homer complete, check output homer_runtimes.df")

############################################################
# Plotting
############################################################
Macs2_runtimes.dfa <- read.csv("macs2_runtimePeaks_byNSamplesresults_parallelruntime_by10samples_10-60.csv")
Macs2_runtimes.dfb <- read.csv("macs2_runtimePeaks_byNSamplesresults_parallelruntime_by10samples_70-90.csv")
Macs2_runtimes.df <- rbind(Macs2_runtimes.dfa, Macs2_runtimes.dfb)

MOCHA_runtimes.dfa <- read.csv("mocha_runtimePeaks_byNSamplesresults_parallelruntime_by10samples_10-60.csv")
MOCHA_runtimes.dfb <- read.csv("mocha_runtimePeaks_byNSamplesresults_parallelruntime_by10samples_70-90.csv")
MOCHA_runtimes.df <- rbind(MOCHA_runtimes.dfa, MOCHA_runtimes.dfb)

homerpaths <- Sys.glob("./homer_runtimePeaks_byNSamples/rep_*_numSamples_*0_results.csv")
homerresultlist <- lapply(homerpaths, function(f) {
    read.csv(f)
})
homer_runtimes.df <- do.call(rbind, homerresultlist)
homer_runtimes.df <- homer_runtimes.df[c("Time", "Nsamples")]

Macs2_runtimes.df$Model <- 'MACS2'
MOCHA_runtimes.df$Model <- 'MOCHA'                  
homer_runtimes.df$Model <- 'Homer'

combined_res <- rbind(Macs2_runtimes.df,
                      MOCHA_runtimes.df,
                      homer_runtimes.df
                     )

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
  # for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable to be summarized
# groupnames : vector of column names to be used as grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
 return(data_sum)
}

data <- data_summary(combined_res, varname="Time", groupnames=c("Model", "Nsamples"))
  
pdf('./ss_runtime_parallelv2.pdf')
p<- ggplot(
    data=data,
    aes( x=Nsamples, y=Time, col=Model, fill=Model)
    ) + 
    geom_line() + 
    geom_errorbar(aes(ymin=Time-sd, ymax=Time+sd), width=.2, position=position_dodge(0.05)) + 
    theme_minimal() + 
    ggtitle('Runtime ~ # Samples') +
    xlab('# Samples') + 
    ylab('Runtime (Seconds)') + 
    geom_point() +
    theme(
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text=element_text(size=15)
    )
                                        
data_ratios <- data.frame(
    "MACS2vMOCHA" = data[data$Model=="MACS2", ]$Time/data[data$Model=="MOCHA", ]$Time,
    "HOMERvMOCHA" = data[data$Model=="Homer", ]$Time/data[data$Model=="MOCHA", ]$Time,
    "Nsamples" = data[data$Model=="MACS2", ]$Nsamples
)
data_ratios <- reshape2::melt(
    data_ratios, 
    measure.vars=c("MACS2vMOCHA", "HOMERvMOCHA"), 
    variable.name="Comparison", value.name="Ratio"
)

pdf('./ss_runtime_parallelv2.pdf')
p2 <- ggplot(
    data = data_ratios,
    aes(x = Nsamples, y = Ratio, col = Comparison, fill = Comparison)
    ) + 
    # scale_x_log10(breaks = trans_breaks("log10", function(x){10^x}), labels = label_number(drop0trailing = TRUE)) +
    geom_bar(stat="identity", position="dodge") + 
    theme_minimal() + 
    scale_color_manual(values=c("#FFD700", "#8B008B")) + 
    scale_fill_manual(values=c("#FFD700", "#8B008B")) +
    theme(
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15),
        axis.text=element_text(size=15)
    ) +
    xlab('Number of Samples')


grid.newpage()
grid.draw(rbind(ggplotGrob(p), ggplotGrob(p2), size = "max"))
dev.off()

# tileResults_nkproliferating <- MOCHA::callOpenTiles(
#   fullCovidArchR,
#   cellPopLabel = "CellSubsets",
#   cellPopulations = c("NK Proliferating"),
#   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#   Org = "org.Hs.eg.db",
#   numCores = 40,
#   studySignal = studySignal,
#   outDir = tempdir()
# )

# tileResults_bright <- MOCHA::callOpenTiles(
#   fullCovidArchR,
#   cellPopLabel = "CellSubsets",
#   cellPopulations = c("NK_CD56bright"),
#   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#   Org = "org.Hs.eg.db",
#   numCores = 40,
#   studySignal = studySignal,
#   outDir = tempdir()
# )

# tileResults_nk <- MOCHA::callOpenTiles(
#   fullCovidArchR,
#   cellPopLabel = "CellSubsets",
#   cellPopulations = c("NK"),
#   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#   Org = "org.Hs.eg.db",
#   numCores = 40,
#   studySignal = studySignal,
#   outDir = tempdir()
# )

# tileResults_cd14 <- MOCHA::callOpenTiles(
#   fullCovidArchR,
#   cellPopLabel = "CellSubsets",
#   cellPopulations = c("CD14 Mono"),
#   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#   Org = "org.Hs.eg.db",
#   numCores = 40,
#   studySignal = studySignal,
#   outDir = tempdir()
# )