setwd("MOCHA_runtime_comparison/")
install.packages("ggbio")
install.packages("tictoc")
devtools::install("../MOCHA")

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
n
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

source("./runtimeHelperFunctions_by_ncells.R")
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
cellQuants <- c(
    min(50000, TotalNCells),
    min(100000, TotalNCells)
)
cellQuants <- c(10, 50, 100, 250, 500, 1000)
cellQuants <- c(
  min(2000, TotalNCells),
  min(3000, TotalNCells),
  min(5000, TotalNCells),
  min(10000, TotalNCells),
  min(15000, TotalNCells)
)
cellQuants <- c(
  min(50000, TotalNCells),
  min(100000, TotalNCells)
)
# Remove repeats (if cell pop) does not have enough cells
# to go beyond its sample size
idxMax <- which.max(cellQuants)
cellQuants <- cellQuants[1:idxMax]
name <- "20k-100k_rep2"
############################################################
#      HOMER
############################################################
homerOutDir <- "./homer_runtimePeaks"

homer_runtimes.df <- runBenchmark(
  peakCaller = homer_measure_runtimes,
  outDir = homerOutDir,
  numRepeats = 10,
  name = name
)
pbPost("note", "Homer complete, check output homer_runtimes.df")

############################################################
#      MACS2
############################################################
macs2OutDir <- "./macs2_runtimePeaks"

Macs2_runtimes.df <- runBenchmark(
  peakCaller = macs2_measure_runtimes,
  outDir = macs2OutDir,
  numRepeats = 10,
  name = name
)
pbPost("note", "MACS2 complete, check output Macs2_runtimes.df")


############################################################
#      MOCHA
############################################################
mochaOutFile <- "mocha_runtimePeaks"

MOCHA_runtimes.df <- runBenchmark(
  peakCaller = MOCHA_measure_runtime,
  outDir = mochaOutFile,
  numRepeats = 5,
  name = name
)
pbPost("note", "MOCHA complete, check output mocha_runtimes.df")

############################################################
# Plotting
############################################################
Macs2_runtimes.df <- read.csv("./results/macs2_runtimePeaksresults.csv")
MOCHA_runtimes.df <- read.csv("./results/mocha_runtimePeaksresults.csv")
homer_runtimes.df <- read.csv("./results/homer_runtimePeaksresults.csv")


MOCHA_runtimes.df = as.data.table(read_xlsx('./MOCHA_runtime_comparison.xlsx',sheet = 1))
Macs2_runtimes.df = as.data.table(read_xlsx('./MOCHA_runtime_comparison.xlsx',sheet = 2))
homer_runtimes.df = as.data.table(read_xlsx('./MOCHA_runtime_comparison.xlsx',sheet = 3))

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

data <- data_summary(combined_res, varname="Time", groupnames=c("Model", "Ncells"))
  
pdf('/home/jupyter/MOCHA_runtime_comparison/runtime_minimal_v3.pdf')
p<- ggplot(
    data=data,
    aes( x=Ncells, y=Time, col=Model, fill=Model)
    ) + 
    geom_line() + 
    theme_minimal() + 
    ggtitle('Runtime ~ # Cells') +
    xlab('# Cells') + 
    ylab('Run time (Seconds)') + 
    geom_point() +
    geom_errorbar(aes(ymin=Time-sd, ymax=Time+sd), width=.2, position=position_dodge(0.05)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x){10^x}), labels = trans_format("log10", math_format(10^.x))) +
    theme(
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text=element_text(size=15)
    )
                                        
data_ratios <- data.frame(
    "MACS2vMOCHA" = data[data$Model=="MACS2", ]$Time/data[data$Model=="MOCHA", ]$Time,
    "HOMERvMOCHA" = data[data$Model=="Homer", ]$Time/data[data$Model=="MOCHA", ]$Time,
    "Ncells" = data[data$Model=="MACS2", ]$Ncells
)
data_ratios <- reshape2::melt(
    data_ratios, 
    measure.vars=c("MACS2vMOCHA", "HOMERvMOCHA"), 
    variable.name="Comparison", value.name="Ratio"
)

pdf('/home/jupyter/MOCHA_runtime_comparison/runtime_minimal_v3.pdf')
p2 <- ggplot(
    data = data_ratios,
    aes(x = Ncells, y = Ratio, col = Comparison, fill = Comparison)
    ) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x){10^x}), labels = trans_format("log10", math_format(10^.x))) +
    geom_bar(stat="identity", position="dodge") + 
    theme_minimal() + 
    scale_color_manual(values=c("#FFD700", "#8B008B")) + 
    scale_fill_manual(values=c("#FFD700", "#8B008B")) +
    theme(
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15),
        axis.text=element_text(size=15)
    ) +
    xlab('# Cells')


grid.newpage()
grid.draw(rbind(ggplotGrob(p), ggplotGrob(p2), size = "first"))
dev.off()
