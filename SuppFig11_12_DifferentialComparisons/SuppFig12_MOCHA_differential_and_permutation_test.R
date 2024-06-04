# -*- coding: utf-8 -*-
# ###########################################################
# ###########################################################
#
# Author: Samir Rachid Zaim
# Date: 11/06/2021
#
# ###########################################################
# ###########################################################

# ###########################################################
# ###########################################################

## Load Libraries
require(MOCHA)
require(ArchR)

## Load the ArchR Project
## and extract the metadata object
ArchRProj <- loadArchRProject("/home/jupyter/FullCovid")
metadata = as.data.table(getCellColData(ArchRProj))
studySignal = median(metadata$nFrags)

# Define your annotation package for TxDb object(s)
# and genome-wide annotation
# Here our samples are human using hg38 as a reference.
# For more info: https://bioconductor.org/packages/3.15/data/annotation/

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

###########################################################
###########################################################

## Get metadata information
## at the sample level
lookup_table <- unique(metadata[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table$Sample[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]

## subset ArchR Project
idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
cellsSample <- ArchRProj$cellNames[idxSample]
ArchRProj <- ArchRProj[cellsSample, ]

###########################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
###########################################################


tileResults <- callOpenTiles( 
    ArchRProj,
    cellPopLabel = "CellSubsets" ,
    cellPopulations = "CD16 Mono",
    TxDb = 'TxDb.Hsapiens.UCSC.hg38.refGene',
    Org = 'org.Hs.eg.db',
    numCores = 20,
    outDir = 'tmp',
    studySignal = studySignal
)

###########################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
###########################################################
groupColumn <- "COVID_status" 
threshold <- 0.2

SampleTileMatrices <- MOCHA::getSampleTileMatrix(
  tileResults,
  cellPopulations = "CD16 Mono",
  groupColumn = groupColumn,
  threshold = threshold,
  verbose = FALSE
)

###########################################################
# 4. Get differential accessibility for specific 
#    cell populations. Here we are comparing CD16  
#    cells between samples where our groupColumn 
#    "COVID_status" is Positive (our foreground) 
#    to Negative samples (our background).
###########################################################

cellPopulation <- "CD16 Mono"
groupColumn <- "COVID_status"
foreground <- "Positive"
background <- "Negative"

# Choose to output a GRanges or data.frame.
# Default is TRUE
outputGRanges <- TRUE

# Optional: Standard output will display the number of tiles found 
# below a false-discovery rate threshold.
# This parameter does not filter results and only affects the 
# afforementioned message. 
fdrToDisplay <- 0.2

differentials <- MOCHA::getDifferentialAccessibleTiles(
  SampleTileObj = SampleTileMatrices,
  cellPopulation = cellPopulation,
  groupColumn = groupColumn,
  foreground = foreground,
  background = background,
  fdrToDisplay = fdrToDisplay,
  outputGRanges = outputGRanges,
  numCores = 20
)

differentials = differentials[!is.na(differentials$FDR)] 
differentials = differentials[differentials$FDR < 0.2]
true_set = differentials
##########################################################
##########################################################

permuted_differentials <- function(seed){
    
    set.seed(seed)
    ### Permute label column
    colData(SampleTileMatrices)$PermutedLabel = colData(SampleTileMatrices)$COVID_status
    colData(SampleTileMatrices)$PermutedLabel = sample(colData(SampleTileMatrices)$COVID_status)


    permuted_differentials <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileObj = SampleTileMatrices,
      cellPopulation = cellPopulation,
      groupColumn = 'PermutedLabel',
      foreground = foreground,
      background = background,
      fdrToDisplay = fdrToDisplay,
      outputGRanges = outputGRanges,
      numCores = 20
    )
    
    permuted_differentials = permuted_differentials[!is.na(permuted_differentials$FDR)]
    ### Get results 
    res = data.table(Iteration =1 ,
               Pvalue = permuted_differentials$P_value,
               FDR = permuted_differentials$FDR)
    
    metadata = colData(SampleTileMatrices)
    metadata = metadata[, c('COVID_status','PermutedLabel','Sample')]
  
    return(list(Results=res,
               metadata=metadata))
}

res1 <- lapply(1:50,
               function(x) permuted_differentials(x)
               )

saveRDS(res1,
        file = 'MOCHA_Revision/Re-run Differential test/permutation.RDS')


sapply(res1,
       function(x) 
                sum(x$Results$FDR < 0.02)
       )


pvalue_reslts <- lapply(1:length(res1),
       function(x){
           y = res1[[x]]
           data.table(
               Pvals = y$Results$Pvalue,
               Iteration = x)}
       )


res = rbindlist(pvalue_reslts)

png('MOCHA_Revision/Re-run Differential test/permutation.png')
ggplot(res,
       aes(x=Pvals))+geom_histogram()+facet_wrap(~Iteration)
dev.off()



###########################################################
###########################################################
validation_differentials= as.data.table(true_set)


# ### Get label column
# metadata = colData(SampleTileMatrices)
# covid_subjects <- metadata$Sample[metadata$COVID_status=='Positive']
# covid_neg_subjects <- metadata$Sample[metadata$COVID_status!='Positive']
# allSamples = c(covid_neg_subjects, covid_subjects)

# ### Generate Sample List 
# sampleList <- 
# lapply(38:25,
#        function(x)
#            lapply(1:15,
#                function(y)
#                    sample(allSamples, x)
#        )
#        )

sampleList_downsampling = readRDS('MOCHA_Revision/Re-run Differential test/sampleList_downsampling.RDS')




########################################################################
########################################################################
### Run analysis on chromosome 4 
differentials = differentials[!is.na(differentials$FDR)] 
differentials = differentials[differentials$FDR < 0.2]
true_set = differentials

chr4= SampleTileMatrices[grep('chr4', rownames(SampleTileMatrices))]

fdrToDisplay <- 0.2

chr4_differentials <- MOCHA::getDifferentialAccessibleTiles(
  SampleTileObj = chr4,
  cellPopulation = cellPopulation,
  groupColumn = groupColumn,
  foreground = foreground,
  background = background,
  fdrToDisplay = fdrToDisplay,
  outputGRanges = outputGRanges,
  numCores = 20
)
verified_dats_chr4 = chr4_differentials[!is.na(chr4_differentials$FDR)] 
verified_dats_chr4 = verified_dats_chr4[verified_dats_chr4$FDR < 0.2]


#### Splitting training and confirmatory cohort 
splitting_differentials <- function(sampleList){

    ### Run 
    discovery_differentials <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileObj = chr4[,SampleTileMatrices$Sample %in% sampleList],
      cellPopulation = cellPopulation,
      groupColumn = 'COVID_status',
      foreground = foreground,
      background = background,
      fdrToDisplay = fdrToDisplay,
      outputGRanges = outputGRanges,
      numCores = 50
    )
       
    discovery_differentials = as.data.table(discovery_differentials)
    predicted_dats= discovery_differentials[FDR < 0.2]
        
    common = length(intersect(predicted_dats$Tile, verified_dats_chr4$Tile))
    recall = common / length(verified_dats_chr4)
    total_predicted = nrow(predicted_dats)
    
    res = data.frame(
        N = length(sampleList),
        Common_Dats = common,
        Precision = common / nrow(predicted_dats),
        FalsePostive = (nrow(predicted_dats) -common) / nrow(predicted_dats),
        Recall=recall)
     
    return(res)
}

### Downsample and assess
### precision, recall, FPR
### And assess recall 
downsample =  lapply(1:14,
                        function(x)
                            lapply(sampleList_downsampling[[x]],
                                function(y)
                                    splitting_differentials(y))
                            )

                 )
                 
system.time(splitting_differentials(seed=1, 38))


total_res <- rbindlist(lapply(downsample, 
           function(x) 
           rbindlist(x)
            
                        )
                       )

       
total_res = rbind(total_res,
                  data.frame(
                      N=39,
                      Common_Dats = 6211,
                      Precision=1,
                      FalsePostive=0,
                      Recall = 1))

png('MOCHA_recall.png')
ggplot(total_res,
       aes(x=as.character(N), 
           y=Recall))+ggrain::geom_rain(scale='width')+
xlab('Total Downsampled')+
ylab('Recall (compared to N=39)')
dev.off()


saveRDS(sampleList,
        file='MOCHA_Revision/Re-run Differential test/sampleList_downsampling.RDS')

write.csv(total_res,
          file= 'MOCHA_Revision/Re-run Differential test/downsamplingRes.csv')