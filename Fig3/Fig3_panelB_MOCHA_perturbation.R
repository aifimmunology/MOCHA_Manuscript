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
require(scMACS)
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
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
    Org = "org.Hs.eg.db",
    numCores = 20,
    studySignal = studySignal
)


###########################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
###########################################################

SampleTileMatrices <- getSampleTileMatrix( 
    tileResults,
    cellPopulations = "CD16 Mono",
    groupColumn = "COVID_status",
    threshold = 0.2
)

###########################################################
# 4. Get differential accessibility for specific 
#    cell populations. Here we are comparing MAIT  
#    cells between samples where our groupColumn 
#    "COVID_status" is Positive (our foreground) 
#    to Negative samples (our background).
###########################################################

differentials <- getDifferentialAccessibleTiles(
    SampleTileObj = SampleTileMatrices,
    cellPopulation = "CD16 Mono",
    groupColumn = "COVID_status",
    foreground =  "Positive",
    background =  "Negative",
    fdrToDisplay = 0.2,
    outputGRanges = FALSE,
    numCores = 50
)

# ##########################################################
# ##########################################################

###########################################################
###########################################################
### Remove one sample from the 
# Plot sample-specific reproducible peaks

fdr_threshold = 0.2
nCores=50
true_DAPs = differentials[FDR <= fdr_threshold]
samples = colnames(SampleTileMatrices)

remove_one_sample <- function(xx){
    
    ### Remove one sample
    ### from analysis 
    tmp_STM = SampleTileMatrices[, !colnames(SampleTileMatrices) %in% xx]

    ### re-run differentials 
    tmp_differentials <- getDifferentialAccessibleTiles(
        SampleTileObj = tmp_STM,
        cellPopulation = "CD16 Mono",
        groupColumn = "COVID_status",
        foreground =  "Positive",
        background =  "Negative",
        fdrToDisplay = 0.2,
        outputGRanges = FALSE,
        numCores = 50
    )

    tmp_differentials= tmp_differentials[FDR <= fdr_threshold]
    ### return final results 
    pryr::mem_used()
    return(tmp_differentials)

}

results = lapply(samples, 
                   function(x){
                        print(x)
                        remove_one_sample(x)
                        }
                 )

# ##########################################################
# ##########################################################

overlap_results<- function(tmp, true_DAPs){
    
    n_minus1 = tmp
    
    common = intersect(n_minus1$Tile, true_DAPs$Tile)
    unique = setdiff(n_minus1$Tile,   true_DAPs$Tile)
    missed = setdiff(true_DAPs$Tile,  n_minus1$Tile)
    
    res = data.frame(Common = length(common),
                     Unique = length(unique),
                     missed = length(missed)                     
                    )
    res$Recall = res$Common / nrow(true_DAPs)
    res$Precision = res$Common / nrow(n_minus1)  
    
    res2 = list(Common=common,
                Unique=unique,
                Missed=missed)
    
    
    return(list(Values=res,
                Regions=res2)
           )
}



# ##########################################################
# ##########################################################
# ## Get cumulative 
# ## new peaks 

get_new_peaks <- function(fdr.2){
     
    unique_cumulative=lapply(fdr.2, function(x) x$Regions$Unique
               )

    cumulative_new_peaks <- character()
    count = numeric(39)
    for(i in 1:39){
    
        curr_new_peaks <- unique_cumulative[[i]]
        cumulative_new_peaks = unique(c(cumulative_new_peaks, curr_new_peaks))
        count[i] = length(cumulative_new_peaks)  
                                    
    }
    return(count) 

} 

get_conserved_peaks <- function(fdr.2){
     
    unique_cumulative=lapply(fdr.2, function(x) x$Regions$Common
               )
    fullset = true_DAPs$Tile
    count = numeric(0)
    for(i in 1:39){
    
        curr_conserved <- unique_cumulative[[i]]
        fullset = intersect(fullset, curr_conserved)
        count = c(count, length(fullset))
                                    
    }
    return(count) 

}    


###########################################################
###########################################################

### Threshold 1 
results2 = lapply(results, 
                   function(x){
                        overlap_results(x, true_DAPs)
                        }
                 )

results3 = rbindlist(lapply(results2, function(x) x$Values
                              )
                      )

results3$CountUnique <- get_new_peaks(results2)
results3$RelativeUnique <- results3$CountUnique / nrow(true_DAPs)
results3$F1 = 2 * (results3$Precision*results3$Recall)/(results3$Precision+results3$Recall)

results3$Conserved_Recall = get_conserved_peaks(results2)
results3$Perc = fdr.2.vals$Conserved_Recall/nrow(true_DAPs)
results3
summary(results3)

row1 = data.frame(
    Common = nrow(true_DAPs),
    Unique = 0,
    missed = 0,
    Recall = 1,
    Precision = 1,
    CountUnique=0,
    RelativeUnique=0,
    F1=1,
    Conserved_Recall=nrow(true_DAPs)
)

results3 = rbind(row1, results3)

setwd('/home/jupyter/MOCHA_Manuscript/Fig3/panelB_Stability/')
write.csv(results3,
          file='n-1.daps.MOCHA.csv',
         row.names=F)

# ##########################################################
# ##########################################################                  
