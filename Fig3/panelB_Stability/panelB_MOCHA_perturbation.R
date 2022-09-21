# -*- coding: utf-8 -*-
# ###########################################################
# ###########################################################
#
# Author: Samir Rachid Zaim
# Date: 11/06/2021
#
# ###########################################################
# ###########################################################
## Load Libraries
require(scMACS)
### Set directory
setwd('/home/jupyter/covid/scMACS_manuscript_analyses/diff_analysis')
source('/home/jupyter/scMACS/R/get_reproducible_peaks.R')
source('/home/jupyter/scMACS/R/create_sample_peak_matrix.R')
source('/home/jupyter/theme.R')
require(ggpubr)

runTime = F
## Load the ArchR Project
## and extract the metadata
## object
covidArchR <- loadArchRProject("/home/jupyter/FullCovid")
metadf <- getCellColData(covidArchR)
metadf_dt <- as.data.table(getCellColData(covidArchR)) 

sample_metadata = metadf_dt[, list(Visit=first(Visit),
                             Label=first(COVID_status),
                             Days_since_symptoms=days_since_symptoms,                  
                             Sample_CellType=first(new_cellType_sample)), 
                            by=list(Sample, CellSubsets)
                 ]
sample_metadata

cellsubsets <- unique(sample_metadata$CellSubsets)
# ###########################################################
# ###########################################################

## Subset to visit one 
lookup_table <- unique(metadf[,c('Sample',
                                 'COVID_status',
                                 'new_cellType_sample',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

visit_1 <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms),]

## filter out the cd14 monocytes population 
metaFile <- visit_1[grep('CD16 Mono',visit_1$new_cellType_sample),]
metaFile$SampleCellType <- metaFile$new_cellType_sample
metaFile$Class =  metaFile$COVID_status


# ##########################################################
# ##########################################################


# ######################################################################
# ######################################################################

## set parameters for differential analyses 
#cellType_sample <- metadf[,'new_cellType_sample']
cellTypesToExport <- metaFile$SampleCellType 
cellType_Samples= cellTypesToExport
numCores=20
metaColumn='new_cellType_sample'

## call sample specific peaks 
sample_specific_peaks <- callPeaks_by_sample(covidArchR, 
                                metadf,
                                cellType_Samples=cellTypesToExport,
                                metaColumn=metaColumn,         
                                returnAllPeaks=TRUE,
                                numCores=numCores,
                                returnFrags=T
                     
                     )
NCells = sapply(sample_specific_peaks,
                function(x)
                    x[[1]]$numCells[1]
                )

reproduciblePeaks <- get_reproducible_peaks(sample_specific_peaks,
                       metaFile,
                       reproducibility=0.2,
                       fname='cd16s_peaks.png')

sample_peak_matrix <- create_peak_sampleMatrix(sample_specific_peaks,
                         reproduciblePeaks$Union)


###########################################################
###########################################################
if(runTime){
peaks_to_test <- c(3, 229,
                   4210, 20966,
                   55251, 110009,
                   162531,190779,
                   204103)


timeRuns <- lapply(peaks_to_test,
       function(x){
           
        tmp_mat = sample_peak_matrix[ sample(nrow(sample_peak_matrix),x)]
        time <- system.time(
           differential_regions <- get_differential_accessible_regions(
            tmp_mat,
            metaFile, 
            fdr_control=0.2, 
            nCores=50)
            )
        time
           
           }
       )

timeRuns = data.table(do.call(rbind, timeRuns))
timeRuns$PeaksTested = peaks_to_test
timeRuns$RunTime = timeRuns$elapsed
timeRuns$Model = 'scMACS'
timeRuns =  timeRuns[, c('PeaksTested','RunTime','Model')]
write.csv(timeRuns,
          file='/home/jupyter/covid/scMACS_manuscript_analyses/diff_analysis/scMACS_runtime.csv',
          row.names=F
          )
          
    
}

###########################################################
###########################################################
       
differential_regions <- get_differential_accessible_regions(
    sample_peak_matrix,
    metaFile, 
    nCores=50)


###########################################################
###########################################################

###########################################################
###########################################################
### Remove one sample from the 
# Plot sample-specific reproducible peaks
fdr_threshold = 0.2
nCores=50
true_DAPs = differential_regions[FDR <= fdr_threshold]
samples <- metaFile$SampleCellType

remove_one_sample <- function(xx){
    tmp_meta <-  metaFile[metaFile$SampleCellType != xx ,]
    tmp_matrix = sample_peak_matrix[, 
                                    !names(sample_peak_matrix) 
                                    %in% xx, with=F]
    

    ## Get group labels 
    positive_samples <- tmp_meta$SampleCellType[tmp_meta$Class=='Positive']
    group <- ifelse(names(tmp_matrix)[2:ncol(tmp_matrix)] %in% positive_samples,1,0)

    res_pvals = get_differential_accessible_regions(tmp_matrix,
                                        tmp_meta, 55)
    
    # ## find significance 
    sig_pvals = res_pvals[FDR <= fdr_threshold ]
    sig_pvals[!sig_pvals$Peak %in% true_DAPs$Peak] -> tmp
    sig_pvals[sig_pvals$Peak %in% true_DAPs$Peak] -> tmp2
    
    rm(tmp_meta, tmp_matrix)
    pryr::mem_used()
    
    return(res_pvals)
    
}

results = mclapply(samples, function(x){
    print(x)
    remove_one_sample(x)
    },
                   mc.cores=4
                 )

###########################################################
###########################################################

# try 2 different threhsolds
results_by_threshold <- function(tmp, alpha_thresh=0.2, true_DAPs){
    
    curr_res_pvals = tmp
    n_minus1 = curr_res_pvals[FDR <= alpha_thresh]
    
    common = intersect(n_minus1$Peak,true_DAPs$Peak)
    unique = setdiff(n_minus1$Peak,true_DAPs$Peak)
    missed = setdiff(true_DAPs$Peak,n_minus1$Peak)
    
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

###########################################################
###########################################################
# ## Get cumulative 
### new peaks 

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
    fullset = true_DAPs$Peak
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
thresh = fdr_threshold
true_DAPs = differential_regions[FDR <= thresh]
### Threshold 1 
fdr.2 = mclapply(1:length(results),
               function(x) 
                   results_by_threshold(results[[x]], thresh,true_DAPs ),
                 mc.cores=10
               )
fdr.2.vals = rbindlist(lapply(fdr.2, function(x) x$Values
                              )
                      )
fdr.2.vals$CountUnique <- get_new_peaks(fdr.2)
fdr.2.vals$RelativeUnique <- fdr.2.vals$CountUnique / nrow(true_DAPs)
fdr.2.vals$F1 = 2 * (fdr.2.vals$Precision*fdr.2.vals$Recall)/(fdr.2.vals$Precision+fdr.2.vals$Recall)
fdr.2.vals$Conserved_Recall = get_conserved_peaks(fdr.2)
fdr.2.vals$Perc = fdr.2.vals$Conserved_Recall/nrow(true_DAPs)
fdr.2.vals
summary(fdr.2.vals)

row1 = data.frame(
    Common = nrow(true_DAPs),
    Unique = 0,
    missed = 0,
    Recall = 1,
    Precision = 1,
    CountUnique=0,
    RelativeUnique=0,
    F1=1,
    Conserved_Recall=nrow(true_DAPs),
    Perc = 1)

fdr.2.vals = rbind(row1, fdr.2.vals)

setwd('/home/jupyter/covid/scMACS_manuscript_analyses/sample_specific_analyses/dap_cd16')
write.csv(fdr.2.vals,
          file='n-1.daps.scmacs.csv',
         row.names=F)

write.csv(true_DAPs,
          file='cd16_scmacs.csv',
          row.names=F)

###########################################################
###########################################################                  