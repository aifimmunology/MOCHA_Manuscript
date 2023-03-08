# ###########################################################
# ###########################################################

# Author: Samir Rachid Zaim
# Date: 11/06/2021

# Desc: 
#     This script generates the 
#     bedfiles used as inputs for
#     macs2, homer, and hmmratac 
#     to create those comparisons 
#     with scMACS.

############################################################
############################################################
setwd('/home/jupyter/scATAC_Supplements')
source('ArchR_Export_ATAC_Data.R')
source('/home/jupyter/getPopFrags_v3.R')

require(data.table)
require(ggplot2)
require(ggpubr)
require(scMACS)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)

# Load the ArchR Project
covidArchR <- loadArchRProject("/home/jupyter/FullCovid")
metadf <- getCellColData(covidArchR) 


#######################################################################
metadf_dt <- as.data.table(getCellColData(covidArchR)) 
###########################################################
###########################################################

## Subset to visit one 
lookup_table <- unique(metadf[,c('Sample',
                                 'COVID_status',
                                 'new_cellType_sample',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

visit_1 <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms),]

callMacs2Homer = TRUE

###########################################################
###########################################################
# Export "bulk" 
# by poooling cells
# across samples by cell types

numCores=10
metaColumn='new_cellType_sample'

celltypes =unique(gsub('[A-Z0-9_]*-[A-Z0-9]*__','', 
                       lookup_table$new_cellType_sample))
homeDir = '/home/jupyter/covid/scMACS_manuscript_analyses/sample_specific_analyses/'

get_bedgraphs_and_calls <- function(cell){
    print(cell)
    
    ## filter out the cd14 monocytes population 
    metaFile <- visit_1[grep(cell,visit_1$new_cellType_sample),]
    metaFile$SampleCellType <- metaFile$new_cellType_sample
    metaFile$Class =  metaFile$COVID_status
    metaColumn='new_cellType_sample'
    cellTypesToExport <- metaFile$SampleCellType 
 
    ## create directory for cell type
    cellType_dir <- paste(homeDir,'/data/',
                     cell, 
                     sep='')    
    dir.create(cellType_dir)
    setwd(cellType_dir)
    ############################################################
    ############################################################

    ############################################################
    ############################################################
    dir.create('scMACS')
    setwd('scMACS')

    ## call scMACS peaks 
    ## call sample specific peaks 
    sample_specific_peaks <- callPeaks_by_sample(covidArchR, 
                                metadf,
                                cellType_Samples=cellTypesToExport,
                                metaColumn=metaColumn,         
                                returnAllPeaks=TRUE,
                                numCores=numCores,
                                returnFrags=T
                     
                     )
    fnames <- sapply(sample_specific_peaks, function(x) names(x))
    ### and export peaks 
    ### to scMACS flat files 
    mclapply(1:length(sample_specific_peaks),
             function(x)
                  saveRDS(sample_specific_peaks[[x]], 
                          file=paste(fnames[x],'_peaks.rds', sep='')
                          ),
             mc.cores=numCores
             )
    rm( sample_specific_peaks)
    pryr::mem_used()
    
    ############################################################
    ############################################################
    if(callMacs2Homer){
        
        ## identify cell population to export 
        ## across samples by cell types
        setwd('..')
        cellTypesToExport <- unique(cellTypesToExport)
        frags <- getPopFrags(covidArchR, 
                         metaColumn = metaColumn, 
                         cellSubsets = cellTypesToExport, 
                         numCores= numCores)
    
        fnames <- names(frags)
        cell_transformed <- gsub(' ','.', cell)
    
        cellType_Frags <- frags
        rm(frags)
        ## setwd to bulk directories 
        dir.create('rawFiles')
        setwd('rawFiles')
        ## Export "sample-specific" 
        ## cells by cell types
        frags_by_sample2 <- cellType_Frags

        for(i in 1:length(cellType_Frags)){
                frags_by_sample2[i] <- cellType_Frags[i]
        }

        ### Export fragment
        ### files to bedgraph
        FragsToCoverageFiles(frags_by_sample2,
                             '/home/jupyter/scATAC_Supplements/hg38_chromInfo.txt', 
                             '/home/jupyter/scATAC_Supplements/hg38.chrom.sizes', 
                             files = 'BedGraph',
                             numCores= 50)
        rm(cellType_Frags, frags_by_sample2)
        pryr::mem_used()
        
    }

}
lapply(celltypes[2:length(celltypes)], function(x)
    get_bedgraphs_and_calls(x)
       )