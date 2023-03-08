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
rm(getPopFrags)

require(data.table)
require(ggplot2)
require(ggpubr)
require(scMACS)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)

# Load the ArchR Project
lp_ArchRProject <- loadArchRProject("/home/jupyter/longPilot")

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(lp_ArchRProject) 

## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup_Col2.5'
sample_celltypeName <- 'sample_celltype'
cellType <- metadf[,cellColName]


#######################################################################
metadf_dt <- as.data.table(getCellColData(lp_ArchRProject)) 


numCores=10
metaColumn='new_sample_cellType'

celltypes = c('B naive','CD16 Mono', 'CD8 TEM')
homeDir = '/home/jupyter/longPilotValidation/sample_specific_analyses/'
setwd(homeDir)
dir.create('results'); dir.create('data')
cellTypes_sample <- metadf_dt$new_sample_cellType
callMacs2Homer=T
callscMACS=T


extract_info_by_cell <- function(cell){
    
    cellTypesToExport = unique(cellTypes_sample[grep(paste(cell,'$',sep=''), cellTypes_sample)])
    frags <- scMACS::getPopFrags(lp_ArchRProject, 
                                 metaColumn = "new_sample_cellType", 
                         cellSubsets = cellTypesToExport, 
                         numCores= 30)
    
    cell_transformed <- gsub(' ','.', cell)
    
    cellType_Frags <- frags
    rm(frags)
    ## create directory for cell type
    cellType_dir <- paste('/home/jupyter/longPilotValidation/sample_specific_analyses/data/',
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

    if(callscMACS){
        ## call scMACS peaks 
        ## call sample specific peaks 
        sample_specific_peaks <- callPeaks_by_sample(lp_ArchRProject, 
                                    metadf,
                                    cellType_Samples=cellTypesToExport,
                                    metaColumn=metaColumn,         
                                    returnAllPeaks=TRUE,
                                    numCores=numCores,
                                    returnFrags=T

                         )

        ### filter any possible errors 
        fnames <- sapply(sample_specific_peaks, function(x) names(x))
        ### and export peaks 
        ### to scMACS flat files 
        mclapply(1:length(sample_specific_peaks),
                 function(x)
                      saveRDS(sample_specific_peaks[[x]], 
                              file=paste(fnames[x],'_peaks.rds', sep='')
                              ),
                 mc.cores=20
                 )
        rm( sample_specific_peaks,peak_lists)
        pryr::mem_used()

    }
    ############################################################
    ############################################################
    if(callMacs2Homer){
        
        ## identify cell population to export 
        ## across samples by cell types
        setwd('..')
   
        cell_transformed <- gsub(' ','.', cell)
    
        #cellType_Frags <- frags
        #rm(frags)
        ## setwd to bulk directories 
        dir.create('rawFiles')
        setwd('rawFiles')
        ## Export "sample-specific" 
        ## cells by cell types
        frags_by_sample2 <- cellType_Frags


        for(i in 1:length(cellType_Frags)){
                frags_by_sample2[i] <- cellType_Frags[i]
        }
        names(frags_by_sample2) <- gsub('_','#', 
                                        names(frags_by_sample2))

        #Export that fragment info into BigWigs for Insertions and Coverage. 
        FragsToCoverageFiles(frags_by_sample2, 
                             '/home/jupyter/scATAC_Supplements/hg38_chromInfo.txt', 
                             '/home/jupyter/scATAC_Supplements/hg38.chrom.sizes', 
                             files = 'BedGraph',
                             numCores= 50)
        rm(cellType_Frags, frags_by_sample2)
        pryr::mem_used()
        
    }

}

tictoc::tic()                     
lapply(celltypes,
       function(x)
           extract_info_by_cell(x)
       )
tictoc::toc()                     