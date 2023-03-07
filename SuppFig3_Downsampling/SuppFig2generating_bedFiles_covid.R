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
source('/home/jupyter/scMACS/R/utils.R')

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

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(covidArchR) 

## export cell types 
## as "bulk" 
cellType <- metadf$CellSubsets

## create new 
## sample X cell type column
## to do sample-specific 
## peak-calls 
cellType_sample <- paste(metadf$predictedGroup_Co2, metadf$Sample)
sample_specific=FALSE

# ###########################################################
# ###########################################################

cell_counts <- data.frame(table(cellType))
cell_counts <- cell_counts[order(cell_counts$Freq, decreasing=T),]

# ###########################################################
# ###########################################################

# # Export "bulk" 
# # by poooling cells
# # across samples by cell types

celltypes <- c('DC','NK_CD56Bright',
               'B naive','CD14 Mono', 'CD16 Mono',
              'CD4 CTL TEM')

export_bulk_bedfiles = TRUE
sample_specific = FALSE

for(cell in celltypes){
    pryr::mem_used()
    print(cell)
    cellTypesToExport = cell
    frags <- getPopFrags(covidArchR, metaColumn = "CellSubsets", 
                         cellSubsets = cellTypesToExport, 
                         numCores= 10)
    
    cell_transformed <- gsub(' ','.', cell)
    
    cellType_Frags <- frags
    rm(frags)
    ## create directory for cell type
    cellType_dir <- paste('/home/jupyter/covid/scMACS_manuscript_analyses/ManuscriptData/',
                     cell, 
                     sep='')    
    dir.create(cellType_dir)
    setwd(cellType_dir)
    
    
    ############################################################
    ############################################################
    ## setwd to bulk directories 
    dir.create('bulk')
    setwd('bulk')

    
    ## identify cell barcodes 
    ## from fragment files 
    cellnames <- unique(cellType_Frags[[1]]$RG)
    NCells = length(cellnames)
    
    ### determine quantity
    ### for downsampling 
    cellQuants <- c(rep(5, 10),
                    rep(10, 10),
                    rep(20, 10),
                    rep(30, 10),
                    rep(50, 10),
                    rep(100, 10),
                    rep(200, 10),
                    rep(300, 10),                    
                    rep(500, 10),
                    rep(1000,10),
                    rep(min(2000, NCells),5),
                    rep(min(3000, NCells),5),                    
                    rep(min(5000, NCells),5),
                    rep(min(7500, NCells),5),                    
                    rep(min(10000, NCells),3),
                    rep(min(15000, NCells),3),                    
                    rep(min(20000, NCells),3),
                    rep(min(30000, NCells),3),                    
                    rep(min(40000, NCells),3),                                        
                    rep(min(50000, NCells),3),
                    rep(min(75000, NCells),3),                    
                    rep(min(100000, NCells),3),
                    rep(min(115000, NCells),3),                    
                    rep(NCells,1)
                    )
    
    ### remove repeats (if cell pop)
    ### does not have enough cells 
    ### to go beyond its sample size 
    cellQuants_dt <- data.table::as.data.table(cellQuants)
    cellQuants_dt[,index := 1:.N, by=cellQuants]
    
    ### avoid sampling at 'max count'
    idx <- which.max(cellQuants) 
    cellQuants_dt = cellQuants_dt[1:idx,]
    
    fname <- paste(cellQuants_dt$cellQuants, 'cells_sample',cellQuants_dt$index,sep='_')

    ### use cell quantities to 
    ### sample the specific barcodes
    ### for the downsampling 
    
    barcodes <- lapply(cellQuants_dt$cellQuants, function (x) sample(cellnames, size=x)
           )
    
    ### use cell barcodes 
    ### to subset fragments 
    ### per downsample 
    subset_Frag <- function(cellNames, tmp){
            tmp_df <- GenomicRanges::as.data.frame(tmp)
            idx <- which(tmp_df$RG %in% cellNames)
            tmp[ idx]
           
    }
    
    subsetFragsList <- mclapply(barcodes,
                                function(x)
                                    subset_Frag(x, cellType_Frags[[1]]),
                                mc.cores=5
                                )
    
    normalization_factors <- as.numeric(sapply(subsetFragsList, length))
              
    ############################################################
    ############################################################
    ### call MOCHA peaks 
    
    studySignal = median(metadf$nFrags)
    ### subset archr projects using 
    ### the barcodes to create 
    ### different projects per downsample 
    
    subset_ArchR_projects <- lapply(barcodes, function(x) 
           subsetCells(covidArchR, x)
           )
    pryr::mem_used()

    ### call peaks using 
    ### scMACS functionalities
    peak_lists <- mclapply(1:length(subset_ArchR_projects),
                           function(x){
                               
                           
        tmpfrags <-subsetFragsList[[x]]
        subArchRProj <- subset_ArchR_projects[x][[1]]
                               
        study_prefactor <- 3668 / studySignal
        normalization_factors <- length(tmpfrags)
        
        tilesGRangesList <- MOCHA:::callTilesBySample(
              blackList = ArchR::getBlacklist(subArchRProj),
              returnAllTiles = TRUE,
              totalFrags = normalization_factors,
              fragsList = tmpfrags,
              verbose = T,
              StudypreFactor = study_prefactor
            )
        preds <- as.data.table(tilesGRangesList)

                               
    },
             
                           mc.cores=5)
        
    
    ### create directory
    dir.create('MOCHA')
    setwd('MOCHA/')
  
    ### and export peaks 
    ### to scMACS flat files 
    mclapply(1:length(peak_lists),
             function(x)
                  saveRDS(peak_lists[[x]], 
                          file=paste(fname[x],'_peaks.rds', sep='')
                          ),
             mc.cores=10
             )
                        
    ############################################################
    ############################################################
    if(export_bulk_bedfiles){                    
    ### Export subsets of fragment files 
    ### to generate bed files for 
    ### Homer & Macs2 peak calling 
    
                        
    ### subsetFrag function takes 
    ### the whole fragment list and 
    ### extracts fragments belonging to
    ### user-inputted barcodes 
    subsetFrag <- function(frags, cells ){
        
        tmp_frags  <- frags
        tmp <- frags[[1]]
        tmp <- tmp[tmp$RG %in% cells]
        tmp_frags[[1]] <- tmp
        
        return(tmp_frags)
        }
    
    ### loop through all barcodes 
    ### to create subset fragment
    ### files 
                        
    fragList <- lapply(1:length(barcodes),
                       function(x)
                           subsetFrag(cellType_Frags,barcodes[[x]])
                       )                        
                        
    setwd('../')
    dir.create('downsamples')
    setwd('downsamples')
    
    #Export that fragment info into BigWigs for Insertions and Coverage. 
        mclapply(1:length(fragList),
               function(x)
               
               FragsToCoverageFiles(fragList[[x]], 
                         '/home/jupyter/scATAC_Supplements/hg38_chromInfo.txt', 
                         '/home/jupyter/scATAC_Supplements/hg38.chrom.sizes', 
                         files = 'bedGraph', 
                         numCores= 50,
                         fname=fname[x]
                        ),
                 mc.cores=20
               )
    }

    ############################################################
    ############################################################ 
    rm(fragList,peak_lists,cellType_Frags,subsetFragsList)
    pryr::mem_used()
#    setwd(cellType_dir)
    if(sample_specific){

            ## export bulk datasets
            dir.create('sample_specific')
            setwd('sample_specific')

            ## Export "sample-specific" 
            ## cells by cell types
            idx <- grep(cell, cellType_sample) 
            cellTypes_Samples_ToExport = unique(cellType_sample[idx])

            frags_by_sample <- getPopFrags(covidArchR, metaColumn = "cellType_sample", 
                                 cellSubsets = cellTypes_Samples_ToExport, 
                                 numCores= 50)

            frags_by_sample2 <- frags_by_sample

            for(i in 1:length(frags_by_sample)){
                print(i)
                    print(frags_by_sample[i])
                frags_by_sample2[i] <- frags_by_sample[[i]]
            }

            #Export that fragment info into BigWigs for Insertions and Coverage. 
            FragsToCoverageFiles(frags_by_sample2, 
                                 '/home/jupyter/scATAC_Supplements/hg38_chromInfo.txt', 
                                 '/home/jupyter/scATAC_Supplements/hg38.chrom.sizes', 
                                 files = 'BedGraph',
                                 numCores= 50)
            rm(frags_by_sample, frags_by_sample2)
            pryr::mem_used()
               
    }

}
