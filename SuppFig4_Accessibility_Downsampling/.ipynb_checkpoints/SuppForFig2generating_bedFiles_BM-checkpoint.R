############################################################
############################################################

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
#source('/home/jupyter/getPopFrags_v3.R')

require(data.table)
require(ggplot2)
require(ggpubr)
require(scMACS)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)
# Load the ArchR Project
boneMarrow_ArchRProject <- loadArchRProject("/home/jupyter/DoubletFreeBoneMarrow")

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(boneMarrow_ArchRProject) 
metadf$cellName = row.names(metadf)
## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup'
metadf$cellType <- metadf[,cellColName]
cellType = metadf$cellType <- metadf[,cellColName]
metadf$sample_celltype <- paste(metadf$Sample,metadf$cellType,sep='_')

## create new 
## sample X cell type column
## to do sample-specific 
## peak-calls 
cellType_sample <- metadf$sample_celltypeName
sample_specific=FALSE


### get fragment counts
metadf_dt = as.data.table(metadf)
fragments_counts <-  metadf_dt[, list(Count=.N, 
                                      nFrags=median(nFrags)),
                               by=new_cellType]
############################################################
############################################################

## Add these new cell groups as a 
## metadata column in the ArchR project

# boneMarrow_ArchRProject <- addCellColData(boneMarrow_ArchRProject, 
#                              data =cellType, 
#                              name = "new_cellType", 
#                              cells = getCellNames(boneMarrow_ArchRProject),
#                                           force=TRUE
#                                          )
           
# saveArchRProject(boneMarrow_ArchRProject)
## Add genome hg19
## Export using hg19 chrome size
## export using hg19 chrome info 
        
#############################################
#############################################           
            
############################################################
############################################################

## Export "bulk" 
## by poooling cells
## across samples by cell types

celltypes <- c('20_CD4.N1',
                '12_CD14.Mono.2','10_cDC'
               )


callPeaks = T 
cell=celltypes[1]

setwd('/home/jupyter/BoneMarrowValidation')
dir.create('ManuscriptData')
setwd('ManuscriptData')

for(cell in celltypes){
    print(cell)
    cellTypesToExport = unique(cellType[cellType == cell])
    
    frags <- scMACS::getPopFrags(boneMarrow_ArchRProject, 
                                 metaColumn = "predictedGroup", 
                         cellSubsets = cellTypesToExport, 
                         numCores= 30)
    
    cell_transformed <- gsub(' ','.', cell)
    
    cellType_Frags <- frags
    rm(frags)
    ## create directory for cell type
    cellType_dir <- paste('/home/jupyter/BoneMarrowValidation/ManuscriptData/',
                     cell, 
                     sep='')    
    dir.create(cellType_dir)
    setwd(cellType_dir)
    
    
    ############################################################
    dir.create('bulk')
    setwd('bulk')

    ## identify cell barcodes 
    ## from fragment files 
    cellnames <- unique(cellType_Frags[[1]]$RG)
    NCells = length(cellnames)
    
    ### determine quantity
    ### for downsampling 
    cellQuants <- c(rep(5, 10),
                    rep(10,10),
                    rep(25, 10),
                    rep(50, 10),
                    rep(100, 10),
                    rep(250, 10),
                    rep(500, 10),
                    rep(1000,5),
                    rep(min(2000, NCells),5),
                    rep(min(3000, NCells),5),                    
                    rep(min(5000, NCells),5),
                    rep(min(10000, NCells),3),
                    rep(min(15000, NCells),3),
                    rep(NCells,5)
                    )
    
    ### remove repeats (if cell pop)
    ### does not have enough cells 
    ### to go beyond its sample size 
    idx <- which.max(cellQuants)
    cellQuants <- c(cellQuants[cellQuants < NCells], NCells)
    cellQuants
    
    
    cellQuants_dt <- data.table::as.data.table(cellQuants)
    cellQuants_dt[,index := 1:.N, by=cellQuants]
    fname <- paste(cellQuants_dt$cellQuants, 'cells_sample',cellQuants_dt$index,sep='_')

    ### use cell quantities to 
    ### sample the specific barcodes
    ### for the downsampling 
    
    barcodes <- lapply(cellQuants, function (x) sample(cellnames, size=x)
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
    
              
    ############################################################
    ############################################################
    ### call MOCHA peaks 
    
    studySignal = median(metadf$nFrags)
    ### subset archr projects using 
    ### the barcodes to create 
    ### different projects per downsample 
    
    subset_ArchR_projects <- lapply(barcodes, function(x) 
           subsetCells(boneMarrow_ArchRProject, x)
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
                   
    if(callPeaks){
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
                         '/home/jupyter/scATAC_Supplements/hg19_chromInfo.txt', 
                         '/home/jupyter/scATAC_Supplements/hg19.chrom.sizes', 
                         files = 'BedGraph', 
                         numCores= 50,
                         fname=fname[x]
                        ),
                 mc.cores=30
               )

    ############################################################
    ############################################################ 
    rm(fragList,peak_lists,cellType_Frags)
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
                                 '/home/jupyter/scATAC_Supplements/hg19_chromInfo.txt', 
                                 '/home/jupyter/scATAC_Supplements/hg19.chrom.sizes', 
                                 files = 'BedGraph',
                                 numCores= 50)
            rm(frags_by_sample, frags_by_sample2)
            pryr::mem_used()
               
    }
 }   
}