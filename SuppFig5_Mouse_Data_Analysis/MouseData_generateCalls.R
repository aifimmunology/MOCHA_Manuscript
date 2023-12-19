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
setwd('/home/jupyter/MOCHA_Revision/')

source('/home/jupyter/scATAC_Supplements/ArchR_Export_ATAC_Data.R')
require(data.table)
require(ggplot2)
require(ggpubr)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)
require(parallel)

### Load the ArchR Project
mouseArchR <- loadArchRProject("/home/jupyter/MOCHA_Revision/MouseData")

######################################################################
####### get TSS for mouse
#bioconductor::install('https://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm9.knownGene.html')

######################################################################
######################################################################

### removal doublet
#mouseArchR = addDoubletScores(mouseArchR)
#mouseArchR = filterDoublets(mouseArchR)
#saveArchRProject(mouseArchR)

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(mouseArchR) 
metadf$CellID = row.names(metadf)
extended_metadata = fread('MouseDataAnalysis/cell_metadata.tissue_freq_filtered.txt')

### Create Cell ID 
extended_metadata$CellID = paste(extended_metadata$tissue.replicate,
                                          extended_metadata$cell, sep='#')

### join different metadata
metadf = dplyr::left_join(as.data.frame(metadf), extended_metadata[,c('cell_label','CellID','tissue.replicate')], multiple='first')

## export cell types 
## as "bulk" 
cellType <- metadf$Sample

## create cell type X organ analysis
metadf$CellType_Organ = paste(metadf$cell_label,'#',
                              metadf$tissue.replicate, sep='')

### replace empty space with '_'
metadf$CellType_Organ = gsub(' ', '_',metadf$CellType_Organ)

### add columns back to ArchR
mouseArchR = addCellColData(mouseArchR, data = metadf$cell_label, name= 'CellPop', cells = metadf$CellID, force=T)
mouseArchR = addCellColData(mouseArchR, data = metadf$CellType_Organ, name= 'CellType_Organ', cells = metadf$CellID, force=T)

saveArchRProject(mouseArchR)

### calculate counts 
metadt = as.data.table(metadf)
tmp = metadt[, list(nFrags=median(nFrags), 
                    cellCount = .N),
             by=list(cell_label,tissue.replicate)]

tmp$CellOrgan = paste(tmp$cell_label, 
                      tmp$tissue.replicate,
                      sep='#')
                      

tmp = tmp %>% arrange(nFrags)

tmp2 = metadt[, list(nFrags=median(nFrags), 
                    cellCount = .N),
             by=list(tissue.replicate)]
###########################################################
###########################################################

#### Choose samples 
#astrocytes = prefrontal
#enterocytes = gut
#granule = cerebellum
#cardio = heart 
#collecting duct / tubule  = kidney 
#   = liver 
#'Type I pneumocytes' = Lung

cellList = c('Astrocytes#PreFrontalCortex_62216',
             'Cardiomyocytes#HeartA_62816', 
             'Enterocytes#LargeIntestineA_62816',
             'Cerebellar_granule_cells#Cerebellum_62216',
             'Collecting_duct#Kidney_62016', 
             'Endothelial_I_cells#Liver_62016',
             'Podocytes#Kidney_62016',
             'Sperm#Testes_62016',
             'Enterocytes#SmallIntestine_62816',
             'T_cells#BoneMarrow_62216', 
             'NK_cells#BoneMarrow_62216', 
             'Microglia#WholeBrainA_62216',
             'Purkinje_cells#WholeBrainA_62216',
             'Hematopoietic_progenitors#BoneMarrow_62216',
             'Regulatory_T_cells#BoneMarrow_62216',
             'Proximal_tubule#Kidney_62016',
             'Alveolar_macrophages#Lung1_62216',
             'Alveolar_macrophages#Lung2_62216',             
             'Type_I_pneumocytes#Lung2_62216',
             'Type_I_pneumocytes#Lung1_62216')
            
cell_counts <- data.frame(table(cellType))
cell_counts <- cell_counts[order(cell_counts$Freq, decreasing=T),]
tmp$CellOrgan = gsub(' ','_',tmp$CellOrgan)

tmp[CellOrgan %in% cellList]

# ###########################################################
# ###########################################################

# # Export "bulk" 
# # by poooling cells
# # across samples by cell types

export_bulk_bedfiles = FALSE
sample_specific = FALSE

for(cell in cellList){
    pryr::mem_used()
    print(cell)
    cellTypesToExport = cell
    frags <- MOCHA::getPopFrags(mouseArchR, cellPopLabel = "CellType_Organ", 
                         cellSubsets = cell, 
                         numCores= 30)
    
    cell_transformed <- gsub(' ','.', cell)
    
    cellType_Frags <- frags
    rm(frags)
    ## create directory for cell type
    cellType_dir <- paste('/home/jupyter/MOCHA_Revision/MouseDataAnalysis/Data/',
                     cell, 
                     sep='')    
    dir.create(cellType_dir)
    setwd(cellType_dir)
    
    ############################################################
    ############################################################
    ## setwd to bulk directories 
    dir.create('bulk')
    setwd('bulk')

    normalization_factors <- length(cellType_Frags[[1]])
    
    ### Calculate study signal
    ### by different organs
    ### to reflect different biology
    organ = unlist(strsplit(cell, '#'))[2]
    studySignal_median = median(metadt[CellType_Organ ==cell]$nFrags, na.rm=T)
    studySignal_mean = mean(metadt[CellType_Organ ==cell]$nFrags, na.rm=T)
    
    combinedSignal = mean(c(studySignal_median, studySignal_mean))

    subset_ArchR_projects = mouseArchR
    
    ### call peaks using 
    ### scMACS functionalities
    tmpfrags <-cellType_Frags
    subArchRProj <- subset_ArchR_projects
    
    ### calculate adjustment prefactor
    study_prefactor <- 3668 / combinedSignal
        
    ### call open tiles 
    tilesGRangesList <- MOCHA:::callTilesBySample(
              blackList = ArchR::getBlacklist(subArchRProj),
              returnAllTiles = TRUE,
              totalFrags = normalization_factors,
              fragsList = tmpfrags[[1]],
              verbose = T,
              StudypreFactor = study_prefactor
            )
    preds <- as.data.table(tilesGRangesList)
    
    ### create directory
    dir.create('MOCHA')
    setwd('MOCHA/')
  
    ### and export peaks 
    ### to scMACS flat files 
    saveRDS(preds, 
            file=paste(cell,'_peaks.rds', sep='')
                          )
                
    ############################################################
    ############################################################
    if(export_bulk_bedfiles){                    
    ### Export subsets of fragment files 
    ### to generate bed files for 
    ### Homer & Macs2 peak calling 
                                               
    fragList <- cellType_Frags           
                        
    setwd('../')
    dir.create('bedfiles')
    setwd('bedfiles')
    
    #Export that fragment info into BigWigs for Insertions and Coverage. 
    FragsToCoverageFiles(fragList, 
                         '/home/jupyter/MOCHA_Revision/MouseDataAnalysis/mm9_chromInfo.txt', 
                         '/home/jupyter/MOCHA_Revision/MouseDataAnalysis/mm9.chrom.sizes', 
                         files = 'bedGraph', 
                         numCores= 50,
                         fname=cell
                        )

 } 
}
