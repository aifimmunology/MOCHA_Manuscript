# ###########################################################
# ###########################################################

# Author: Samir Rachid Zaim
# Date: 11/06/2021

# Desc: 
#     This script calls peaks using
#     macs2, homer, and hmmratac 
#     to create those comparisons 
#     with scMACS.

# ###########################################################
# ###########################################################

require(parallel)
require(GenomicRanges)
require(plyranges)
require(data.table)
source('/home/jupyter/covid/scMACS_manuscript_analyses/helper_granges.R')
bm_ArchRProject <- ArchR::loadArchRProject('/home/jupyter/DoubletFreeBoneMarrow/')
require(ArchR)

## boolean whether or not 
## to call peaks 
callPeaks = T

## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup_Col2.5'
sample_celltypeName <- 'sample_celltype'

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(bm_ArchRProject) 

## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup'
metadf$cellType <- metadf[,cellColName]
cellType = metadf$cellType <- metadf[,cellColName]
sample_specific=FALSE


# ###########################################################
# ## call macs2 peaks 
# without dynamic lambda

celltypes <- c('18_Plasma','20_CD4.N1',
                '25_NK','10_cDC',
               '17_B','12_CD14.Mono.2'
               )

cell = celltypes[1]

for(cell in celltypes){
    print(paste('analyzing', cell))

        setwd('/home/jupyter/BoneMarrowValidation/data')
        setwd(cell)

        ############################################################
        ############################################################

        ## change directory
        ## to location where the 
        ## data exist 
    
        setwd('bulk/downsamples')

        fnames <- dir()
        fnames <- fnames[grep('bedGraph', fnames)]
        
        archLog_idx <- grep('ArchR', fnames)
        if(length(archLog_idx) >0){
            fnames <- fnames[- grep('ArchR', fnames)]
            
        }
        samples <- gsub('_perThousandCells.bedGraph', '', fnames)

        ############################################################
        ############################################################

        ## MACS2 Peaks 
        dir.create('../downsample_macs2_peaks')   
            
        cmd_sample<- lapply(1:length(samples),

           function(x)
                  paste('macs2 callpeak -t', fnames[x],
                
           '-g hs -f BED --nolambda --shift -75 --extsize 150 --broad',
                         '--nomodel -n',

                    paste('../downsample_macs2_peaks/',samples[x],sep='')
                       )
                                )

        mclapply(cmd_sample,
                 function(x) system(x),
                 mc.cores=30
                 )
       
        ############################################################
        ############################################################

        ## Homer Peaks 
        dir.create('../downsample_homer_peaks')   
        
        dir.create('../downsample_homer_peaks/homer_dirs')
        outdir <- paste('../downsample_homer_peaks/homer_dirs/', samples, sep='')

        ## create command list
        cmd_list <- mclapply(1:length(fnames), 
             function(i)
                     paste('makeTagDirectory',
                           outdir[i],fnames[i])
                         )
        
                
        mclapply(cmd_list,
                 function(x)
                     system(x),
                 mc.cores=20
                 )

        ## call peaks 
        setwd('../downsample_homer_peaks/homer_dirs')
        tagDirs <- dir()

        peak_file_names <- paste(tagDirs,'/peaks.txt', sep='')

        non_tag_dirs <- grep('NA|Rplots|ArchR|union',peak_file_names)

        if(length(peak_file_names[-non_tag_dirs]) >0 ){
            peak_file_names <- peak_file_names[-non_tag_dirs]
            tagDirs <- tagDirs[-non_tag_dirs]
        }

        ## make broad peak calls
        ## with "-style histone" 
        ## to identify enriched 
        ## regions rather than
        ## peaks 
        callPeaks_cmds <- mclapply(1:length(tagDirs),
                                   function(x)
                                      paste('findPeaks',tagDirs[x],'>',
                                            peak_file_names[x],
                                            '-style histone',
                                            sep=' '
                                            )
                                   )
        
        
        mclapply(callPeaks_cmds,
             function(ZZ)
                 system(ZZ),
             mc.cores=20
             )   

}
