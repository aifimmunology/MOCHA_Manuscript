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
require(ArchR)
source('/home/jupyter/covid/scMACS_manuscript_analyses/helper_granges.R')

# Load the ArchR Project
lp_ArchRProject <- loadArchRProject("/home/jupyter/longPilot")

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(lp_ArchRProject) 

## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup_Col2.5'
sample_celltypeName <- 'sample_celltype'
cellType <- metadf[,cellColName]

metadf_dt <- as.data.table(getCellColData(lp_ArchRProject)) 
### 

numCores=10
metaColumn='new_sample_cellType'

celltypes = c('B naive','CD16 Mono', 'CD8 TEM')

## boolean whether or not 
## to call peaks 
callPeaks = TRUE

for(cell in celltypes){
    print(paste('analyzing', cell))

            setwd('/home/jupyter/longPilotValidation/sample_specific_analyses/data')
        setwd(cell)

        ############################################################
        ############################################################
        setwd('rawFiles')
        ## change directory
        ## to location where the 
        ## data exist 
    
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
        dir.create('../sample_specific_macs2_peaks')   
            
        cmd_sample<- lapply(1:length(samples),

           function(x)
                  paste('macs2 callpeak -t', fnames[x],
                
           '-g hs -f BED --nolambda --shift -75 --extsize 150 --broad',
                         '--nomodel -n',

                    paste('../sample_specific_macs2_peaks/',samples[x],sep='')
                       )
                                )

        mclapply(cmd_sample,
                 function(x) system(x),
                 mc.cores=30
                 )
        
        
        ############################################################
        ############################################################

        ## Homer Peaks 
        dir.create('../sample_specific_homer_peaks')         
        dir.create('../sample_specific_homer_peaks/dirs')
        outdir <- paste('sample_specific_homer_peaks/dirs/', samples, sep='')
    
        fnames = paste('rawFiles', fnames, sep='/')

        setwd('..')
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
        setwd('sample_specific_homer_peaks/dirs')
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
