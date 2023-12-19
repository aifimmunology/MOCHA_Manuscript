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
setwd('/home/jupyter/scATAC_Supplements')
source('ArchR_Export_ATAC_Data.R')

require(parallel)
require(GenomicRanges)
require(plyranges)
require(data.table)


## boolean whether or not 
## to call peaks 
callPeaks = TRUE
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

cell = cellList[1]
callMacs2Peaks =T

for(cell in cellList){
    print(paste('analyzing', cell))

        setwd('/home/jupyter/MOCHA_Revision/MouseDataAnalysis/Data')
        setwd(cell)

        ############################################################
        ############################################################

        ## change directory
        ## to location where the 
        ## data exist 
    
        setwd('bulk/bedfiles')

        fnames <- dir()
        fnames <- fnames[grep('bedGraph', fnames)]
        
        archLog_idx <- grep('ArchR', fnames)
        if(length(archLog_idx) >0){
            fnames <- fnames[- grep('ArchR', fnames)]
            
        }
        samples <- gsub('_perThousandCells.bedGraph', '', fnames)

        ############################################################
        ############################################################

        if(callMacs2Peaks){
        ## MACS2 Peaks 
        dir.create('../macs2')   
            
        cmd_sample<- lapply(1:length(samples),

           function(x)
                  paste('macs2 callpeak -t', fnames[x],
                
           '-g hs -f BED --nolambda --shift -75 --extsize 150 --broad',
                         '--nomodel -n',

                    paste('../macs2/',samples[x],sep='')
                       )
                                )

        mclapply(cmd_sample,
                 function(x) system(x),
                 mc.cores=30
                 )
        }
        
        ############################################################
        ############################################################

        ## Homer Peaks 
        dir.create('../homer_peaks')   
        
        dir.create('../homer_peaks/homer_dirs')
        outdir <- paste('../homer_peaks/homer_dirs/', samples, sep='')

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
        setwd('../homer_peaks/homer_dirs')
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
