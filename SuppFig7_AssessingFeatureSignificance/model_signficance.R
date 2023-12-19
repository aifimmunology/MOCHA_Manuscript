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
source('/home/jupyter/covid/scMACS_manuscript_analyses/trainModel/helper_functions.R')
#devtools::install('/home/jupyter/scMACS')
require(data.table)
require(ggplot2)
require(ggpubr)
require(scMACS)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)
require(parallel)
source('/home/jupyter/theme.R')

# Load the ArchR Project
covidArchR <- loadArchRProject("/home/jupyter/FullCovid")

# Create a metadata column that has your groups of interest. 
metadf <- getCellColData(covidArchR) 
df = as.data.table(metadf)

## get cell types to export
cellType <- metadf$predictedGroup_Co2
cellType_sample <- paste(metadf$Sample, metadf$predictedGroup_Co2,sep='__')

### load NK fragment files 
### for entire NK population
cell <- c('NK')
frags <- getPopFrags(covidArchR, metaColumn = "predictedGroup_Co2", 
                     cellSubsets = cell, 
                     numCores= 30)
NK_population <- frags[[1]]

###########################################################
###########################################################

## create directory 
setwd('/home/jupyter/covid/scMACS_manuscript_analyses/trainModel/training_data')

## load training data
ground_truth <- fread('/home/jupyter/covid/scMACS_manuscript_analyses/trainModel/NK_trainingset.csv')

## generate the genomic 
## region of interest for the
## analyses 
FinalBins <- scMACS::determine_dynamic_range(SimpleList(frags),
                                             covidArchR,
                                             binSize=500, 
                                             doBin=FALSE)
#rm(frags)
pryr::mem_used()

###########################################################
###########################################################

## obtain_coefficients_for_subsample
# This function estimates model
# coefficients for estimating 
# accessibility at a given simulated
# cell population. This process is 
# replicated X times to ensure 
# robust coefficients and modeling.

obtain_coefficients_for_subsample <- function(NK_population,
                                              numCells,
                                              numReplicates,
                                              numCores=10,
                                              normScale=10^9                                             ){
    
    ### Sample a set of cells for each 
    ### simulated population and extract
    ### the corresponding fragments
    cell_subsets <- subsample_cells(NK_population, 
                                              cell_size=numCells,
                                              numReplicates=numReplicates,
                                    numCores=numCores
                                   )

    ### Calculate the normalizing factor 
    ### for each simulated population
    totalFrags_list <- sapply(cell_subsets, length)
    

    ### Calculate the normalized 
    ### intensities for each subsample
    
    intensities_list <- mclapply(1:length(cell_subsets),
           function(x) 
               scMACS::calculate_intensities(fragMat=cell_subsets[[x]],
                                      candidatePeaks=FinalBins,
                                      totalFrags=as.numeric(totalFrags_list[x])
                                     ),
                                 mc.cores=numCores
                                 )
    

    ### Train model only on nonzero
    ### bins with reads 
    intensities_list <- mclapply(intensities_list,
                                 function(x){
                                    x[TotalIntensity>0]
                                     }
                                 ,
                                 mc.cores=numCores
                                 )
    
    intensities_list <- intensities_list[sapply(intensities_list, 
                                                function(x) 
            paste(class(x), collapse='_')
                                                ) =='data.table_data.frame']
    
    
    #print(intensities_list[[1]])

    ### add the ground truth from 
    ### Macs2' NK pop prediction

    intensities_list <- mclapply(1:length(intensities_list),
                                 function(x) 
                                     add_label(intensities_list[[x]], ground_truth),
                                 mc.cores=numCores
                                 )
    
    ### estimate model coefficients for 
    ### accessibility across all replicates
    
    model_raw <- mclapply(intensities_list, 
               function(x)
                   trainModel_2(x),
                 mc.cores=numCores
               )
   
    ### Extract Coefficients 
    model_coefficients_raw <- lapply(model_raw,
                                     function(x) coef(x)
                                     )
    model_coefficients_raw = data.table(do.call(rbind, model_coefficients_raw))
    
    model_coefficients_raw = melt(model_coefficients_raw)
    Z = model_coefficients_raw[,mean(value) / sd(value), by=variable]
    Z = data.frame(t(Z))
    colnames(Z) = Z[1,]
    Z = Z[2,]
    Z$NumCells = numCells
    Z$NumReplicates = numReplicates
 
    return(list(Coefficients=model_coefficients_raw,
                ZScore =Z))
    
}

###########################################################
###########################################################

# Repeat Analysis for different cell subsamples 
# Conduct 50 different replicates
# For each 'subsampled' model to ensure
# model stability

model10 = obtain_coefficients_for_subsample(NK_population,
                                              numCells=10,
                                              numReplicates=100,
                                              numCores=5)


model50 = obtain_coefficients_for_subsample(NK_population,
                                              numCells=50,
                                              numReplicates=100,
                                              numCores=5)


model100 = obtain_coefficients_for_subsample(NK_population,
                                              numCells=100,
                                              numReplicates=100,
                                              numCores=5)


setwd('/home/jupyter/MOCHA_Revision/model_coefficients')

plotFeatureSignificance = function(feature, coef_df){
    df = coef_df$Coefficients
    df = df[variable==feature]
    
    xmin = min(df$value)-sd(df$value)
    xmax = max(df$value)+sd(df$value)
    
    p = ggplot(df,
           aes(x=value))+geom_histogram()+
        xlim(c(xmin,xmax))+ggtitle(paste(feature, '95% Conf Interval'))+
        theme_minimal()+
        theme(text=element_text(size=14))+
        geom_vline(xintercept=quantile(df$value, 0.025), linewidth=2, lty=2)+
        geom_vline(xintercept=quantile(df$value, 0.975), linewidth=2, lty=2)
    return(p)
}

pdf('10coefficients.pdf')

plotlist10 = lapply(unique(model10$Coefficients$variable),
                    function(x) 
                        plotFeatureSignificance(x, model10))

print(ggpubr::ggarrange(plotlist=plotlist10, ncol=1))
dev.off()
                    
pdf('50coefficients.pdf')

plotlist50 = lapply(unique(model50$Coefficients$variable),
                    function(x) 
                        plotFeatureSignificance(x, model50))

print(ggpubr::ggarrange(plotlist=plotlist50, ncol=1))
dev.off()                    
                    
pdf('100coefficients.pdf')

plotlist100 = lapply(unique(model100$Coefficients$variable),
                    function(x) 
                        plotFeatureSignificance(x, model100))

print(ggpubr::ggarrange(plotlist=plotlist100, ncol=1))
dev.off()
                     
       
                     
                     
                     

model500 = obtain_coefficients_for_subsample(NK_population,
                                              numCells=500,
                                              numReplicates=100,
                                              numCores=5)                     
pdf('500coefficients.pdf')
plotlist500 = lapply(unique(model500$Coefficients$variable),
                    function(x) 
                        plotFeatureSignificance(x, model500))

print(ggpubr::ggarrange(plotlist=plotlist500, ncol=1))
dev.off()

                                

           