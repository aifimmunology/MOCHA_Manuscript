############################################################

# Author: Samir Rachid Zaim
# Script: Fig. 1 C Normalization

############################################################
############################################################
rm(list=ls())

require(data.table)
require(ggplot2)
require(ggpubr)
require(MOCHA)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(plyranges)
require(parallel)
require(pbapply)

# Load the ArchR Project
ArchRProject <- loadArchRProject("/home/jupyter/FullCovid")
source('/home/jupyter/theme.R')

## export cell types 
## as "bulk" 
cellColName <- 'predictedGroup_Col2.5'

setwd('/home/jupyter/longPilotValidation/sample_specific_analyses/data')

#######################################################################
metadf_dt <- as.data.table(getCellColData(ArchRProject)) 
numCores=10
metaColumn='new_sample_cellType'

### Select 10 cell types
### To Show normalization 
celltypes = unique(metadf_dt$predictedGroup_Co2)[1:10]
blackList = getBlacklist(ArchRProject)

#### subset early visit

lookup_table <- unique(metadf_dt[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]

#### and extract only 
#### 10 from each group 
set.seed(101)
covid_neg <- sample(samplesToKeep[COVID_status =='Negative']$Sample, 10)
covid_pos <- sample(samplesToKeep[COVID_status =='Positive']$Sample, 10)

### 
samplesForAnalysis <- c(covid_neg, covid_pos)
samplesForAnalysis

calculate_lambdas <- function(cellFrags, FinalBins,
                              sampleName
                             ){
    ### Calculate raw and 
    ### Normalized counts 
    
    ### Raw Counts matrix: 
    ### fixed scale across samples
    RawCountsMatrix <- MOCHA:::calculate_intensities(cellFrags, 
                                          FinalBins,
                                          totalFrags= as.integer(10^9)
                                                )
    
    ### Raw Counts matrix: 
    ### fixed scale across samples
    nFrags = as.numeric(length(cellFrags))
    NormalizedMatrix2 <- MOCHA:::calculate_intensities(cellFrags, 
                                          FinalBins,
                                          totalFrags= as.integer(nFrags * 100)
                                          )
    df = data.table::data.table(
        Sample= sampleName,
        Raw=RawCountsMatrix$TotalIntensity,
        Norm=NormalizedMatrix2$TotalIntensity)
    return(df)

}


calculate_intensities <- function(cell, FinalBins=FinalBins){
    
        print(cell)
    
        # Get our fragments for this cellPop
        frags <- getPopFrags(
          ArchRProj = ArchRProject,
          cellPopLabel = 'new_cellType',
          cellSubsets = cell,
          numCores = 45,
          poolSamples = F,
          verbose = FALSE
        )
   
        ### Filter Samples
        ### from Visit 1 
        idx <- grep(paste(samplesForAnalysis, collapse='|'), names(frags))
        
        ### Subset Fragments
        ### for Visit one
        frags = frags[idx]
        
        ### Determine dynamic range
        ### of bins that have fragments 
        ### across all fragment files 
    
        FinalBins <-  MOCHA:::determine_dynamic_range(
            AllFragmentsList = stack(GenomicRanges::GRangesList(frags)),
            blackList = blackList,
            binSize = 500,
            doBin = FALSE
        )
    
        ### Create Cluster object
        ### With exported Items to parallelize
        df_list = mclapply(1:length(frags),
               function(x) calculate_lambdas(frags[[x]],
                                             FinalBins,
                                            names(frags[x])),
                           mc.cores=5
               )

        ### Extract raw and normalized
        ### Intensities for final 
        ### graphical analyses 
    
        full_df = rbindlist(df_list)
        full_df = melt(full_df)
        full_df$Sample_Numeric = as.numeric(factor(full_df$Sample))
    
        ### Filter 0s from TSAM
        full_df = full_df[value > 0]
    
        ### Calculate total normalized
        ### intensities for each sample
        tmp <- full_df[, sum(value), by=list(variable, Sample)] 
        
        ### Create final object
        resObj = list(full_df,
                    tmp)
               
        rm(frags, df_list, Final, full_df, tmp)
        pryr::mem_used()
        gc()
    
        return(resObj)

}

cells = celltypes[1:10]

values = lapply(cells, 
             function(x) calculate_intensities(x)
                )

sample_ints <- rbindlist(lapply(1:10, 
                               function(x)
                                    cbind(values[[x]][[1]], cells[x])
                               )
                               )

total_ints <- rbindlist(lapply(1:10, 
                               function(x)
                                    cbind(values[[x]][[2]], cells[x])
                               )
                               )

total_ints$Count <- ifelse(total_ints$variable=='Raw', 
                           'Raw',
                           'Normalized')
total_ints$Count <- factor(total_ints$Count ,
                           levels=c('Raw', 'Normalized'))

total_ints$Batch = gsub(paste(paste(celltypes[1:5], collapse='|'),'#',sep=''), '', total_ints$Sample)

total_ints$Batch = sapply(total_ints$Batch,
                          function(x)
                              unlist(strsplit(x,'_'))[[1]]
                                                      )

total_ints$Sample = gsub('__[0-9\\.]+','', total_ints$Sample)
total_ints$Sample = gsub(
    paste(unique(total_ints$V2),collapse='|'),
    replacement='', 
    total_ints$Sample)

total_ints$Sample = gsub('#','', total_ints$Sample)

total_ints = dplyr::left_join(total_ints, metadf_dt[,c('Sample','Batch.ID')], multiple='first')

total_ints$Batch = total_ints$Batch.ID

pdf('/home/jupyter/MOCHA_Revision/batch_effect_by_normalization.pdf', height=24, width=8)
p= ggplot(total_ints,
       aes(x=Batch,
           y=V1,
          fill=V2))+geom_point()+geom_boxplot()+facet_wrap(V2~Count, ncol=2 )+
theme_minimal()+
theme(axis.text.x=element_text(size=10, angle=90),
     text=element_text(size=12))+
ylab('Number of total Fragments by Sample')+
ggtitle('Comparison of total fragments by sample Pre and Post Normalization')
print(p)
dev.off()