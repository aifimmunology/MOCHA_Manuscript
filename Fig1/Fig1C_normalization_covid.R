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
          metaColumn = 'new_cellType',
          cellSubsets = cell,
          region = NULL,
          numCores = 45,
          sampleSpecific = TRUE,
          NormMethod = "nfrags",
          blackList = NULL,
          verbose = FALSE,
          overlapList = 50
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
        cl <- makeCluster(20)
        clusterExport(cl, c("calculate_lambdas", "frags",
                      "FinalBins"),
                      envir=environment())

        df_list = parLapply(1:length(frags),
               function(x) calculate_lambdas(frags[[x]],
                                             FinalBins,
                                            names(frags[x])),
                           cl=cl
               )
        stopCluster(cl)

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

################################################################
################################################################
#### plot the results 
setwd('/home/jupyter/MOCHA_Manuscript/Fig1/')

### Plot 1 = total across cell types 
pdf('plots/Fig1A_cellTypeNormalization.pdf', width=5, height=4)
ggplot(total_ints,
       aes(x=V2,
       y=V1))+geom_point()+geom_violin()+
            facet_wrap(~Count, scales='free_y', ncol=1)+
         xlab('')+ylab('')+    theme_minimal()+

        theme(strip.text = element_text(size=14),
              axis.text.x = element_text(size=14, angle=90),
              axis.text.y = element_text(size=14))
dev.off()

write.csv(total_ints, 
          file='supplementalFiles/Fig1C_celltype_normalizations.csv',
         row.names=F)

cv_sample_celltype = total_ints[, sd(V1)/mean(V1), 
                                by=list(Count,V2)]
cv_sample_celltype[, mean(V1), by=Count]


### Set Factor Levels
sample_ints$Count <- ifelse(sample_ints$variable=='Raw', 
                           'Raw',
                           'Normalized')
sample_ints$Count <- factor(sample_ints$Count ,
                           levels=c('Raw', 'Normalized'))

### Choose One Cell Type
### To Show Sample Normalization
cd16_mat <- sample_ints[V2=='CD16 Mono']
cd16_sum <- cd16_mat[, sum(value), by=list(Sample, variable)]
cd16_sum$Sample <- gsub('NK#X001_|NK#X002_|__[0-9]\\.[0-9]+', '', cd16_sum$Sample)

cd16_sum$variable <- ifelse(cd16_sum$variable=='Raw',
                          'Raw',
                          'Normalized')

cd16_sum$variable <- factor(cd16_sum$variable ,
                           levels=c('Raw', 'Normalized'))

### Plot 2 = total across samples in a cell type 
pdf('plots/Fig1C_SampleNormalization.pdf', width=4, height=8)
cd16_sum$Sample <- gsub('CD16_Mono#|\\#' ,'', cd16_sum$Sample)

ggplot(cd16_sum,
       aes(y=as.character(Sample),
           x=V1
          ))+geom_bar(scale='width', stat='identity')+
            facet_wrap(~variable, scales='free_x', ncol=2)+
         xlab('')+ylab('')+theme_minimal()+
        theme(strip.text = element_text(size=25),
              axis.text.x = element_text(size=14, angle=90),
              axis.text.y = element_blank())
    
dev.off()

write.csv(cd16_sum,
          file='supplementalFiles/Fig1C_sample_normalizations.csv',
          row.names=F)

## Quantify Coefficient of Variation
## in total intensity
cd16_sum[, sd(V1)/mean(V1), by=variable]