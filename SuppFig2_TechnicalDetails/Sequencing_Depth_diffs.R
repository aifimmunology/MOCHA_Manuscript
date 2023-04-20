############################################################
############################################################

# Author: Samir Rachid Zaim
# Script: Assessing dropout across groups

############################################################
############################################################

require(TxDb.Hsapiens.UCSC.hg38.refGene)
require(data.table)
require(ggplot2)
require(MOCHA)
library(ArchR)
library(GenomicRanges)
require(dplyr)
require(parallel)
library(plyranges)
source('/home/jupyter/theme.R')
source('/home/jupyter/MOCHA_Manuscript/Fig2/helper_granges.R')

############################################################
############################################################

## load ArchR Project 
ArchRProj <- ArchR::loadArchRProject('/home/jupyter/FullCovid/')

## Extract Metadata and 
## required params for 
## calling open tiles 
metadf <- getCellColData(ArchRProj) 
studySignal = median(metadf$nFrags)
celltypes <- c('CD16 Mono')
############################################################
############################################################

### load scMACS peaks
### for plotting and 
### analyses 
setwd('/home/jupyter/covid/scMACS_manuscript_analyses/data')
metadf=as.data.table(metadf)

sample_allocation = metadf[metadf$predictedGroup_Co2=='CD16 Mono', 
       list(Group=first(COVID_status),
            Visit=first(Visit),
            nFrags=sum(nFrags)), 
            by=Sample]

sample_allocation = sample_allocation %>% arrange(nFrags)

### Controls have 
### more fragments 
sample_allocation[, median(nFrags),by=Group]
sample_allocation$Group2 = ifelse(sample_allocation$Group=='Positive','+','-')
sample_allocation$GroupStatus = paste('COVID', sample_allocation$Group2, sep='')

###
pdf('/home/jupyter/MOCHA_Manuscript/SuppFig4_ZI_differentials/nFrags_by_group.pdf')
ggplot(sample_allocation,
       aes(x=nFrags))+geom_histogram(bins=40)+
        facet_wrap(~GroupStatus, ncol=1,
                  scales='free_y')+
        theme_minimal()+
        xlab('Fragments per Sample')+
        ylab('# of Samples')+
        theme(text=element_text(size=14))
dev.off()

write.csv(sample_allocation,
          file='/home/jupyter/MOCHA_Manuscript/SuppFig4_ZI_differentials/fragment_counts.csv')