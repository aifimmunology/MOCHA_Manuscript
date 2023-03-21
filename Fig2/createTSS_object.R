################################################################################
################################################################################

### create code to generate 
### the TSS sites 
### for figure 2 analyses 

################################################################################
################################################################################
require(dplyr)
require(plyranges)

setwd('/home/jupyter/MOCHA_Manuscript/Fig2')

load('tss_reorganized.RDS')

### HG38 
library(TxDb.Hsapiens.UCSC.hg38.refGene)
TxDb_HG38 <- TxDb.Hsapiens.UCSC.hg38.refGene

tss_hg38 <- ensembldb::transcriptsBy(TxDb_HG38, by = ("gene")) %>% IRanges::stack(.)

tss_hg38 = promoters(tss_hg38, upstream = 1, downstream = 1)
save(tss_hg38,
     file='TSS_HG38.RDS')


### HG19 
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.refGene)

TxDb_HG19 <- TxDb.Hsapiens.UCSC.hg19.refGene
tss_hg19 <- ensembldb::transcriptsBy(TxDb_HG19, by = ("gene")) %>% IRanges::stack(.)

tss_hg19 = promoters(tss_hg19, upstream = 1, downstream = 1)

save(tss_hg19,
     file='TSS_HG19.RDS')