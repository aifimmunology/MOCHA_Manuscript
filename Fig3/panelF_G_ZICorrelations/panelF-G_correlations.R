#################################################################
#################################################################
###
### Script: cd16_scMACS_comparison
### Description: Conducts venn diagrams
###              and peak quality comparisons 
### Author: Samir Rachid Zaim, PhD
###
#################################################################
#################################################################

## Load Libraries
require(scMACS)
require(data.table)
require(VennDiagram)
require(ggplot2)
require(ArchR)
require(ggdist)
require(ggpubr)
source('/home/jupyter/theme.R')

covidArchR <- loadArchRProject("/home/jupyter/FullCovid")
metadf <- getCellColData(covidArchR)
metadf_dt <- as.data.table(getCellColData(covidArchR)) 


#################################################################
#################################################################
lookup_table <- unique(metadf[,c('Sample',
                                 'COVID_status',
                                 'new_cellType_sample',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

visit_1 <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms),]

## filter out the cd16 monocytes population 
metaFile <- visit_1[grep('CD16 Mono',visit_1$new_cellType_sample),]
metaFile$SampleCellType <- metaFile$new_cellType_sample
metaFile$Class =  metaFile$COVID_status
covidPats= metaFile$SampleCellType[metaFile$Class =='Positive']

#################################################################
#################################################################
setwd('/home/jupyter/covid/scMACS_manuscript_analyses/sample_specific_analyses/')
setwd('dap_cd16')
archr=fread('CD16_DAPs_ArchR_on_scMACSTiles.csv')
seurat = fread('Seurat_204K_PeaksTest.csv')
scmacs = fread('cd16_scmacs.csv')


setwd('../daps/backup')
load('CD16 Mono.RDS')
tileMatrix = results$TileMatrix
seurat = seurat[p_val_adj <= 0.05]
twoPart =results$DAPs
group = ifelse(colnames(tileMatrix)[2:ncol(tileMatrix)] %in% covidPats,
               1,0)

#################################################################
#################################################################
numCores=50
chr4 = tileMatrix[grep('chr4', tileMatrix$tileID),]
start_end <- gsub('chr4:','' , chr4$tileID)
start = as.numeric(gsub('-[0-9]*','',start_end))
idx = order(start, decreasing=F)
ordered_chr4 = chr4[idx,]

#################################################################
#################################################################
mat = as.matrix(ordered_chr4[1:200,-1])

### generate all pairwise combinations
N = nrow(mat)
pairwise_combos = expand.grid(1:N, 1:N)

zi_spearman <- co_accessibility(mat, numCores=50)
### Loop through all pairwise 
### combinations of peaks 
zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
         function(x)
             weightedZISpearman(x=mat[pairwise_combos$Var1[x],],
                             y=mat[pairwise_combos$Var2[x],]
                             ),
         mc.cores=numCores
         ))


spearman_correlation <- unlist(mclapply(1:nrow(pairwise_combos),
         function(x)
             cor(x=mat[pairwise_combos$Var1[x],],
                 y=mat[pairwise_combos$Var2[x],],
                 method='spearman'
                             ),
         mc.cores=numCores
         ))
### Create zero-inflated correlation matrix
### from correlation values, 
zi_spear_mat <- data.frame(ZI_Correlation=zero_inflated_spearman,
                           Spearman =spearman_correlation,
                           Peak1= ordered_chr4$tileID[1:200][pairwise_combos$Var1],
                           Peak2= ordered_chr4$tileID[1:200][pairwise_combos$Var2]
                           )

#################################################################
#################################################################
setwd('../correlation_results')


zi_spear_mat$diff = zi_spear_mat$ZI_Correlation - zi_spear_mat$Spearman
    
idx = which(abs(zi_spear_mat$diff) > 0.45)
plot_coaccessibility(zi_spear_mat, idx[10], tileMatrix, fname='large_diff1',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx[2], tileMatrix, fname='large_diff2',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx[3], tileMatrix, fname='large_diff3',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx[4], tileMatrix, fname='large_diff4',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx[15], tileMatrix, fname='large_diff5',addLabel=T)



idx2 = which(zi_spear_mat$ZI_Correlation < -.2 & zi_spear_mat$Spearman >0 |
             zi_spear_mat$ZI_Correlation > 0.2 & zi_spear_mat$Spearman < 0)
plot_coaccessibility(zi_spear_mat, idx2[1], tileMatrix, fname='sign_change1',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx2[2], tileMatrix, fname='sign_change2',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx2[10], tileMatrix, fname='sign_change3',addLabel=T)
plot_coaccessibility(zi_spear_mat, idx2[6], tileMatrix, fname='sign_change4',addLabel=T)


zi_spear_mat$Group <- ''
zi_spear_mat$Group[idx] <- 'Spurrious_Correlations'
zi_spear_mat$Group[idx2] <- 'Sign_Change'

png('ZI_vs_spearman.png')
ggplot(zi_spear_mat,
       aes(x=ZI_Correlation,
       y=Spearman))+
        geom_point(color='gray')+
            geom_point(data = zi_spear_mat %>% filter(Group == "Spurrious_Correlations"), 
               color = "red",size=3)+
        geom_point(data = zi_spear_mat %>% filter(Group == "Sign_Change"), 
               color = "blue", size=3)
dev.off()
#################################################################
#################################################################


### Quantify the results 
N = nrow(zi_spear_mat)
a <- sum(zi_spear_mat$ZI_Correlation >0 & zi_spear_mat$Spearman <0)
b <- sum(zi_spear_mat$ZI_Correlation <0 & zi_spear_mat$Spearman >0)


c <- sum( abs(zi_spear_mat$ZI_Correlation - zi_spear_mat$Spearman)>0.25)
d <- sum( abs(zi_spear_mat$ZI_Correlation - zi_spear_mat$Spearman)>0.3)