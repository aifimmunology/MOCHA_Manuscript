# ################################################################
# ################################################################
# ##
# ## Script: cd16_scMACS_comparison
# ## Description: Conducts venn diagrams
# ##              and peak quality comparisons 
# ## Author: Samir Rachid Zaim, PhD
# ##
# ################################################################
# ################################################################

## Load Libraries
#require(scMACS)
require(MOCHA)
require(data.table)
require(VennDiagram)
require(ggplot2)
require(ArchR)
require(parallel)
source('/home/jupyter/theme.R')
source('/home/jupyter/scMACS/R/plot_coaccessibility.R')

setwd('/home/jupyter/MOCHA_Manuscript/Fig3')
tileMatrix = fread('generateData/tileMatrix.csv', data.table=F)

# ################################################################
# ################################################################

numCores=50
chr4 = tileMatrix[grep('chr4', tileMatrix$tileID),]
start_end <- gsub('chr4:','' , chr4$tileID)
start = as.numeric(gsub('-[0-9]*','',start_end))
idx = order(start, decreasing=F)
ordered_chr4 = chr4[idx,]

# ################################################################
# ################################################################

tileID = ordered_chr4$tileID
start_end = sapply(tileID,
                function(x)
                    strsplit(x,':')[[1]][2]
                   )
start = as.numeric(sapply(start_end,
               function(x)
                   strsplit(x,'-')[[1]][1]
               ))

## find all tiles greater than 1mbp
idx = which(start - start[1] > 1000000)

## Select 1st one to establish
## 1 MBP neighborhood
mat = ordered_chr4[1:idx[1], ]
mat = mat[,!colnames(mat) %in% 'tileID']
mat = as.data.frame(mat)
row.names(mat) <-  ordered_chr4$tileID[1:idx[1]]

### generate all pairwise combinations
N = nrow(mat)
pairwise_combos = expand.grid(1:N, 1:N)
mat = mat[,!colnames(mat) %in% 'tileID']

### Loop through all pairwise 
### combinations of peaks 
zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
         function(x)
             MOCHA:::weightedZISpearman(x=mat[pairwise_combos$Var1[x],1:39],
                             y=mat[pairwise_combos$Var2[x],1:39]
                             ),
         mc.cores=numCores
         ))


spearman_correlation <- unlist(mclapply(1:nrow(pairwise_combos),
         function(x)
             cor(x=as.numeric(mat[pairwise_combos$Var1[x],]),
                 y=as.numeric(mat[pairwise_combos$Var2[x],]),
                 method='spearman'
                             ),
         mc.cores=numCores
         ))

# ################################################################
# ################################################################

# ## Create zero-inflated correlation matrix
# ## from correlation values,

zi_mat <- data.frame(ZI_Correlation=zero_inflated_spearman,
                           Spearman =spearman_correlation,
                           Peak1= ordered_chr4$tileID[1:idx[1]][pairwise_combos$Var1],
                           Peak2= ordered_chr4$tileID[1:idx[1]][pairwise_combos$Var2]
                           )

cor.test(zi_mat$ZI_Correlation, zi_mat$Spearman,
         method='spearman')
#################################################################
#################################################################
setwd('/home/jupyter/MOCHA_Manuscript/SuppFig6_Correlations')

weightedZISpearman = MOCHA:::weightedZISpearman
zi_mat$diff = zi_mat$ZI_Correlation - zi_mat$Spearman

### plot major differences 
### in magnitude 
idx = which(abs(zi_mat$diff) > 0.45)
plot_coaccessibility(zi_mat, idx[9], tileMatrix, fname='large_diff1',addLabel=T)
plot_coaccessibility(zi_mat, idx[2], tileMatrix, fname='large_diff2',addLabel=T)
plot_coaccessibility(zi_mat, idx[3], tileMatrix, fname='large_diff3',addLabel=T)
plot_coaccessibility(zi_mat, idx[4], tileMatrix, fname='large_diff4',addLabel=T)
plot_coaccessibility(zi_mat, idx[15], tileMatrix, fname='large_diff5',addLabel=T)



### plot major differences 
### in reversing direction  
idx2 = which(zi_mat$ZI_Correlation < -.2 & zi_mat$Spearman >0 |
             zi_mat$ZI_Correlation > 0.2 & zi_mat$Spearman < 0)
plot_coaccessibility(zi_mat, idx2[1], tileMatrix, fname='sign_change1',addLabel=T)
plot_coaccessibility(zi_mat, idx2[2], tileMatrix, fname='sign_change2',addLabel=T)
plot_coaccessibility(zi_mat, idx2[10], tileMatrix, fname='sign_change3',addLabel=T)
plot_coaccessibility(zi_mat, idx2[6], tileMatrix, fname='sign_change4',addLabel=T)


zi_mat$Group <- ''
zi_mat$Group[idx] <- 'Spurrious_Correlations'
zi_mat$Group[idx2] <- 'Sign_Change'

pdf('ZI_vs_spearman.pdf')
ggplot(zi_mat,
       aes(x=ZI_Correlation,
       y=Spearman))+
        geom_point(color='gray')+
            geom_point(data = zi_mat[zi_mat$Group == "Spurrious_Correlations",], 
               color = "red",size=3)+
        geom_point(data = zi_mat[zi_mat$Group == "Sign_Change",], 
               color = "blue", size=3)
dev.off()
#################################################################
#################################################################


### Quantify the results 
N = nrow(zi_mat)
flipped_a <- sum(zi_mat$ZI_Correlation >0 & zi_mat$Spearman <0)
flipped_b <- sum(zi_mat$ZI_Correlation <0 & zi_mat$Spearman >0)

c(flipped_a+flipped_b)/N

magnitude <- sum( abs(zi_mat$ZI_Correlation - zi_mat$Spearman)>0.25)
magnitude/N

#################################################################
#################################################################
setwd('/home/jupyter/MOCHA_Manuscript/SuppFig6_Correlations')

write.csv(zi_mat,
          file='cd16_within_correlations.csv')

final_zi_indices <- c(idx[c(2:4, 9,15)],
                      idx2[c(1,2,6,10)])

peak_names <- c(zi_mat[final_zi_indices,]$Peak1,
                zi_mat[final_zi_indices,]$Peak2)

write.csv(tileMatrix[tileMatrix$tileID %in% peak_names,],
          file='exemplar_correlations.csv')
