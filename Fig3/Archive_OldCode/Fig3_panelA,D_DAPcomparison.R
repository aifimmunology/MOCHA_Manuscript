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
archr$Peak = paste(archr$seqnames,
                   ':',
                   archr$start,
                   '-',
                   archr$end,
                   sep='')

seurat_tiles <- sapply(seurat$V1, 
                       function(x){
                           elements = unlist(strsplit(x, '-'))
                           peak = paste(elements[1],
                                        ':',
                                        elements[2],
                                        '-',
                                        elements[3],
                                        sep='')
                           peak
                           }
                       )

seurat$Peak = seurat_tiles

#################################################################
#################################################################

#################################################################
#################################################################

calculateMeanDiff <- function(tmp_mat, group){
    
    a = log2(tmp_mat[, which(group==1)+1, with=F]+1)
    b = log2(tmp_mat[, which(group==0)+1, with=F]+1)
    
    mean_diff = rowMeans(a) - rowMeans(b)
    mean_diff  
}

calculateMeanDiff2 <- function(tmp_mat, group){
    
    a = log2(tmp_mat[, which(group==1)+1, with=F]+1)
    b = log2(tmp_mat[, which(group==0)+1, with=F]+1)
    
    a[a==0] <- NA
    b[b==0] <- NA
    
    mean_diff = rowMeans(a, na.rm=T) - rowMeans(b, na.rm=T)
    mean_diff  
}





common = intersect(intersect(seurat$Peak, scmacs$Peak), archr$Peak)
union =  union(union(seurat$Peak, scmacs$Peak), archr$Peak)

unique_archr = setdiff(archr$Peak, union(seurat$Peak, scmacs$Peak))


unique_scmacs = setdiff(scmacs$Peak, union(seurat$Peak, archr$Peak))



unique_seurat= setdiff(seurat$Peak, union(archr$Peak, scmacs$Peak))




unique_archr <- twoPart[Peak %in% unique_archr]; unique_archr$Model='ArchR'
unique_archr$MeanDiff = calculateMeanDiff2(tileMatrix[tileID %in% unique_archr$Peak], group)

unique_scmacs <- twoPart[Peak %in% unique_scmacs];unique_scmacs$Model='scMACS'
unique_scmacs$MeanDiff = calculateMeanDiff2(tileMatrix[tileID %in% unique_scmacs$Peak], group)

unique_seurat <- twoPart[Peak %in% unique_seurat];unique_seurat$Model='Seurat'
unique_seurat$MeanDiff = calculateMeanDiff2(tileMatrix[tileID %in% unique_seurat$Peak], group)


df = rbind(unique_archr,unique_scmacs,
           unique_seurat)

df$DiffRho = df$Case_rho - df$Control_rho
           
library(VennDiagram)
setwd('../../figure3')
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- c('blue','green','red')

# Chart
venn.diagram(
        x = list(scmacs$Peak, 
                 seurat$Peak, 
                 archr$Peak),
        category.names = c("MOCHA" , "Signac" , "ArchR"),
        filename = 'venn.png',
        output=TRUE,

        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",

        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,

        # Numbers
        cex = .3,
        fontface = "bold",
        fontfamily = "sans",

        # Set names
        cat.cex = 0.3,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)    

df$Model[df$Model=='Seurat'] <-'Signac'
df$Model[df$Model=='scMACS'] <-'MOCHA'

    
df$MedianDiff = df$Case_mu - df$Control_mu
df$AvgIntensity = (df$Case_mu + df$Control_mu)/2

# ## Volcano plot 
png('panelD_qualitative.png')
df$Model = factor(df$Model,
                  levels=c('MOCHA', 'ArchR','Signac'))
ggplot(df,
       aes(x=MeanDiff,
           y= abs(DiffRho),
          col=Model))+
geom_point(size=1) +
facet_wrap(~Model, ncol=1)+
scale_color_manual(values = c('MOCHA'= 'blue',
                              'ArchR' = 'red',
                              'Signac'= 'green'))+ThemeMain+
xlab('Mean Intensity Difference') + ylab('Difference in Proportion of 0s')+
xlim(-3,3)+

theme(legend.position='none')
dev.off()

#################################################################
## quantify regions
## for qualitative analyses 

df[abs(DiffRho)>0.5, .N, by=Model]
df[(MeanDiff)> 0 & abs(DiffRho) < 0.5, .N, by=Model]
df[(MeanDiff)< 0 & abs(DiffRho) < 0.5, .N, by=Model]


#################################################################
#################################################################

### Plot Unique Regions 
### By Each Method
require(ggpubr)
unique_scmacs = unique_scmacs[order(unique_scmacs$P_value,
                                    decreasing=F),]
unique_archr = unique_archr[order(unique_archr$P_value,
                                    decreasing=F),]
unique_seurat = unique_seurat[order(unique_seurat$P_value,
                                    decreasing=F),]
### get top regions 
plot_top_regions <- function(peaks, fname){
    
        p_list <- lapply(1:5,
                         function(x) 
                             plot_differential_region(tileMatrix, 
                                                      tileID=peaks[x],
                                                      group)
                         )
    
        png(fname, width=1200, height=400)

        p = ggarrange(plotlist=p_list,
               ncol=5)
        print(p)
        dev.off()
    
    }

### Get top 5 regions 
plot_top_regions(unique_scmacs$Peak[1:5], fname='panelE_MOCHA_top.png')
plot_top_regions(unique_archr$Peak[1:5], fname='panelE_archr_top.png')
plot_top_regions(unique_seurat$Peak[1:5], fname='panelE_signac_top.png')

#################################################################
#################################################################

### calculate CV on the peaks 
tM = tileMatrix[,2:ncol(tileMatrix)]
tM = log2(tM+1)

tM = tM[tileMatrix$tileID %in% scmacs$Peak]

calculate_CV <- function(tM, removeZeroes=F){
    
    ## transform Zeroes To Nas 
    if(removeZeroes){
        tM[tM==0] <- NA   
    }
    
    ## transform to matrix
    tM = as.matrix(tM)
    
    ## calculate mean & sd vectors 
    mu_vec = rowMeans(tM, na.rm=T)
    sd_vec = rowSds(tM, na.rm=T)
    
    ## return coefficient of variation
    CVs = sd_vec / mu_vec
    return(CVs)
    
}
cvAll = calculate_CV(tM, F)
cvNonzeroes=calculate_CV(tM, T)

png('CVs.png')

par(mfrow=c(2,1))
hist(cvAll, main='CV All values',xlab='CV', breaks=50)
hist(cvNonzeroes, main='CV Nonzero values',xlab='CV', breaks=50)
dev.off()

scmacs$CV = cvAll
scmacs$CVNonzero=cvNonzeroes

write.csv(scmacs,
          file='cd16_daps_withCV.csv',
          row.names=F)

