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

#################################################################
#################################################################
setwd('/home/jupyter/MOCHA_Manuscript/Fig3')
setwd('generateData')

archr =fread('cd16_ArchR.csv')
signac = fread('cd16_signac.csv');signac = signac[p_val_adj <= 0.05]
mocha = fread('cd16_mocha.csv')
allTiles = fread('cd16_allTiles.csv')
tileMatrix = fread('tileMatrix.csv')

#################################################################
#################################################################
archr$Tile = paste(archr$seqnames,
                   ':',
                   archr$start,
                   '-',
                   archr$end,
                   sep='')

signac$Tile <- sapply(signac$V1, 
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


#################################################################
#################################################################

common = intersect(intersect(signac$Tile, mocha$Tile), archr$Tile)
union =  union(union(signac$Tile, mocha$Tile), archr$Tile)

unique_archr = setdiff(archr$Tile, union(signac$Tile, mocha$Tile))
unique_mocha = setdiff(mocha$Tile, union(signac$Tile, archr$Tile))
unique_signac= setdiff(signac$Tile, union(archr$Tile, mocha$Tile))


### 

unique_archr <- allTiles[Tile %in% unique_archr]; unique_archr$Model='ArchR'

unique_mocha <- allTiles[Tile %in% unique_mocha];
unique_mocha$Model='scMACS'

unique_signac <- allTiles[Tile %in% unique_signac];
unique_signac$Model='signac'


df = rbind(unique_archr,unique_mocha,
           unique_signac)

df$DiffRho = df$Case_rho - df$Control_rho
           
library(VennDiagram)
setwd('../panelA_Venn/')
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- c('blue','green','red')

# Chart
venn.diagram(
        x = list(mocha$Tile, 
                 signac$Tile, 
                 archr$Tile),
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





















