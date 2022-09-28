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


### get unique regions all methods 

unique_archr <- allTiles[Tile %in% unique_archr]; 
unique_archr$Model='ArchR'

unique_mocha <- allTiles[Tile %in% unique_mocha];
unique_mocha$Model='MOCHA'

unique_signac <- allTiles[Tile %in% unique_signac];
unique_signac$Model='Signac'


df = rbind(unique_archr,unique_mocha,
           unique_signac)

df$DiffRho = abs(df$Pct0_Case - df$Pct0_Control)
df$MedianDiff = df$Avg_Intensity_Case - df$Avg_Intensity_Control
df$AvgIntensity = (df$Avg_Intensity_Case + df$Avg_Intensity_Control)/2


#################################################################
# ## Volcano plot 
png('../panelD-E_uniqueRegions/panelD_qualitative.png')
df$Model = factor(df$Model,
                  levels=c('MOCHA', 'ArchR','Signac'))

ggplot(df,
       aes(x=Log2FC_C,
           y= abs(DiffRho),
          col=Model))+
    geom_point(size=1) +
    facet_wrap(~Model, ncol=1)+
    scale_color_manual(values = c('MOCHA'= 'blue',
                              'ArchR' = 'red',
                              'Signac'= 'green'))+
    ThemeMain+
        xlab('Median Intensity Difference') + 
        ylab('Difference in Proportion of 0s')+
        ylim(0,1)+
    theme(legend.position='none')
dev.off()

#################################################################
## quantify regions
## for qualitative analyses 

df[abs(DiffRho)>0.5, .N, by=Model]
df[(MedianDiff)> 0 & abs(DiffRho) < 0.5, .N, by=Model]
df[(MedianDiff)< 0 & abs(DiffRho) < 0.5, .N, by=Model]


#################################################################
#################################################################
ArchRProj = loadArchRProject('/home/jupyter/FullCovid')
metadata = as.data.table(ArchRProj@cellColData)


## Get metadata information
## at the sample level
lookup_table <- unique(metadata[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
lookup_table <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]

## subset ArchR Project
group = ifelse(colnames(tileMatrix)[1:39] %in% 
               lookup_table$Sample[lookup_table$COVID_status=='Positive'],
               1,
               0)

#################################################################
#################################################################



### Plot Unique Regions 
### By Each Method
require(ggpubr)
unique_mocha = unique_mocha[order(unique_mocha$P_value,
                                    decreasing=F),]
unique_archr = unique_archr[order(unique_archr$P_value,
                                    decreasing=F),]
unique_seurat = unique_signac[order(unique_signac$P_value,
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
setwd('../panelD-E_uniqueRegions/')
plot_top_regions(unique_mocha$Tile[1:5], fname='panelE_MOCHA_top.png')
plot_top_regions(unique_archr$Tile[1:5], fname='panelE_archr_top.png')
plot_top_regions(unique_signac$Tile[1:5], fname='panelE_signac_top.png')

#################################################################
#################################################################
