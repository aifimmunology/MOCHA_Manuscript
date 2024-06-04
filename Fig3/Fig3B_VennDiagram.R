#################################################################
#################################################################
###
### Script: 
### Description: Conducts venn diagrams
###              and peak quality comparisons 
###
#################################################################
#################################################################

## Load Libraries
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

