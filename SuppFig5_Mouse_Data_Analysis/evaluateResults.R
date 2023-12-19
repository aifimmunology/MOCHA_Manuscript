# ###########################################################
# ###########################################################

# Author: Samir Rachid Zaim

# Desc: 
#     This script calls peaks using
#     macs2, homer, and hmmratac 
#     to create those comparisons 
#     with scMACS.

# ###########################################################
# ###########################################################
## load libraries 
require(MOCHA)
require(parallel)
require(ggpubr)
require(data.table)
require(ggplot2)
require(ArchR)
require(MultiAssayExperiment)
require(RaggedExperiment)
require(parallel)
require(GenomicRanges)
require(plyranges)
require(data.table)

## Identify TSS Sites 
## to analyze 
require('TxDb.Mmusculus.UCSC.mm9.knownGene')
txdb = TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
tss = promoters(genes(txdb), upstream = 1, downstream = 1)

## boolean whether or not 
## to call peaks 
callPeaks = TRUE
cellList = c('Astrocytes#PreFrontalCortex_62216',
             'Cardiomyocytes#HeartA_62816', 
             'Enterocytes#LargeIntestineA_62816',
             'Cerebellar_granule_cells#Cerebellum_62216',
             'Collecting_duct#Kidney_62016', 
             'Endothelial_I_cells#Liver_62016',
             'Podocytes#Kidney_62016',
             'Sperm#Testes_62016',
             'Enterocytes#SmallIntestine_62816',
             'T_cells#BoneMarrow_62216', 
             'NK_cells#BoneMarrow_62216', 
             'Microglia#WholeBrainA_62216',
             'Purkinje_cells#WholeBrainA_62216',
             'Hematopoietic_progenitors#BoneMarrow_62216',
             'Regulatory_T_cells#BoneMarrow_62216',
             'Proximal_tubule#Kidney_62016',
             'Alveolar_macrophages#Lung1_62216',
             'Alveolar_macrophages#Lung2_62216',             
             'Type_I_pneumocytes#Lung2_62216',
             'Type_I_pneumocytes#Lung1_62216')



cell = cellList[1]

### set directory 
setwd('/home/jupyter/MOCHA_Manuscript2/Fig2/')
source('../theme.R')
source('helper_granges.R')
source('utils.R')
setwd('../../MOCHA_Revision/MouseDataAnalysis')
ctcf <- plyranges::read_bed('All_Mouse_CTCF.bed')


### set directory to revisions 
homeDir = '/home/jupyter/MOCHA_Revision/MouseDataAnalysis/'
setwd(homeDir)
## load ArchR project
ArchRProj = loadArchRProject('/home/jupyter/MOCHA_Revision/MouseData')
metadata = as.data.table(ArchRProj@cellColData)

################################################################
### load MOCHA Tiles 
datasetDir = "/home/jupyter/MOCHA_Revision/MouseDataAnalysis/Data"
setwd(datasetDir)
cell_dirs <-  dir()

cell_dir = cell_dirs[1]

################################################################
################################################################
### panelA cell counts 

cells_per_sample <- metadata[
    CellType_Organ %in% cell_dirs, 
    list(Count= .N), 
    by= list(CellType_Organ)]

cells_per_sample = cells_per_sample %>% arrange(Count)

cells_per_sample$CellType_Organ = factor(
    cells_per_sample$CellType_Organ,
    levels=cells_per_sample$CellType_Organ)

dir.create('../results')
pdf('../results/cellsCounts_violin_h.pdf', height=6,width=7)
    p=ggplot(cells_per_sample,
       aes(x=CellType_Organ,
           y=Count))+geom_bar(stat='identity')+ 
            xlab('Cell and Organ Sample')+ylab('# of Cells')+
            scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                  axis.text.x=element_text(size=12, angle=90),
                 axis.text.y=element_text(size=16))+
            scale_fill_cell()
   print(p)
   
dev.off()

sizeText = 17
lwd=2

################################################################
################################################################

### summarize across 
### the selected cell types 
setwd(datasetDir)
total_res <- lapply(cell_dirs, 
                    function(x){
                        print(x)
                        setwd(x)
                        setwd('bulk')
                        setwd('MOCHA')
                        fnameMocha = paste(x, '_peaks.rds',sep='')
                        tsam <- readRDS(fnameMocha)
                        mocha_tile = tsam[peak==T]
                        
                        setwd('../macs2')
                        ### Macs2 Peaks
                        fnames <- dir()
                        fnames <- fnames[grep('broadPeak', fnames)]

                        ## read summits list
                        macs2_tile_list <- fread(fnames, data.table=T)
                        macs2_tile_list =  peak_to_tile_macs2_bulk(1, 
                                                                   list(macs2_tile_list))
                        
                        ### remove tiles extending 
                        ### onto tiles with no fragment 
                        macs2_tile_list = macs2_tile_list[macs2_tile_list$tileID %in% tsam$tileID]
                        
                        setwd('../homer_peaks/homer_dirs')
                        ### HOMER Peaks
                        fnames <- dir()
                        peakNames <- paste(fnames,'/peaks.txt', sep='')
                        peakNames <- peakNames[!peakNames %in% c("homer_peaks/peaks.txt",
                                                                 "ArchRLogs/peaks.txt",
                                                                 '-style/peaks.txt')]

                        ### extract tiles 
                        homer_peak_list <- fread(peakNames)
                        homer_peak_list = peak_to_tile_homer_bulk(1, 
                                                                  list(homer_peak_list))
                        
                        ### remove tiles extending 
                        ### onto tiles with no fragment 
                        homer_peak_list = homer_peak_list[homer_peak_list$tileID %in% tsam$tileID]
	                               
                        ### convert objects to Granges
                        mocha_gr <- makeGRangesFromDataFrame(mocha_tile, 
                                                             keep.extra.columns=T)

                        macs2_gr <- makeGRangesFromDataFrame(macs2_tile_list, 
                                                             keep.extra.columns=T)

                        homer_gr <- makeGRangesFromDataFrame(homer_peak_list, 
                                                             keep.extra.columns=T)    
                        
                        tile_res = data.frame(
                            MOCHA= length(mocha_gr),
                            MACS2= length(macs2_gr),
                            HOMER= length(homer_gr)
                        )
                        
                        tile_res$CellPopulation <- x
                        tile_res = melt(tile_res)
                        ## In this section we run the 
                        ## CTCF Chip-Seq data base hits 
                                               
                        ctcf_res = data.frame(
                            MOCHA= length(subsetByOverlaps(mocha_gr, ctcf,)),
                            MACS2= length(subsetByOverlaps(macs2_gr, ctcf)),
                            HOMER= length(subsetByOverlaps(homer_gr, ctcf))
                        )
                        ctcf_res$CellPopulation <- x
                        ctcf_res = melt(ctcf_res)   
                        
                        ### Avg Peakset  
                        tss_res = data.table(MOCHA=  length(subsetByOverlaps(mocha_gr, tss)),
                                             MACS2= length(subsetByOverlaps(macs2_gr, tss )),
                                             HOMER= length(subsetByOverlaps(homer_gr, tss))
                        )

                        tss_res$CellPopulation <- x
                        tss_res = melt(tss_res)   

                        results = list(
                            tiles = tile_res,
                            CTCF = ctcf_res,
                            TSS = tss_res,
                            MOCHA = mocha_tile$tileID,
                            MACS2 = macs2_tile_list$tileID,
                            HOMER = homer_peak_list$tileID)
                        setwd(datasetDir)
                        results
                        
                        }
                    )


tile_counts = rbindlist(lapply(total_res,
               function(x)
                   x$tiles
                               )
                        )

ctcf_counts = rbindlist(lapply(total_res,
               function(x)
                   x$CTCF
                               )
                               
                       )

tss_counts = rbindlist(lapply(total_res,
               function(x)
                   x$TSS
                               )
                               
                       )

models=c('HOMER',
              'MACS2',
              'MOCHA')
             


cells_per_sample$Index = 1:20
cells_per_sample$CellPopulation = cells_per_sample$CellType_Organ

tile_counts$Model = factor(tile_counts$variable, levels=models)
tss_counts$Model  = factor(tss_counts$variable, levels=models)
ctcf_counts$Model = factor(ctcf_counts$variable, levels=models)

tile_counts$CellPopulation=factor(tile_counts$CellPopulation,
                            levels=cells_per_sample$CellType_Organ)

tss_counts$CellPopulation=factor(tss_counts$CellPopulation,
                            levels=cells_per_sample$CellType_Organ)

ctcf_counts$CellPopulation=factor(ctcf_counts$CellPopulation,
                            levels=cells_per_sample$CellType_Organ)

tile_counts = dplyr::left_join(tile_counts, cells_per_sample, multiple='first')
tss_counts = dplyr::left_join(tss_counts, cells_per_sample, multiple='first')
ctcf_counts = dplyr::left_join(ctcf_counts, cells_per_sample, multiple='first')

setwd('../results')
pdf('mouse_cells.pdf')    

  p0=ggplot(cells_per_sample,
       aes(x=CellType_Organ,
           y=Count))+geom_bar(stat='identity')+ 
            xlab('Cell and Organ Sample')+ylab('# of Cells')+
            scale_y_continuous(trans='log10')+  
            theme_minimal()+
            theme(legend.position='none',
                  axis.text.x=element_text(size=12, angle=90),
                 axis.text.y=element_text(size=16))+
            scale_fill_cell()
   print(p0)

dev.off()

pdf('mouse_tiles.pdf')    

p = ggplot(tile_counts,
       aes(x=CellPopulation,
           y=value,
           fill=Model))+
        #geom_bar(stat='identity', position='dodge')+
        theme_minimal() + 
        theme(text=element_text(size=16),
            axis.text.x = element_text(size=12,
                                         angle=90))+
    ggtitle('Tile Counts')+
    xlab('Cell & Organ')+
    geom_point(aes(x=Index, y=value, color=Model), stat='identity')+
    geom_smooth(method='lm',
               aes(x=Index, y=value, color=Model), stat='identity')
    print(p)
dev.off()

pdf('mouse_ctcf.pdf')    
p2 = ggplot(ctcf_counts,
       aes(x=CellPopulation,
           y=value,
           fill=Model))+geom_bar(stat='identity',
                                    position='dodge')+
        theme_minimal() + 
       theme(text=element_text(size=16),
            axis.text.x = element_text(size=12,
                                         angle=90))+
    ggtitle('CTCF Counts')+
    xlab('Cell & Organ')
    print(p2)
dev.off()

pdf('mouse_tss.pdf')
p3 = ggplot(tss_counts,
       aes(x=CellPopulation,
           y=value,
           fill=Model))+geom_bar(stat='identity',
                                    position='dodge')+
        theme_minimal() + 
        theme(text=element_text(size=16),
            axis.text.x = element_text(size=12,
                                         angle=90))+
    ggtitle('TSS Counts')+
    xlab('Cell & Organ')
    print(p3)
dev.off()
                        

                        
                     
# Load the ArchR Project
mouseArchR <- loadArchRProject("/home/jupyter/MOCHA_Revision/MouseData")
metadf <- getCellColData(mouseArchR) 
metadt = as.data.table(metadf)
metadt$Organ = sapply(metadt$CellType_Organ,
                      function(x)
                          unlist(strsplit(x, '#'))[2]
                      )

### calculate counts 
tmp = metadt[, list(nFrags=median(nFrags), 
                    cellCount = .N),
             by=list(CellPop,Organ)]
                     
tmp = tmp %>% arrange(nFrags)

pdf('frags_cellOrgan.pdf')
ggplot(tmp,
       aes(y=reorder(CellPop, nFrags, max),
           x=nFrags,
          fill=Organ))+geom_bar(stat='identity', position='dodge')+
        theme_minimal()+
        theme(text=element_text(size=12),
             axis.text = element_text(size=14))+
        xlab('median Fragments per cell')+
        ylab('Cell Population')+
        ggtitle('# Fragments per Cell & Organ')+
        geom_vline(xintercept=7804)
dev.off()

tmp2 = metadt[, list(nFrags=median(nFrags), 
                     Organ = first(Organ),
                     CellPop=first(CellPop),
                    cellCount = .N),
             by=list(CellType_Organ)]


### calculate counts 
tmp = metadt[, list(nFrags=mean(nFrags), 
                    cellCount = .N),
             by=list(CellPop,Organ)]
                     
tmp = tmp %>% arrange(nFrags)

png('frags_cellOrgan2.png', width=700, height=900)
ggplot(tmp,
       aes(y=reorder(CellPop, nFrags, max),
           x=nFrags,
          fill=Organ))+geom_bar(stat='identity', position='dodge')+
        theme_minimal()+
        theme(text=element_text(size=12),
             axis.text = element_text(size=14))+
        xlab('median Fragments per cell')+
        ylab('Cell Population')+
        ggtitle('# Fragments per Cell & Organ')+
        geom_vline(xintercept=7804)
dev.off()

tmp2 = metadt[, list(nFrags=median(nFrags), 
                     Organ = first(Organ),
                     CellPop=first(CellPop),
                    cellCount = .N),
             by=list(CellType_Organ)]



mocha_tiles = lapply(total_res,
               function(x)
                   x$MOCHA
                               )
macs2_tiles = lapply(total_res,
               function(x)
                   x$MACS2
                               )
                                                
homer_tiles = lapply(total_res,
               function(x)
                   x$HOMER
                               )


overlapLengths = lapply(1:20,
                        function(x)
                            c(  length(intersect(mocha_tiles[[x]], 
                                                 macs2_tiles[[x]])),
                                length(intersect(mocha_tiles[[x]],
                                                 homer_tiles[[x]])),                          
                                       
                                length(Reduce(intersect, list(mocha_tiles[[x]],
                                                          macs2_tiles[[x]],
                                                          homer_tiles[[x]]))),
                                length(intersect(macs2_tiles[[x]],
                                                 homer_tiles[[x]]))
                                  )
                        )


df = rbindlist(lapply(1:20,
       function(x) data.table(
                       Intersection_MACS2_MOCHA = overlapLengths[[x]][1],
                       Intersection_HOMER_MOCHA = overlapLengths[[x]][2],  
                       Intersection_HOMER_MACS2 = overlapLengths[[x]][4],  
                       ThreeWay_Intersection = overlapLengths[[x]][3],  
                       MOCHA_Tiles = length(mocha_tiles[[x]]),
                       MACS2_Tiles = length(macs2_tiles[[x]]),
                       HOMER_Tiles = length(homer_tiles[[x]])
       )
       )
      )
df %>% arrange(MOCHA_Tiles)

summary(df$Intersection_MACS2_MOCHA / df$MOCHA_Tiles)
summary(df$Intersection_HOMER_MOCHA / df$MOCHA_Tiles)
summary(df$Intersection_HOMER_MACS2 / df$MACS2_Tiles)
summary(df$Intersection_HOMER_MACS2 / df$HOMER_Tiles)


summary(df$MOCHA_Tiles)
summary(df$MACS2_Tiles)
summary(df$HOMER_Tiles)
