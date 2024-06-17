###################################################
## hematopoiesis data comparisons 
###################################################

# install.packages("remotes")
# remotes::install_github("wcstcyx/TxDb.Hsapiens.UCSC.hg19.refGene")
setwd('/home/jupyter/MOCHA_Manuscript/Fig2/BoneMarrow')
## load libraries 
require(MOCHA)
# require(ggpubr)
require(UpSetR)
require(data.table)
require(ggplot2)
require(ArchR)
require(MultiAssayExperiment)
require(RaggedExperiment)
source('../../theme.R')
source('../helper_granges.R')
source('../utils.R')
require(parallel)

homeDir = '/home/jupyter/MOCHA_Manuscript2/Fig2/BoneMarrow'
setwd(homeDir)

# # Identify Cell populations
# # to analyze 

cells <- c('10_cDC','20_CD4.N1','12_CD14.Mono.2')
cell_dirs <-  dir('Macs2/')
cell_dirs = cell_dirs[c(1,4,2)]


## Load CTCF & TSS Datasets
ctcf <- plyranges::read_bed('../All_Blood_CTCF_hg19.bed')
load('../TSS_HG19.RDS')

## load ArchR project
ArchRProj = loadArchRProject('/home/jupyter/DoubletFreeBoneMarrow')
metadata = as.data.table(ArchRProj@cellColData)

### load databases 
#BiocManager::install('BSgenome.Hsapiens.UCSC.hg19') 
#remotes::install_github("wcstcyx/TxDb.Hsapiens.UCSC.hg19.refGene")

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.refGene)

library(org.Hs.eg.db)

TxDb <- TxDb.Hsapiens.UCSC.hg19.refGene
Org <- org.Hs.eg.db


################################################################
### load MOCHA Tiles 
tileResults <- readRDS('MOCHA.RDS')
datasetDir='/home/jupyter/MOCHA_Manuscript2/Fig2/BoneMarrow/results'
setwd(datasetDir)
# ###############################################################
# ###############################################################
# ## panelA 

cells_per_sample <- metadata[
    predictedGroup %in% cells, 
    list(Count= .N), 
    by= list(predictedGroup)]

cell_levels = c('10_cDC','20_CD4.N1', '12_CD14.Mono.2')

cells_per_sample$CellSubsets = factor(cells_per_sample$predictedGroup,
                                      levels=cell_levels)


# ###############################################################
# ###############################################################

extract_tiles <- function(i){
    setwd(homeDir)

    #### Extract MOCHA tiles 
    subsetTileResults <- tileResults[[cells[i]]]
    tsam <- RaggedExperiment::compactAssay(subsetTileResults, i='TotalIntensity')
    tsam = data.frame(tsam[, which(colnames(tsam)==cells[i])])

    ### 
    subsetTileResults <- subsetTileResults@assays[[grep(cells[i], 
                                                        names(subsetTileResults@assays))]]
                                                 
    
    MOCHA_tileList <- list(subsetTileResults)

    sampleTileMatrix <- as.data.table(subsetTileResults)

    ################################################################
    ################################################################
    ### 

    #### Extract Macs2 Tiles 
    setwd('Macs2')
    setwd(cell_dirs[i])

    ### Macs2 Peaks
    fnames <- dir()
    fnames <- fnames[grep('broadPeak', fnames)]

    ## read summits list
    macs2_tile_list <- fread(fnames, data.table=T)
                
    ### convert broad peaks into accessible tiles
    macs2_tile_list2 <- peak_to_tile_macs2_bulk(1, list(macs2_tile_list))
    macs2_tile_list2 = macs2_tile_list2[tileID %in% row.names(tsam)]
    ################################################################
    ################################################################


    ################################################################
    ################################################################
    #### Extract Homer Tiles 
    setwd('../../Homer')
    setwd(cell_dirs[i])

    ### Homer Peaks
    fnames <- dir()
    peakNames <- paste(fnames,'/peaks.txt', sep='')
    peakNames <- peakNames[!peakNames %in% c("ArchRLogs/peaks.txt",'-style/peaks.txt')]

    ### extract tiles 
    homer_peak_list <- fread(peakNames)

    ### Transform to tiles 
    homer_peak_list2 <- peak_to_tile_homer_bulk(1, list(homer_peak_list))
    homer_peak_list2 = homer_peak_list2[tileID %in% row.names(tsam)]
                                  
    ################################################################
    ################################################################
    
    ### convert objects to Granges
    mocha_gr <- makeGRangesFromDataFrame(MOCHA_tileList[[1]], 
                                         keep.extra.columns=T)
    
    mocha_gr = mocha_gr[mocha_gr$peak==T]
    
    macs2_gr <- makeGRangesFromDataFrame(macs2_tile_list2, 
                                         keep.extra.columns=T)
    
    homer_gr <- makeGRangesFromDataFrame(homer_peak_list2, 
                                         keep.extra.columns=T)    

    N =  ncol(sampleTileMatrix)

    tilesPerMethod<- data.table(
                MACS2 = nrow(macs2_tile_list2),
                HOMER = nrow(homer_peak_list2),
                MOCHA = length(mocha_gr)
            )
    tilesPerMethod$CellPopulation = cells[i]

    tilesPerMethod =  melt(tilesPerMethod)
    
    
    ## In this section we run the 
    ## CTCF Chip-Seq data base hits 
    ctcf_res = data.frame(
              MOCHA= length(subsetByOverlaps(mocha_gr, ctcf,)),
              MACS2= length(subsetByOverlaps(macs2_gr, ctcf)),
              HOMER= length(subsetByOverlaps(homer_gr, ctcf))
        )
    ctcf_res$CellPopulation <- cells[i]   
    ctcf_res = melt(ctcf_res)       
    
    ################################################################   
    ################################################################
        
    ### to calculate tss overlap 

    ### Avg Peakset  
    tss_res = data.table(MOCHA=  length(subsetByOverlaps(mocha_gr, tss_hg19)),
                          MACS2= length(subsetByOverlaps(macs2_gr, tss_hg19 )),
                          HOMER= length(subsetByOverlaps(homer_gr,tss_hg19))
                                        )
    
    tss_res$CellPopulation <- cells[i]   
    tss_res = melt(tss_res)
    
    ################################################################
    ################################################################
    
    
    macs2_gr <- makeGRangesFromDataFrame(macs2_tile_list2, 
                                         keep.extra.columns=T)
    
    homer_gr <- makeGRangesFromDataFrame(homer_peak_list2, 
                                         keep.extra.columns=T)    
    
    
    ### create annotations 
    mocha_gr <- annotateTiles(mocha_gr, TxDb=TxDb, Org=Org)
    macs2_gr <- annotateTiles(macs2_gr, TxDb=TxDb, Org=Org)
    homer_gr <- annotateTiles(homer_gr, TxDb=TxDb, Org=Org)
    
    macs2_homer_gr <- unique(c(macs2_gr, homer_gr))
    
    promoter_res = data.frame(
        MOCHA= sum(mocha_gr$tileType=='Promoter'),
        MACS2= sum(macs2_gr$tileType=='Promoter'),
        HOMER= sum(homer_gr$tileType=='Promoter'))
    
    promoter_res$CellPopulation <- cells[i]   
    promoter_res = melt(promoter_res)
    
    ### Percentage by Tiles 
    distributions_tiles <- data.frame(rbind(table(mocha_gr$tileType),
          table(macs2_gr$tileType),
          table(homer_gr$tileType)))
    
    distributions_tiles$variable=c('MOCHA','MACS2','HOMER')
    distributions_tiles$CellPopulation <- cells[i]   
    distributions_tiles = melt(distributions_tiles)
    colnames(distributions_tiles)[3] <- 'tileType'
    
    
    mocha_unique_idx <- which(!mocha_gr$tileID %in% union(macs2_gr$tileID,
                                                         homer_gr$tileID))
    mocha_missed_idx <- which(!macs2_homer_gr$tileID %in%
                             mocha_gr$tileID)
    common_all_three <- which(mocha_gr$tileID %in% intersect(macs2_gr$tileID,
                                                         homer_gr$tileID))   
    
    mocha_unique_tileType =  as.data.table(mocha_gr[mocha_unique_idx])
    mocha_missed_tileType =  as.data.table(macs2_homer_gr[mocha_missed_idx])
    common_three_tileType =  as.data.table(mocha_gr[common_all_three])
    
    ### Get Counts 
    mocha_unique_tileType =  mocha_unique_tileType[, list(value=.N), by=tileType]
    mocha_missed_tileType =  mocha_missed_tileType[, list(value=.N), by=tileType]
    common_three_tileType =  common_three_tileType[, list(value=.N), by=tileType]
         
    ### unique lambda
    tsam = log2(tsam+1)
    unique_lambda= tsam[row.names(tsam) %in% mocha_gr$tileID[mocha_unique_idx],]
    
    ### missed lambda
    missedLambda= tsam[row.names(tsam) %in% macs2_homer_gr$tileID[mocha_missed_idx],]                  
    ### all three 
    #all3Lambda= tsam
    #names(all3Lambda) = row.names(tsam)
     
    ################################################################
    ################################################################
    
    ## summarize results by 
    ## different groups 
    
    results <- list(
        Tiles=tilesPerMethod,
        TSS= tss_res,
        CTCF=ctcf_res,
        Promoter = promoter_res,
        TileDistributions = distributions_tiles,
        MOCHA_tiles = mocha_gr$tileID,
        MACS2_tiles = macs2_gr$tileID,
        Homer_tiles = homer_gr$tileID,
        TileTypesUniqueMOCHA=mocha_unique_tileType,
                TileTypesAllThree = common_three_tileType,
                TileTypesMissedMOCHA=mocha_missed_tileType,
                MochaUnique_lambda=unique_lambda,
                MochaMissed_lambda= missedLambda

    )
    
    setwd(homeDir)
    return(results)

}

################################################################
################################################################
tileCountsPerMethod_list = lapply(1:3,
                           function(x) extract_tiles(x)
                )

model_levels <- c('MACS2','HOMER','MOCHA')

setwd('/home/jupyter/MOCHA_Manuscript2/Fig2/BoneMarrow/results')

#### plot results by method 
tiles <- rbindlist(lapply(tileCountsPerMethod_list,
                function(x)
                    x$Tiles
                          )
                          )
tiles$CellPopulation = factor(tiles$CellPopulation,
                            levels=cell_levels)
tiles$CellPopulation = gsub('20_|10_|12_','',tiles$CellPopulation)
tiles$CellPopulation = factor(tiles$CellPopulation,
                              levels=c('cDC','CD4.N1','CD14.Mono.2'))

tiles$variable <- factor(tiles$variable,
                       levels=model_levels)
                   
                          
                          
pdf('tiles.pdf', height=3.5, width=6,)
ggplot(tiles,
       aes(x=variable,
           y=value,
           fill=variable))+
                          geom_bar(stat='identity', position='dodge')+
  scale_fill_MOCHA()+theme_minimal()+
                        theme(text=element_text(size=20),
                              strip.text.x=element_blank(),
                              axis.text.x=element_blank(),
                              legend.position='none')+
                              xlab('')+ylab('')+
    ggtitle('')+facet_wrap(~CellPopulation, ncol=3)

      
dev.off()           

### 
dtiles = dcast(tiles,  CellPopulation~variable)
dtiles$Macs2_Homer <- dtiles$MOCHA / dtiles$MACS2
dtiles$MOCHA_Homer <- dtiles$MOCHA / dtiles$HOMER
dtiles$Combined_Ratio <- (dtiles$Macs2_Homer + dtiles$MOCHA_Homer)/2


################################################################
################################################################                       

#### plot results by method 
tss <- rbindlist(lapply(tileCountsPerMethod_list,
                function(x)
                    x$TSS
                        )
                )

tss$CellPopulation = factor(tss$CellPopulation,
                            levels=cell_levels)
tss$CellPopulation = gsub('20_|10_|12_','',tss$CellPopulation)
tss$variable <- factor(tss$variable,
                       levels=model_levels)
                          
pdf('tss.pdf',
   width=6,height=7)
tss_plot <- ggplot(tss,
       aes(x=reorder(CellPopulation, value, median),
           y=value,
           fill=variable))+
                          geom_bar(stat='identity', position='dodge')+
                        scale_fill_MOCHA()+theme_minimal()+
                        theme(text=element_text(size=22),
                              legend.position='none',
                             axis.text.x=element_blank())+
                              xlab('')+ylab('')+
            ggtitle('')
print(tss_plot)
dev.off()             
          
################################################################                    ################################################################

ctcf_res <- rbindlist(lapply(tileCountsPerMethod_list,
                function(x)
                    x$CTCF
                         )       
                 )
ctcf_res$CellPopulation = factor(ctcf_res$CellPopulation,
                            levels=cell_levels)
ctcf_res$CellPopulation = gsub('20_|10_|12_','',ctcf_res$CellPopulation)
ctcf_res$variable <- factor(ctcf_res$variable,
                       levels=model_levels)
                        
pdf('ctcf.pdf',
   width=6, height=7)
                        
ctcf_plot <- ggplot(ctcf_res,
       aes(x=reorder(CellPopulation, value, median),
           y=value,
           fill=variable))+
        geom_bar(stat='identity', position='dodge')+
        scale_fill_MOCHA()+theme_minimal()+
                        theme(text=element_text(size=22),
                              legend.position='none',
                             axis.text.x=element_blank())+
                              xlab('')+ylab('')+
            ggtitle('')
print(ctcf_plot)
dev.off()      
         
summarize_ctcf_tss <- function(tss){
    tmp = tss
    tmp = tmp[, first(value)/value,by=list(CellPopulation)]
    tmp = tmp[-c(1,4,7),]
    tmp = tmp[, mean(V1), by=CellPopulation]
    c(1-tmp$V1[1],  mean(tmp$V1[2:3]))
}
summarize_ctcf_tss(tss) 
summarize_ctcf_tss(ctcf_res) 

# ################################################################                    ################################################################
# promoter <- rbindlist(lapply(tileCountsPerMethod_list,
#                 function(x)
#                     x$Promoter
#                              )
#                      )
# promoter$CellPopulation = factor(promoter$CellPopulation,
#                             levels=cell_levels)
# promoter$variable <- factor(promoter$variable,
#                        levels=model_levels)
                          

# promoter_plot <- ggplot(promoter,
#        aes(x=reorder(CellPopulation, value, median),
#            y=value,
#            fill=variable))+
#                           geom_bar(stat='identity', position='dodge')+
#      scale_fill_MOCHA()+
#                         theme(text=element_text(size=15),
#                              axis.text.x=element_text(size=15, angle=90))+
#                               xlab('')+
#         ggtitle('Promoters')
      
          

# ################################################################                    ################################################################
# distals <- rbindlist(lapply(tileCountsPerMethod_list,
#                 function(x)
#                     x$TileDistributions
#                             )
#                     )

# distals$CellPopulation = factor(distals$CellPopulation,
#                             levels=cell_levels)
# distals$variable <- factor(distals$variable,
#                        levels=model_levels)
                          
# pdf('distal.pdf')
# ggplot(distals[variable=='MOCHA'],
#        aes(x=CellPopulation,
#            y=value,
#            fill=tileType))+
#                           geom_bar(stat='identity', position='stack')+
#                         theme(text=element_text(size=15),
#                              axis.text.x=element_text(size=15, angle=90))+
#                               xlab('')+
#     ggtitle('MOCHA Tile Distributions')
      
# dev.off()             
          

################################################################           
################################################################

# venn diagrams 
venns  <- lapply(tileCountsPerMethod_list,
                function(x){
                    
                    res= data.table(
                        Tiles = c(x[[6]],
                                  x[[7]],
                                  x[[8]]),
                        Method= c(rep('MOCHA', length(x[[6]])),
                                  rep('MACS2', length(x[[7]])),
                                  rep('HOMER', length(x[[8]])))
                                  
                        )
                    res$CellPopulation <- x$TSS$CellPopulation[1]
                    res
                    }
                    )

require(UpSetR)
### cDC UpSet
pdf(paste(cells[1],'.pdf',sep=''), width=8,height=6)
listInput <- list(
        "MOCHA"=venns[[1]][Method=='MOCHA']$Tiles,
        "MACS2"=venns[[1]][Method=='MACS2']$Tiles,
        "HOMER"=venns[[1]][Method=='HOMER']$Tiles
)
       
p <- upset(fromList(listInput), order.by = "freq",
   point.size = 3.5, line.size = 2, 
mainbar.y.label = "Common", sets.x.label = "Method", 
text.scale = c(0, 2.2, 2.2, 2.2,2.2 , 2.2),
                      keep.order=T,
                      sets=c('MOCHA','HOMER','MACS2'),
#    main.bar.color=c('black','blue',
#                     'red','turquoise','purple', 'orange','green'),
    sets.bar.color=c('blue','red','green')
   )
print(p)
dev.off()


### CD Naives 
pdf(paste(cells[2],'.pdf',sep=''), width=8, height=6)
listInput <- list(
        "MOCHA"=venns[[2]][Method=='MOCHA']$Tiles,
        "MACS2"=venns[[2]][Method=='MACS2']$Tiles,
        "HOMER"=venns[[2]][Method=='HOMER']$Tiles
)
       
p <- upset(fromList(listInput), order.by = "freq",
   point.size = 3.5, line.size = 2, 
mainbar.y.label = "Common", sets.x.label = "Method", 
text.scale = c(0, 2.2, 2.2, 2.2,2.2 , 2.2),
                      keep.order=T,
                      sets=c('MOCHA','HOMER','MACS2'),
#    main.bar.color=c('black','blue',
#                     'red','turquoise','purple', 'orange','green'),
    sets.bar.color=c('blue','red','green'))
print(p)
dev.off()



### CD14 Monocytes 
pdf(paste(cells[3],'.pdf',sep=''), width=8, height=6)
listInput <- list(
        "MOCHA"=venns[[3]][Method=='MOCHA']$Tiles,
        "MACS2"=venns[[3]][Method=='MACS2']$Tiles,
        "HOMER"=venns[[3]][Method=='HOMER']$Tiles
)
       
p <- upset(fromList(listInput), order.by = "freq",
   point.size = 3.5, line.size = 2, 
mainbar.y.label = "Common", sets.x.label = "Method", 
text.scale = c(0, 2.2, 2.2, 2.2,2.2 , 2.2),
                      keep.order=T,
                      sets=c('MOCHA','HOMER','MACS2'),
#    main.bar.color=c('black','blue',
#                     'red','turquoise','purple', 'orange','green'),
    sets.bar.color=c('blue','red','green'))
print(p)
dev.off()



                    
################################################################
################################################################

cells_per_sample$CellGroup <- 'Monocytes'
cells_per_sample[cells_per_sample$CellSubsets=='10_cDC']$CellGroup <- 'Myeloid'
cells_per_sample[cells_per_sample$CellSubsets=='20_CD4.N1']$CellGroup <- 'T'
cells_per_sample$CellSubsets = gsub('20_|10_|12_','',cells_per_sample$CellSubsets)                    
png('BoneMarrow_cellCounts_h.pdf', height=2.5, units='in',width=7, res=600)

    ggplot(cells_per_sample,
       aes(y=Count,
           x=reorder(CellSubsets,  Count, mean),
           fill=CellGroup
           ))+geom_bar(stat='identity')+ ThemeMain+
            xlab('')+ylab('')+
            scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                 axis.text.y=element_text(size=16),
                 axis.text.x=element_blank(),
                 title=element_text(size=20))+
        scale_fill_cell()
   
dev.off()


################################################################
################################################################
vennResults <- rbindlist(lapply(1:3, 
                 function(x) summarize_upset_plots(x)
                 ))

vennResults

vennResults$MOCHA_MACS2 = vennResults$MOCHA / vennResults$MACS2
vennResults$MOCHA_HOMER = vennResults$MOCHA / vennResults$HOMER
vennResults$Combined_Ratio <- (vennResults$MOCHA_MACS2 + vennResults$MOCHA_HOMER)/2

# TSAM Increase
mean(vennResults$Combined_Ratio)

# Recall Results 

vennResults$Combined_Percentage <- (vennResults$MACS2_Detected_P + vennResults$HOMER_Detected_C)/2



################################################################
################################################################
pieDfs <- lapply(tileCountsPerMethod_list,
                 function(x)
                     x$TileTypesUniqueMOCHA
                 )
    
    
blank_theme <- theme_minimal()+
           theme(
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.border = element_blank(),
           panel.grid=element_blank(),
           axis.ticks = element_blank(),
           plot.title=element_text(size=14, face="bold")
           )
    

draw_pie <- function(x){
    
    df = pieDfs[[x]]
    totalN = sum(df$value)
    
    pdf(paste(cells[x], 'pie.pdf', sep='_'))

         p=ggplot(df, aes(x = "", y = value, fill = tileType)) +
               geom_col(color = "black") +
               geom_text(aes(label = paste(round(value/totalN,2)*100,'%',sep='')),
                         position = position_stack(vjust = 0.5),
                        size=16) +
               coord_polar(theta = "y") +blank_theme+
              theme(legend.position='none',
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank())+
                 scale_fill_brewer()
         print(p)
     dev.off()
    
}

lapply(1:3, function(x) draw_pie(x)
       )


################################################################           
################################################################
break
# pie charts 

# pie charts 

tileTypes_mocha_unique <- lapply(1:3,
                 function(x){
                     y=data.frame(cbind(tileCountsPerMethod_list[[x]]$TileTypesUniqueMOCHA,
                           cells[x]))
                     y$percentage = y$value / sum(y$value)
                     y
                     }
                 )

mocha_unique= rbindlist(tileTypes_mocha_unique)
    

tileTypes_mocha_missed <- lapply(1:3,
                 function(x){
                     y=data.frame(cbind(tileCountsPerMethod_list[[x]]$TileTypesMissedMOCHA,
                           cells[x]))
                     y$percentage = y$value / sum(y$value)
                     y
                     }
                                 )
mocha_missed = rbindlist(tileTypes_mocha_missed)
    
tileTypes_allThree <- lapply(1:3,
                 function(x){
                      y=data.frame(cbind(tileCountsPerMethod_list[[x]]$TileTypesAllThree,
                           cells[x]))
                     y$percentage = y$value / sum(y$value)
                     y
                     }                                    
                 )
allThree = rbindlist(tileTypes_allThree)

    
mocha_unique[tileType=='Promoter', range(percentage),]
mocha_missed[tileType=='Promoter', range(percentage),]
allThree[tileType=='Promoter', range(percentage),]


mocha_unique[tileType=='Promoter', range(value),]
mocha_missed[tileType=='Promoter', range(value),]
allThree[tileType=='Promoter', range(value),]

################################################################           
################################################################

## unique region distributions 

summarize_dist <- function(i){
    
    dt= data.table(
        values= c(tileCountsPerMethod_list[[i]]$MochaUnique_lambda,
                  tileCountsPerMethod_list[[i]]$MochaMissed_lambda),
        tile = c(rep('MOCHA Unique', length(tileCountsPerMethod_list[[i]]$MochaUnique_lambda)),
                 rep('Missed by MOCHA',length(tileCountsPerMethod_list[[i]]$MochaMissed_lambda))
                 )
        )
    
    ggplot(dt,
           aes(x=tile,
               y=values,
               fill=tile))+geom_boxplot(scale='width')+
            ggpubr::stat_compare_means()+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          text = element_text(size=14))+
    xlab(cells[i])+ylab('Avg Log2 Intensity')+
       ggtitle("Bone Marrow Dataset")
    
    
}

plotlist <- lapply(1:3,
       function(x) summarize_dist(x)
                  )

ggpubr::ggarrange(plotlist=plotlist, ncol=3, common.legend=T)
ggsave(filename = "bm_intensity_comparisons.pdf", height=8, width=12)



### Save results to File
intensity_dt = rbindlist(lapply(1:3,
       function(i) 
           data.table(
        values= c(tileCountsPerMethod_list[[i]]$MochaUnique_lambda,
                  tileCountsPerMethod_list[[i]]$MochaMissed_lambda),
        tile = c(rep('MOCHA Unique', length(tileCountsPerMethod_list[[i]]$MochaUnique_lambda)),
                 rep('Missed by MOCHA',length(tileCountsPerMethod_list[[i]]$MochaMissed_lambda))
                 ),
               Cell = cells[i]
        )
       )
                         )



cells_per_sample$Dataset='BoneMarrow'
tiles$Dataset='BoneMarrow'
vennResults$Dataset='BoneMarrow'
intensity_dt$Dataset='BoneMarrow'
ctcf_res$Dataset ='BoneMarrow'
ctcf_res$TileType='CTCF'
tss$Dataset ='BoneMarrow'
tss$TileType='TSS'


write.csv(cells_per_sample,
          file='cells_bm.csv')
write.csv(tiles,
          file='tiles_bm.csv')
write.csv(vennResults,
          file='venn_bm.csv')
write.csv(intensity_dt,
          file='intensity_bm.csv')
write.csv(rbind(tss, ctcf_res),
          file='tss_ctcf.csv')          
          
##### Save Peaksets with intensities 
summarize_peak_intensities <- function(x){
    
        mochaPeaks = StringsToGRanges(tileCountsPerMethod_list[[x]]$MOCHA_tiles)
        macs2Peaks = StringsToGRanges(tileCountsPerMethod_list[[x]]$MACS2_tiles)
        homerPeaks = StringsToGRanges(tileCountsPerMethod_list[[x]]$Homer_tiles)
    
        fullPeakset <- MOCHA::StringsToGRanges(names(tileCountsPerMethod_list[[x]]$All3Lambda))
        fullPeakset$score <- tileCountsPerMethod_list[[x]]$All3Lambda
        fullPeakset$name = ''
        
        mocha_idx <-  queryHits(findOverlaps(fullPeakset, mochaPeaks))
        macs2_idx <-  queryHits(findOverlaps(fullPeakset, macs2Peaks))    
        homer_idx <-  queryHits(findOverlaps(fullPeakset, homerPeaks))  
    
        mocha_and_macs2 = intersect(mocha_idx, macs2_idx)
        mocha_and_homer = intersect(mocha_idx, homer_idx)
        macs2_and_homer = intersect(macs2_idx, homer_idx)
        common = intersect(mocha_and_homer, mocha_and_macs2)
        union_peaks = unique(c(mocha_idx, macs2_idx,homer_idx))
        
        fullPeakset$name[union_peaks] <- 'Common'
    
        fullPeakset$name[setdiff(homer_idx, mocha_and_macs2)]  <- 'HOMER_Unique'
        fullPeakset$name[setdiff(mocha_idx, macs2_and_homer)]  <- 'MOCHA_Unique'    
        fullPeakset$name[setdiff(macs2_idx, mocha_and_homer)]  <- 'MACS2_Unique'       
    
        fullPeakset$name[setdiff(mocha_and_macs2, common)]  <- 'MOCHA_MACS2'
        fullPeakset$name[setdiff(mocha_and_homer, common)]  <- 'MOCHA_HOMER'    
        fullPeakset$name[setdiff(macs2_and_homer, common)]  <- 'MACS2_HOMER'        

        table(fullPeakset$name)
        return(fullPeakset[fullPeakset$name!=''])
}

peak_intensities <- mclapply(1:3,
       function(x)
       summarize_peak_intensities(x),
                         mc.cores=3
       )



lapply(1:3,
       function(x)
       plyranges::write_bed(peak_intensities[[x]], file=paste(cells[x],'intensities_bm.bed',sep='_'))
       )
