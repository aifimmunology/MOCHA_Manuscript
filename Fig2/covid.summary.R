################################################################
################################################################
### aggregate results 

## load libraries 
require(MOCHA)
require(parallel)
#require(scMACS)
require(ggpubr)
require(data.table)
require(ggplot2)
require(ArchR)
require(MultiAssayExperiment)
require(RaggedExperiment)

## Identify Cell populations
## to analyze 

cells = c('B naive', 'CD16 Mono','CD4 CTL TEM')

### set directory 
homeDir = '/home/jupyter/MOCHA_Manuscript/Fig2/'
setwd(homeDir)
ctcf <- plyranges::read_bed('All_Blood_CTCF_Sites.bed')
load('tss_reorganized.RDS')

source('../theme.R')
source('helper_granges.R')
source('utils.R')


## load ArchR project
ArchRProj = loadArchRProject('/home/jupyter/FullCovid')
metadata = as.data.table(ArchRProj@cellColData)

################################################################
### load MOCHA Tiles 
datasetDir = "/home/jupyter/MOCHA_Manuscript/Fig2/Covid-19"
setwd(datasetDir)
cell_dirs <-  dir('Macs2/')
cell_dirs = cell_dirs[c(1,3,4)]
tileResults <- readRDS('MOCHA.RDS')

################################################################
################################################################
### panelA cell counts 

cells_per_sample <- metadata[
    CellSubsets %in% cells, 
    list(Count= .N), 
    by= list(CellSubsets, Sample)]

cell_levels = c('CD16 Mono','B naive','CD4 CTL TEM')

cells_per_sample$CellSubsets = factor(cells_per_sample$CellSubsets,
                                      levels=cell_levels)

cells_early_visit <- cells_per_sample[Sample %in% colnames(tileResults[[1]])]





sizeText = 17
lwd=2

################################################################
################################################################

### summarize across 
### the selected cell types 
setwd(datasetDir)
total_res <- lapply(1:3, 
                    function(x) summarize_celltypes(x)
                    )

################################################################
################################################################

total_res_filtered = total_res
cell_levels= c('CD4 CTL TEM','B naive','CD16 Mono')

setwd('results')

################################################################
################################################################
### plot tiles per 
### sample per cell type 
tiles <- lapply(total_res_filtered,
                function(x)
                    x$Tiles
                )
tiles <- rbindlist(tiles)
tiles[tiles$variable=='Macs2']$variable <- 'MACS2'
tiles[tiles$variable=='Homer']$variable <- 'HOMER'


cell_levels = c('CD16 Mono','B naive','CD4 CTL TEM')
tiles$CellPop = factor(tiles$CellPop,
                       levels=cell_levels)
tiles$variable = factor(tiles$variable,
                       levels=c('MACS2','HOMER','MOCHA'))

pdf('TileCounts_h.pdf', width=6,height=2.5)
ggplot(tiles,
       aes(x=reorder(variable, value, median),
           y=value,
          fill=variable))+geom_violin(scale='width')+
        scale_fill_MOCHA()+
        facet_wrap(~CellPop, ncol=3)+
    theme_minimal()+ylab('')+xlab('')+
  theme(legend.position='none',
          text=element_text(size=16),
          strip.text.x=element_blank(),
          axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())

dev.off()

### HYPOTHESIS TEST
pairwise_wilcox <- function(groupA, groupB, cell){

    valsA = tiles[variable==groupA & CellPop==cell]$value
    valsB = tiles[variable==groupB & CellPop==cell]$value
    
    res = data.frame(
        A = groupA,
        ValsA= median(valsA),
        B = groupB,
        ValsB=median(valsB),
        P_value=wilcox.test(valsA, valsB)$p.value
    )
    res
}

options(scipen=999)
tile_pvals <- rbind(
    pairwise_wilcox('MOCHA','MACS2','B naive'),
    pairwise_wilcox('MOCHA','MACS2','CD16 Mono'),
    pairwise_wilcox('MOCHA','MACS2','CD4 CTL TEM'),
    pairwise_wilcox('MOCHA','HOMER','B naive'),
    pairwise_wilcox('MOCHA','HOMER','CD16 Mono'),
    pairwise_wilcox('MOCHA','HOMER','CD4 CTL TEM'))


summarized_tiles = tiles[, median(value), by=list(variable, CellPop)]

## range of num peaks 
summarized_tiles[, list(min(V1),max(V1)), by=variable]

## range of performance 
dtiles = dcast(summarized_tiles,  CellPop~variable)
dtiles$Macs2_Homer <- dtiles$MOCHA / dtiles$MACS2
dtiles$MOCHA_Homer <- dtiles$MOCHA / dtiles$HOMER
dtiles$Combined_Ratio <- (dtiles$Macs2_Homer + dtiles$MOCHA_Homer)/2
################################################################
################################################################
### plot tiles per 
### sample per cell type 

metric_levels = c('Tiles','CTCF','TSS','Promoter')

cumul_res <- lapply(total_res_filtered,
                function(x)
                    x$Cumul
                )

cumul_res <- rbindlist(cumul_res)
cumul_res[cumul_res$Model=='Macs2']$Model <- 'MACS2'
cumul_res[cumul_res$Model=='Homer']$Model <- 'HOMER'
cumul_res$Cell = factor(cumul_res$Cell,
                       levels=cell_levels)
cumul_res$Metric <- factor(cumul_res$Metric,
                           levels=metric_levels)

pdf('CumulativeResults_tiles.pdf',
   width=7, height=2.5)
ggplot(cumul_res[Metric == 'Tiles'],
       aes(x=Pct0,
           y=values,
           col=Model))+
            geom_line(linewidth=lwd)+
        facet_wrap( ~ Cell, ncol=3)+
        theme_classic()+xlab('')+ylab('')+
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size=16),
             strip.background = element_blank(),
             strip.text=element_text(size=0),
             legend.position='none')+ylim(10000,250000)+
        scale_col_MOCHA()
dev.off()



pdf('CumulativeResults_ctcf.pdf',
   width=7, height=2.5)
ggplot(cumul_res[Metric == 'CTCF'],
       aes(x=Pct0,
           y=values,
           col=Model))+
            geom_line(linewidth=lwd)+
        facet_wrap( ~ Cell, ncol=3)+
        theme_classic()+xlab('')+ylab('')+
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size=16),
             strip.background = element_blank(),
             strip.text=element_text(size=0),
             legend.position='none')+ylim(10000,120000)+
        scale_col_MOCHA()

dev.off()

pdf('CumulativeResults_tss.pdf',
   width=7, height=2.5)
ggplot(cumul_res[Metric == 'TSS'],
       aes(x=Pct0,
           y=values,
           col=Model))+
            geom_line(linewidth=lwd)+
        facet_wrap( ~ Cell, ncol=3)+
        theme_classic()+xlab('')+ylab('')+
        theme(axis.text.x = element_text(size=16, angle=90),
              axis.text.y = element_text(size=16),
             strip.background = element_blank(),
             strip.text=element_text(size=0),
             legend.position='none')+ylim(4000,25000)+
        scale_col_MOCHA()

dev.off()

### quantify sparsity at 50% 

thresh = 20/39

get_cell_metric_summary <- function(thresh,metric){
    
    df_list <- lapply(cells,
           function(x){
                d1= cumul_res[Pct0 == thresh & Metric == metric & Cell ==x]
                d2= cumul_res[Cell==x & Metric==metric, max(values), by=Model]
                d3 =cbind(d1,d2)
                d3$Pct50 = d3[,2]/d3$V1
                d3$MOCHA_RATIO = 1/(d3$values/d3$values[1])
                max_zeroes <-  cumul_res[Cell==x & Metric==metric, max(Pct0), by=Model]
                increase_tiles = mean(d3$MOCHA_RATIO[2:3])
                data.frame(
                    Increase_tiles = increase_tiles,
                    MOCHA_largest0 = max_zeroes$V1[1],
                    MACS2_largest0 = max_zeroes$V1[2],
                    HOMER_largest0 = max_zeroes$V1[3])
               }
           )
    df_list = rbindlist(df_list)
    print(df_list)
    
}
get_cell_metric_summary(thresh,'Tiles')
tmp <- get_cell_metric_summary(thresh,'CTCF')
tmp2 =get_cell_metric_summary(thresh,'TSS')

################################################################
################################################################
### plot tiles per 
### sample per cell type 

metric_levels = c('Tiles','CTCF','TSS','Promoter')

hist_res <- lapply(total_res_filtered,
                function(x)
                    x$Hist
                )

hist_res <- rbindlist(hist_res)
hist_res[hist_res$Model=='Macs2']$Model <- 'MACS2'
hist_res[hist_res$Model=='Homer']$Model <- 'HOMER'
hist_res$Cell = factor(hist_res$Cell,
                       levels=cell_levels)
hist_res$Metric <- factor(hist_res$Metric,
                           levels=metric_levels)



### generate violin showing
### each sample's cellular 
### abundance for each cell type

cells_per_sample$CellGroup <- 'Monocytes'
cells_per_sample[cells_per_sample$CellSubsets=='B naive']$CellGroup <- 'B'
cells_per_sample[cells_per_sample$CellSubsets=='CD4 CTL TEM']$CellGroup <- 'T'

pdf('cellsCounts_violin_h.pdf', height=2.5,width=7)
    p=ggplot(cells_per_sample[CellSubsets !='CD14 Mono'],
       aes(x=CellSubsets,
           y=Count,
          fill=CellGroup))+geom_violin(linewidth=1.25)+ ThemeMain+
            xlab('')+ylab('')+
            scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                 axis.text.y=element_text(size=16),
                 axis.text.x=element_blank(),
                 title=element_blank())+
            scale_fill_cell()
   print(p)
   
dev.off()

################################################################           
################################################################

# venn diagrams 
venns  <- lapply(1:3,
                function(x){
                    x = total_res_filtered[[x]]
                    res= data.table(
                        Tiles = c(x$VennList[[1]],
                                  x$VennList[[2]],
                                  x$VennList[[3]]),
                        Method= c(rep('MOCHA', length(x$VennList[[1]])),
                                  rep('MACS2', length(x$VennList[[2]])),
                                  rep('HOMER', length(x$VennList[[3]])))
                                  
                        )
                    res$CellPopulation <- x$Hist$Cell[1]
                    res
                    }
                    )

require(UpSetR)

### B naive
pdf(paste(venns[[1]]$CellPopulation[1],'.pdf',sep=''), width=9, height=6)

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
   # main.bar.color=c('black','blue',
   #                  'red','purple','turquoise','orange','green'),
    sets.bar.color=c('blue','red','green'))

print(p)
dev.off()



### CD 16
pdf(paste(venns[[2]]$CellPopulation[1],'.pdf',sep=''), width=9, height=6)
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
#                 ?s    'red','turquoise','purple', 'orange','green'),
    sets.bar.color=c('blue','red','green'))

print(p)
dev.off()


### CD8 TEMS

pdf(paste(venns[[3]]$CellPopulation[1],'.pdf',sep=''), width=9, height=6)
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
#                     'purple','red','turquoise', 'orange','green'),
    sets.bar.color=c('blue','red','green'))

print(p)
dev.off()

###############################################################
# tabular upset plots 

vennResults <- rbindlist(lapply(1:3, 
                 function(x) summarize_upset_plots(x)
                 ))

vennResults$MOCHA_MACS2 = vennResults$MOCHA / vennResults$MACS2
vennResults$MOCHA_HOMER = vennResults$MOCHA / vennResults$HOMER
vennResults$Combined_Ratio <- (vennResults$MOCHA_MACS2 + vennResults$MOCHA_HOMER)/2

# TSAM Increase
mean(vennResults$Combined_Ratio)

# Recall Results 
vennResults$Combined_Percentage <- (vennResults$MACS2_Detected_P + vennResults$HOMER_Detected_C)/2

################################################################           
################################################################

# pie charts 

pieDfs <- lapply(total_res_filtered,
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

# pie charts 

tileTypes_mocha_unique <- lapply(1:3,
                 function(x){
                     y=data.frame(cbind(total_res[[x]]$TileTypesUniqueMOCHA,
                           cells[x]))
                     y$percentage = y$value / sum(y$value)
                     y
                     }
                 )

mocha_unique= rbindlist(tileTypes_mocha_unique)
    

tileTypes_mocha_missed <- lapply(1:3,
                 function(x){
                     y=data.frame(cbind(total_res[[x]]$TileTypesMissedMOCHA,
                           cells[x]))
                     y$percentage = y$value / sum(y$value)
                     y
                     }
                                 )
mocha_missed = rbindlist(tileTypes_mocha_missed)
    
tileTypes_allThree <- lapply(1:3,
                 function(x){
                      y=data.frame(cbind(total_res[[x]]$TileTypesAllThree,
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

## unique region distributions 

summarize_dist <- function(i){
    
    dt= data.table(
        values= c(total_res[[i]]$MochaUnique_lambda,
                  total_res[[i]]$MochaMissed_lambda),
        tile = c(rep('MOCHA Unique', length(total_res[[i]]$MochaUnique_lambda)),
                 rep('Missed by MOCHA',length(total_res[[i]]$MochaMissed_lambda))
                 )
        )
    
    ggplot(dt,
           aes(x=tile,
               y=values,
               fill=tile))+geom_violin()+
            ggpubr::stat_compare_means()+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          text = element_text(size=14))+
    xlab(cells[i])+ylab('Avg Log2 Values')+
    ggtitle("Covid-19 Dataset")
    
    
}

plotlist <- lapply(1:3,
       function(x) summarize_dist(x)
                  )

ggpubr::ggarrange(plotlist=plotlist, ncol=3, common.legend=T)
ggsave(filename = "covid_intensity_comparisons.pdf", height=8, width=12)


### generate intensity plots
intensity_dt = rbindlist(lapply(1:3,
       function(i) 
           data.table(
        values= c(total_res[[i]]$MochaUnique_lambda,
                  total_res[[i]]$MochaMissed_lambda),
        tile = c(rep('MOCHA Unique', length(total_res[[i]]$MochaUnique_lambda)),
                 rep('Missed by MOCHA',length(total_res[[i]]$MochaMissed_lambda))
                 ),
               Cell = cells[i]
        )
       )
                         )




cells_per_sample$Dataset='Covid'
tiles$Dataset='Covid'
vennResults$Dataset='COVID'
intensity_dt$Dataset='COVID'
cumul_res$Dataset ='COVID'

write.csv(cells_per_sample,
          file='cells_covid.csv')
write.csv(tiles,
          file='tiles_covid.csv')
write.csv(vennResults,
          file='venn_covid.csv')
write.csv(intensity_dt,
          file='intensity_covid.csv')
write.csv(cumul_res,
          file='zeroes_covid.csv')          
          
break

##### Save Peaksets with intensities 
summarize_peak_intensities <- function(x){
    
        mochaPeaks = StringsToGRanges(total_res[[x]]$VennList[[1]])
        macs2Peaks = StringsToGRanges(total_res[[x]]$VennList[[2]])
        homerPeaks = StringsToGRanges(total_res[[x]]$VennList[[3]])
    
        fullPeakset <- MOCHA::StringsToGRanges(names(total_res[[x]]$All3Lambda))
        fullPeakset$score <- total_res_filtered[[x]]$All3Lambda
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
       plyranges::write_bed(peak_intensities[[x]], file=paste(cells[x],'intensities_covid.bed',sep='_'))
       )