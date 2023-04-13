# ##############################################################

## load libraries 
require(MOCHA)
#require(scMACS)
require(ggpubr)
require(parallel)
require(data.table)
require(ggplot2)
require(ArchR)
require(MultiAssayExperiment)
require(RaggedExperiment)

# # Identify Cell populations
# # to analyze 

cells = c('B naive', 'CD16 Mono','CD8 TEM')

### set directory 
homeDir = '/home/jupyter/MOCHA_Manuscript/Fig2/'
setwd(homeDir)

### Load MOCHA tiles,
### TSS and
### CTCF 
ctcf <- plyranges::read_bed('All_Blood_CTCF_hg38.bed')
load('TSS_HG38.RDS')
tileResults <- readRDS('LongPilot/MOCHA.RDS')

### Load Helper Functions 
source('../theme.R')
source('helper_granges.R')
source('utils.R')

## load ArchR project
ArchRProj = loadArchRProject('/home/jupyter/longPilot')
metadata = as.data.table(ArchRProj@cellColData)

summarize_cell <- function(cell){
    
    
    #####################################################################
    ### Load MOCHA tiles
    setwd('/home/jupyter/longPilotValidation/ManuscriptData')
    setwd(cell)
    setwd('bulk')
    setwd('MOCHA')
    
    fnames <- dir()
    tile_list <- mclapply(fnames,
                          function(x) readRDS(x),
                          mc.cores=10
                          )
    
    
    cells_mocha = sapply(tile_list,
                   function(x) 
                       x$numCells[1]
                   )
    tile_list2 <- mclapply(tile_list,
                           function(x)
                               x[peak==T],
                           mc.cores=10
                           )
    
    subsetTileResults <- tileResults[[cell]]
    tsam <- RaggedExperiment::compactAssay(tileResults[[cell]], 
                                           i='TotalIntensity')
    colnames(tsam) <- gsub('_','-', colnames(tsam))
    rm(tile_list)
    
    #####################################################################        
    ### load macs2 tiles 
    setwd('../downsample_macs2_peaks/')
    fnames <- dir()
    fnames <- fnames[grep('broadPeak', fnames)]
    
    ### load macs2 peaks 
    macs2_tile_list <- mclapply(fnames,
                                function(x) fread(x, data.table=T),
                                mc.cores=10
                                )
    
    ### overlay peaks into tiles 
    macs2_tile_list2 <- mclapply(1:length(macs2_tile_list),

                                  function(x) 
                                     peak_to_tile_macs2_bulk(x, 
                                                         macs2_tile_list),
                                  mc.cores=10

                                  )
    
    ### Remove tiles that extended
    ### into empty regions 
    macs2_tile_list3 <- mclapply(macs2_tile_list2,
                                 function(x)
                                     
                                     x[tileID %in% row.names(tsam)],
                                 mc.cores=10
                                 )
    
    ### get cell counts macs2
    cells_macs2 <- sapply(fnames,
                          function(x)
                              as.numeric(strsplit(x,'_')[[1]][1])
                          )    
    rm(macs2_tile_list, macs2_tile_list2)
    #####################################################################        
    ### load homer tiles     
    setwd('../downsample_homer_peaks/homer_dirs')
    fnames <- dir()
    peakNames <- paste(fnames,'/peaks.txt', sep='')
    peakNames <- peakNames[!peakNames %in% c("ArchRLogs/peaks.txt",'-style/peaks.txt')]

    ### extract tiles 
    homer_peak_list <- mclapply(peakNames,
                                        function(x)
                                            fread(x),
                                        mc.cores=10
                                        )

    ### Transform to tiles 
    homer_peak_list2 <- mclapply(1:length(homer_peak_list),
                                  function(x) peak_to_tile_homer_bulk(x, 
                                                                      homer_peak_list),
                                  mc.cores=15
                                  )
    
  
    ### Remove tiles that extended
    ### into empty regions 
    homer_peak_list3 <- mclapply(homer_peak_list2,
                                 function(x)
                                     
                                     x[tileID %in% row.names(tsam)],
                                 mc.cores=10
                                 )
    
    
    ### get cell counts macs2
    cells_homer <- sapply(fnames,
                          function(x)
                              as.numeric(strsplit(x,'_')[[1]][1])
                          )   
    rm(homer_peak_list, homer_peak_list2)
    #####################################################################
    ### Tiles per methods
    
    mocha_tiles = data.table(
        CellCounts = cells_mocha,
        Tiles = sapply(tile_list2, function(x) sum(x$peak)
                       ),
        Model='MOCHA'
        )
    
    macs2_tiles = data.table(
        CellCounts = cells_macs2,
        Tiles = sapply(macs2_tile_list3, function(x) nrow(x)
                       ),
        Model='MACS2'
        )        
    
    homer_tiles = data.table(
        CellCounts = cells_homer,
        Tiles = sapply(homer_peak_list3, function(x) nrow(x)
                       ),
        Model='HOMER'
        
        ) 
    
    tiles_df = rbind(mocha_tiles,
                     macs2_tiles,
                     homer_tiles)
    
    
    #####################################################################
    ### CTCF per methods     
    
    count_overlaps <- function(mat, ctcf){
        
        gr = makeGRangesFromDataFrame(mat)
        overlaps = subsetByOverlaps(gr, ctcf)
        length(overlaps)
        
    }
    
    mocha_ctcf = mclapply(tile_list2,
                          function(x)
                              count_overlaps(x, ctcf),
                          mc.cores=10
                          )

    macs2_ctcf = mclapply(macs2_tile_list3,
                          function(x)
                              count_overlaps(x, ctcf),
                          mc.cores=10
                          )
    
    homer_ctcf = mclapply(homer_peak_list3,
                          function(x)
                              count_overlaps(x, ctcf),
                          mc.cores=10
                          )   
    
    mocha_ctcf = data.table(
        CellCounts = cells_mocha,
        Tiles = mocha_ctcf,
        Model='MOCHA'
        )
    
    macs2_ctcf = data.table(
        CellCounts = cells_macs2,
        Tiles = macs2_ctcf,
        Model='MACS2'
        )        
    
    homer_ctcf = data.table(
        CellCounts = cells_homer,
        Tiles = homer_ctcf,
        Model='HOMER'
        ) 
    
    ctcf_df = rbind(mocha_ctcf,
                     macs2_ctcf,
                     homer_ctcf)
    
    
    #####################################################################
    ### TSS per methods     
        
    mocha_tss = mclapply(tile_list2,
                          function(x)
                              count_overlaps(x, tss_hg38),
                          mc.cores=10
                          )

    macs2_tss = mclapply(macs2_tile_list3,
                          function(x)
                              count_overlaps(x, tss_hg38),
                          mc.cores=10
                          )
    
    homer_tss = mclapply(homer_peak_list3,
                          function(x)
                              count_overlaps(x, tss_hg38),
                          mc.cores=10
                          )   
    
    mocha_tss = data.table(
        CellCounts = cells_mocha,
        Tiles = mocha_tss,
        Model='MOCHA'
        )
    
    macs2_tss = data.table(
        CellCounts = cells_macs2,
        Tiles = macs2_tss,
        Model='MACS2'
        )        
    
    homer_tss = data.table(
        CellCounts = cells_homer,
        Tiles = homer_tss,
        Model='HOMER'
        ) 
    
    tss_df = rbind(mocha_tss,
                     macs2_tss,
                     homer_tss)
    
    tiles_df$Cell = cell
    ctcf_df$Cell =cell
    tss_df$Cell = cell
    
    res_list <- list(Tiles=tiles_df,
                     CTCF=ctcf_df,
                     TSS=tss_df)
    return(res_list)
}

### Extract results from down sampling 
results_list <- lapply(cells,
                       function(x) summarize_cell(x)
                       )

### get tile counts 
results_tiles <- rbindlist(lapply(results_list,
                        function(x)
                            x$Tiles
                        ))

setwd('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling')

a = ggplot(results_tiles,
       aes(x=CellCounts,
           y=Tiles,
           col=Model))+geom_point() +
            geom_line()+
        theme_minimal()+
    facet_wrap(~Cell, ncol=3, scales='free')
Cell, ncol=3)+ylab('Tiles')

### get CTCF counts 
results_ctcf <- rbindlist(lapply(results_list,
                       function(x)
                           x$CTCF
                       )
                          )
results_ctcf$Tiles <- unlist(results_ctcf$Tiles)

b = ggplot(results_ctcf,
       aes(x=CellCounts,
           y=Tiles,
           col=Model))+geom_point() +
            geom_line()+
        theme_minimal()+
        facet_wrap(~Cell, ncol=3, scales='free')+
        ylab('CTCF Counts')
### get TSS counts 
results_tss <- rbindlist(lapply(results_list,
                      function(x)
                          x$TSS
                      )
                         )
results_tss$Tiles <- unlist(results_tss$Tiles)


c = ggplot(results_tss,
       aes(x=CellCounts,
           y=Tiles,
           col=Model))+geom_point() +
            geom_line()+
        theme_minimal()+
        facet_wrap(~Cell, ncol=3, scales='free')+
            ylab('TSS Counts')


ggpubr::ggarrange(a,b,c, ncol=1,
                 common.legend=T)
ggsave('lp_downsampling.pdf', height=9, width=9)
       
               
saveRDS(results_list, file='lp_results.RDS')


######
### write results 
results_tss$Type = 'TSS'
results_ctcf$Type = 'CTCF'
results_tiles$Type ='Tiles'

write.csv(rbind(results_tss, results_ctcf, results_tiles),
          file='results_healthyDonors.csv')