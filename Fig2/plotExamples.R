

library(MOCHA)
library(SummarizedExperiment)
library(plyranges)

setwd('scMACS_Analysis')

STM <- readRDS('CD16_SampleTileObj.rds')

allCalls <- plyranges::read_bed('CD16 Mono_intensities_covid.bed')

## Look at immunologically important genes. NFKB-family, AP-1? TLR?  Myc?

importTSS <- rowRanges(STM) %>% plyranges::filter(grepl('NFK|IKK|TLR|MYC|P53|NOD', Gene)) %>%
                plyranges::filter(tileType != 'Intragenic')

missOpportunity <- plyranges::filter(allCalls, name %in% 'MOCHA_Unique') %>% join_overlap_intersect(importTSS)

plottingWindows <- missOpportunity %>% anchor_centre() %>% stretch(10000) %>% reduce_ranges()

regionGRanges <- plottingWindows[1]


regionList <- lapply(seq_along(plottingWindows), function(x){
                       extractRegion(STM, 'CD16 Mono', 
                            groupColumn = 'COVID_status', sampleSpecific = FALSE,
                            region = plottingWindows[x], numCores = 35)
            })
saveRDS(regionList, 'DifferentRegions.RDS')

specGeneLists <- lapply(seq_along(plottingWindows), function(x){
                            filter_by_overlaps(missOpportunity, plottingWindows[x])
    })


addGRangesTrack <- lapply(seq_along(plottingWindows), function(x){
                            filter_by_overlaps(allCalls, plottingWindows[x])
    })

        
names(regionList) = c('TLR9', 'TLR3', 'MYCBP2', 'NOD2')

pdf('MissingOpportunity_Examples.pdf')

for( i in seq_along(regionList)){
    print(i)
     p1 <- plotRegion(regionList[[i]], whichGene= names(regionList)[i],
                      additionalGRangesTrack = addGRangesTrack[[i]],
                relativeHeights = c(`Chr` = 0.9, `Normalized Counts` = 7, 
                                    `Genes` = 1.5, `AdditionalGRanges` = 3))
                      
    print(p1)
    
    }
dev.off()

## TLR3 in particular is interesting: We need to plot the entire gene?

regionGRanges <- StringsToGRanges('chr4:186067000-186089000')
TLRRegion <-  extractRegion(STM, 'CD16 Mono', 
                            groupColumn = 'COVID_status', sampleSpecific = FALSE,
                            region = regionGRanges , numCores = 35)
addTracks <-   filter_by_overlaps(allCalls, regionGRanges)

pdf('Larger_TLR3_Region.pdf')
plotRegion(TLRRegion, whichGene= 'TLR3',
                      additionalGRangesTrack = addTracks, legend.position = 'none',
                relativeHeights = c(`Chr` = 0.9, `Normalized Counts` = 7, 
                                    `Genes` = 1.5, `AdditionalGRanges` = 3))
dev.off()


