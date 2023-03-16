

library(MOCHA)
library(SummarizedExperiment)

setwd('scMACS_Analysis')

STM <- readRDS('CD16_SampleTileObj.rds')

allCalls <- plyranges::read_bed('CD16 Mono_intensities_covid.bed')

## Look at immunologically important genes. NFKB-family, AP-1? TLR?  Myc?

importTSS <- rowRanges(STM) %>% plyranges::filter(grepl('NFK|IKK|TLR|MYC|P53|NOD', Gene))

missOpportunity <- plyranges::filter(allCalls, name %in% 'MOCHA_Unique') %>% join_overlap_intersect(importTSS)

plottingWindows <- missOpportunity %>% anchor_centre() %>% stretch(4000) %>% reduce_ranges()


tmp <- MOCHA::extractRegion(STM, 'CD16 Mono', groupColumn = 'COVID_status', sampleSpecific = FALSE,
                            region = plottingWindows[3], numCores = 35)

plotRegion(tmp, whichGene = 'TLR9')