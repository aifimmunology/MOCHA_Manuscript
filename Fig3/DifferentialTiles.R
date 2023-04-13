#### Code for generating Figure 4 results and plots
### Additional supporting function are placed at the very end. 

library(ArchR)
library(MOCHA)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(plyranges)
library(ggrastr)
library(ggrepel)
library(RcppAlgos)
library(WebGestaltR)

setwd('scMACS_Analysis')

## Taken from the larger ArchRProject. 
studySignal = 3628
MonoDCE <- loadArchRProject('MonoDC_Edits')

tR <- callOpenTiles( 
    MonoDCE,
    cellPopLabel = 'predictedGroup_Co2',
    cellPopulations= 'CD16 Mono',
    TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
    Org = org.Hs.eg.db,
    outDir = getOutputDirectory(MonoDCE),
    numCores = 45,
    studySignal= studySignal
)

saveRDS(tR, 'CD16_tileResults.rds')
tR <- readRDS( 'CD16_tileResults.rds')

##Filter down to just Early Infection and Uninfected samples
## This also removes repeated measures issues, 
#   because some donors have multiples samples with days_since_symptoms < 15
tR2 <- subsetMOCHAObject(tR, subsetBy = 'visit',
                         groupList = 1, na.rm = FALSE)

plotConsensus(tR, groupColumn = 'InfectionStages',
                     numCores = 25)

STM <- getSampleTileMatrix( 
    tR,
    groupColumn = 'InfectionStages',
    threshold = 0.2,
    numCores = 40
)

STM <- annotateTiles(STM)

saveRDS(STM, 'CD16_SampleTileObj.rds')


####### Run Differential Accessibility

### Let's compare COVID+ early infection (< day 15) vs Uninfected Controls

daps <- getDifferentialAccessibleTiles(STM,cellPopulation = 'CD16 Mono',
                                           groupColumn = 'InfectionStages'r,
                                           foreground = 'Early Infection',
                                           background = 'Uninfected',
                                           signalThreshold = 12,
                                           minZeroDiff = 0.5,
                                           fdrToDisplay = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 3)

write.csv(daps, 'Fig4_AllDAPs_CD16s_EarlyInfection.csv')

################ Analyze all DAPs (Supplemental Figure 4 and Supplemental Figure 6A-B)

p1 <- ggplot(as.data.frame(daps), 
             aes(x = Log2FC_C, y = -log10(P_value), 
                 color = ifelse(FDR <= 0.2, 'Significant', 'Unchanged'))) + 
    rasterise(geom_point(size = 0.025)) +  
    ylab('- Log of P-value') + xlab( 'Log2FC') + theme_bw() +
    theme(legend.position = c(0.15,0.9)) + 
    scale_color_manual(name='Significant DAP', 
                       breaks=c('Significant', 'Unchanged'),
                       values = c('Significant' = 'red',
                                  'Unchanged' = 'grey'))

pdf('Fig3_A.pdf')
p1
dev.off()