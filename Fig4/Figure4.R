
############### Figure 4 Code:

#### Code for generating Figure 4 results and plots
### Additional supporting function are placed at the very end. 

library(ArchR)
library(scMACS)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(plyranges)
library(ggrastr)
library(ggrepel)
library(RcppAlgos)
library(WebGestaltR)

setwd('scMACS_Analysis')

FullCovid <- loadArchRProject('../FullCovidf')
medFrags <- median(FullCovid$nFrags)
MonoDCE <- loadArchRProject('MonoDC_Edits')

tR <- callOpenTiles( 
    MonoDCE,
    cellPopLabel = 'predictedGroup_Co2',
    cellPopulations= 'CD16 Mono',
    TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
    Org = org.Hs.eg.db,
    numCores = 40,
    studySignal= medFrags
)

saveRDS(tR, 'CD16_tileResults.rds')
tR <- readRDS( 'CD16_tileResults.rds')

tR2 <- subsetMOCHAObject(tR, subsetBy = 'InfectionStages',
                         groupList = c('Early Infection','Uninfected'))

tR2 <- subsetMOCHAObject(tR2, subsetBy = 'visit',
                         groupList = 1, na.rm = FALSE)

plotConsensus(tR, groupColumn = 'InfectionStages',
                     numCores = 25)

STM <- getSampleTileMatrix( 
    tR2,
    groupColumn = 'InfectionStages',
    threshold = 0.2,
    numCores = 40
)

STM <- annotateTiles(STM)

saveRDS(STM, 'CD16_SampleTileObj.rds')
STM <- readRDS('CD16_SampleTileObj.rds')

################## Add peakset to ArchR Project

MonoDCE <- addPeakSet(MonoDCE,rowRanges(STM), force = TRUE)
MonoDCE <- addMotifAnnotations(MonoDCE,  motifSet = "cisbp", name = "cisbpMotif", force = TRUE)
saveArchRProject(MonoDCE)


####### Run Differential Accessibility


##Filter down to just Early Infection and Uninfected samples
## This also removes repeated measures issues, 
#   because some donors have multiples samples with days_since_symptoms < 15
SampleData <- colData(STM)[colData(STM)$InfectionStages %in% c('Uninfected','Early Infection'),] %>%
                    as.data.frame()
SampleData2 <- SampleData %>% dplyr::group_by(Subject) %>% arrange(days_since_symptoms) %>%
                    slice(1)

### Let's compare COVID+ early infection (< day 15) vs Uninfected Controls

daps <- getDifferentialAccessibleTiles(STM,cellPopulation = 'CD16 Mono',
                                           groupColumn = 'InfectionStages',
                                           foreground = 'Early Infection',
                                           background = 'Uninfected',
                                           signalThreshold = 12,
                                           minZeroDiff = 0.5,
                                           fdrToDisplay = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 25)

write.csv(daps, 'Fig4_AllDAPs_CD16s_EarlyInfection.csv')
daps <- makeGRangesFromDataFrame(read.csv('Fig4_AllDAPs_CD16s_EarlyInfection.csv'), 
                                 keep.extra.columns = TRUE)

allTSS <- getAltTSS(sort(daps), returnAllTSS = 'TRUE')    
write.csv(allTSS, 'Fig4_AllTSS_CD16s_EarlyInfection.csv')

allTSS <- makeGRangesFromDataFrame(read.csv('AllTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)


################ Analyze all DAPs (Supplemental Figure 4 and Supplemental Figure 6A-B)

p1 <- ggplot(as.data.frame(daps), 
             aes(x = Log2FC_C, y = -log10(P_value), 
                 color = ifelse(FDR < 0.2, 'Significant', 'Unchanged'))) + 
    rasterise(geom_point(size = 0.025)) +  
    ylab('- Log of P-value') + xlab( 'Log2FC') + theme_bw() +
    theme(legend.position = c(0.15,0.9)) + 
    scale_color_manual(name='Significant DAP', 
                       breaks=c('Significant', 'Unchanged'),
                       values = c('Significant' = 'red',
                                  'Unchanged' = 'grey'))
pdf('Fig4_A.pdf')
p1
dev.off()


##DAP pathway analysis: Reactome pathway

enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

allTSS_geneSet <- plyranges::filter(allTSS, FDR < 0.2) %>% 
                    group_by(name) %>% summarize(plyranges::n()) 
allTSSReactome<- simplifiedORA(enrichDataBaseList$name[9],
                               as.character(allTSS_geneSet$name), refList)
allTSSReactome$description = 
                        factor(allTSSReactome$description, levels = 
                               c(grep('Toll|TLR',allTSSReactome$description,value = TRUE),
                                 grep("-",grep('D88',allTSSReactome$description,value = TRUE),value=TRUE,invert = TRUE),
                                 grep('Interleukin',allTSSReactome$description,value = TRUE),
                                 grep('Check|53|AKT|kine|MAPK|MAP kinase ac',allTSSReactome$description,value = TRUE),
                                 grep('BCR|B cells|RHO',allTSSReactome$description,value = TRUE),
                                 grep('degen',allTSSReactome$description,value = TRUE)))
                        
pdf('SuppFig6_A.pdf')

ggplot(allTSSReactome, aes(x = str_wrap(description,40),
                           y = enrichmentRatio, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('Term Enrichment') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Reactome Pathway Enrichment in All Genes with Differential TSS')



dev.off()

#### Run MotifEnrichment

library(chromVAR)
library(chromVARmotifs)
data(human_pwms_v2)
STM <- addMotifSet(STM, human_pwms_v2, motifSetName = 'CISBP')
posList <- STM@metadata$CISBP

enrichDAP_df <- MotifEnrichment(filter(daps, FDR < 0.2),
                                filter(daps, FDR >= 0.2), posList, numCores = 40)
    
write.csv(enrichDAP_df, 'Fig4_CD16s_DAPs_MotifEnrichment.csv')


enrichDAP_df2 <- enrichDAP_df  %>% 
                dplyr::mutate(TranscriptionFactor = gsub("_.*", "", rownames(.)),
                                                mlog10Padj = -log10(adjp_val)) %>%
                 dplyr::mutate(label = ifelse(adjp_val < 0.05, TranscriptionFactor, NA)) 

pdf('SuppFig6_B.pdf')

    ggplot(enrichDAP_df2, aes(x = enrichment, y = -log10(adjp_val), label = label)) + 
        geom_point() + theme_bw() + 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) + 
        ggtitle('Motif Enrichment at DAPs for CD16 in Early Infection') +
        ylab('-Log of Adjusted P-Value')+ 
        xlab('Enrichment')

dev.off()


################################# Introduce alternative TSS analysis

##We need to make a list of all combinations of TSS sites per gene. 
##And all we need to pull out is the Log2FC_C for site 1 vs site 2
                          
sumTSS <- allTSS %>% group_by(tx_id, name) %>% 
                        reduce_ranges(Log2FC_C = max(abs(Log2FC_C))*sign(sum(MeanDiff)),
                                                             TestStatistic = max(Test_Statistic),
                                                             Signal = -log10(min(FDR))*sign(sum(MeanDiff)),
                                                             MeanDiff = max(abs(MeanDiff))*sign(sum(MeanDiff)),
                                                             P_value = min(P_value),
                                           FDR = min(FDR)) %>%
                          group_by(name) %>% filter(any(duplicated(name))) %>%
                          filter(!is.na(name)) %>% ungroup() %>% sort()

pairwise_combos = RcppAlgos::comboGrid(c(1:length(sumTSS$tx_id)),
                                       c(1:length(sumTSS$tx_id)), repetition = FALSE)
#Filter for pairs within the same gene
pairwise_combos = pairwise_combos[sumTSS$name[pairwise_combos[,'Var1']] == 
                                  sumTSS$name[pairwise_combos[,'Var1']],]
pairwise_combos = pairwise_combos[which(as.character(sumTSS$name)[pairwise_combos[,'Var1']] == 
                                  as.character(sumTSS$name)[pairwise_combos[,'Var2']]),]
colsForComp <- c('Log2FC_C', 'TestStatistic', 'Signal', 'MeanDiff', 'P_value', 'FDR')

GR_meta <- as.data.frame(mcols(sumTSS))

allComb <- mclapply(colsForComp, function(x){
    
                tmp <- data.frame(GR_meta[pairwise_combos[,'Var1'], x],
                                  GR_meta[pairwise_combos[,'Var2'], x])
                colnames(tmp) <- paste(c('Site1', 'Site2'), x, sep = "_")
                tmp
    
    }, mc.cores = 15) %>%
    do.call('cbind',.) %>% 
    dplyr::mutate(Significant = ifelse((Site1_FDR < 0.2 | Site2_FDR < 0.2),
                                                   'Differential', 'Unchanged')) %>%
    dplyr::mutate(Significant = ifelse(is.na(Significant), 'Unchanged', Significant)) %>% 
    cbind(data.frame(Gene = GR_meta[pairwise_combos[,'Var1'], 'name']), .)

pdf('Fig4_B.pdf')

ggplot(allComb, aes(x = Site1_Signal, y = Site2_Signal, color = Significant)) + geom_point(alpha = 0.5) + 
    xlab('-log of FDR * Direction of Change at Site 1') + 
    geom_hline(yintercept = c(-log10(0.2), log10(0.2)))+ 
    geom_vline(xintercept = c(-log10(0.2), log10(0.2)))  +
    theme_bw()    +
    theme(legend.position = c(0.9, 0.1)) + 
    ylab('-log of FDR * Direction of Change at Site 2')
       
dev.off()

########################## Analyze Alternative TSS Usage

####### Identify Alternative TSS usage

altTSS <- getAltTSS(sort(daps))    

write.csv(altTSS, 'AlternativeTSS_CD16s_EarlyInfection.csv')
altTSS <- makeGRangesFromDataFrame(read.csv('AlternativeTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)


####### Pathway enrichment at altTSS genes

##DAP pathway analysis: Reactome pathway

enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

altTSS_geneSet <- unique(altTSS$name) 
altTSSReactome<- simplifiedORA(enrichDataBaseList$name[9], as.character(altTSS_geneSet), refList)

altTSSReactome$description = 
    factor(altTSSReactome$description, levels =  c(
    grep('degenerat',altTSSReactome$description,value = TRUE),
    grep('beta|WNT|Cell|GTP|AKT1',altTSSReactome$description,value = TRUE),
    
    grep('Toll|TLR|MAP kinase|MAP3K',altTSSReactome$description,value = TRUE),
grep("-",grep('D88',altTSSReactome$description,value = TRUE),value=TRUE,invert = TRUE),
grep('Interleukin',altTSSReactome$description,value = TRUE)))
          
pdf('Fig4_C.pdf')


ggplot(altTSSReactome, aes(x = str_wrap(description,40),
                           y = enrichmentRatio, fill = FDR)) +
            geom_col() +  coord_flip() + theme_bw() + 
                        xlab('Term Description') + 
                        ylab('Term Enrichment') + 
            theme(legend.position = c(0.8, 0.8))+ 
            ggtitle('Reactome Pathway Enrichment in alt Genes with Differential TSS')


dev.off()



##################### Let's plot all the alternative TSS sites

geneList <- c('CAPN1','ARHGAP9','NFKB1', 'SOCS3', 'P2RY6', 'IFI16')

STM <- readRDS('CD16_SampleTileObj.rds')

pdf('Fig4D.pdf')

for(i in 1:length(geneList)){
    
   tmpG <- altTSS %>% filter(name %in% geneList[i]) %>% 
        arrange(start) %>%
            stretch(., extend = 2000)
    reg = paste(unique(seqnames(tmpG)),":",min(start(tmpG)),"-",
                max(end(tmpG)),sep="")
    tmp1 <- extractRegion(STM, reg, groupColumn = 'InfectionStages',
                          subGroups = c('Early Infection','Uninfected'),
                          numCores = 5, sampleSpecific = FALSE)
    
    tmp1 <- addAccessibilityShift(tmp1, 'CD16 Mono.Early Infection','CD16 Mono.Uninfected',
                     'Accessibility Change')
    print(plotRegion(tmp1, whichGene =geneList[i], 
                              legend.position  = "none"))
    
}

dev.off()

    
################################## Let's run peak co-accessibility for one example region
altTSS <- makeGRangesFromDataFrame(read.csv('AlternativeTSS_CD16s_EarlyInfection.csv'),
                                   keep.extra.columns= TRUE)
STM <- readRDS('CD16_SampleTileObj.rds')

## Get correlations with the SOCS3 TSS sites.
SitesOfInterest <- altTSS %>% filter(name == 'SOCS3' & FDR <= 0.2) 
    

loops1 <- getCoAccessibleLinks(STM, cellPopulation= 'CD16 Mono',
                               regions = SitesOfInterest, 
                                windowSize = 0.5*10^6, numCores = 15)

#Filter loops for those with a min absolute Correlation of 0.5
loopsf <- filterCoAccessibleLinks(loops1, threshold = 0.5)

## Get region that includes the correlation windows. 
chrReg = GRanges(  seqnames = unique(loopsf$chr), 
                     ranges = as_iranges(data.frame(start=min(loopsf$start), end = max(loopsf$end))))

## Extract ands count for that region
SOCS3Count <- extractRegion(STM, groupColumn = 'InfectionStages', chrReg, numCores = 10)
saveRDS(SOCS3Count, 'ExtractedCounts_SOCS3.rds')

pdf('Fig4E.pdf')

plotRegion(SOCS3Count, collapseGenes = 'longestTx', 
           linkdf = loopsf, #plotType = 'line',  
          relativeHeights = c(`Chr` = 0.9, `Normalized Counts` = 7, 
                              `Links` =1.5, `Genes` = 2))
dev.off()

orgdb <-  AnnotationDbi::loadDb(SOCS3Count@metadata$Org)

tmp <- scMACS:::simplifiedOrgDb(TxDb = TxDb, orgdb = orgdb)

Homo.sapiens.hg38 <- simplifiedOrgDb(TxDb = AnnotationDbi::loadDb(SOCS3Count@metadata$TxDb), orgdb =  AnnotationDbi::loadDb(SOCS3Count@metadata$Org))

############### Let's run co-accessibility for just the Alternative TSS sites

altTSSLinks <- getCoAccessibleLinks(STM, cellPopulation = 'CD16 Mono',
                                     regions = plyranges::filter(altTSS, FDR < 0.2), 
                                     numCores = 30)

saveRDS(altTSSLinks , 'Links_AltTSS_Final.rds')
altTSSLinks <- readRDS('Links_AltTSS_Final.rds')

altTSSLinksf <- filterCoAccessibleLinks(altTSSLinks, threshold = 0.5)

altTSS_Network <- c(altTSSLinksf$Peak1, altTSSLinksf$Peak2, 
                           GRangesToString(plyranges::filter(altTSS, FDR < 0.2))) %>%
                   unique() %>%
                  StringsToGRanges(.) %>% 
                plyranges::filter_by_overlaps(daps, .)

############################################# Find motif Enrichment around alt TSS sites

## AltTSS DAPs, vs TSS DAPs
foreGround = altTSS_Network

## Find background network
specTSS <- filter_by_non_overlaps(allTSS, altTSS_Network)
startT <- Sys.time()
tssLinks <- FindCoAccessibleLinks(results[[2]], specTSS, numCores = 50)
endT <- Sys.time() - startT    

saveRDS(tssLinks, 'Links_BackgroundTSS.rds')
tssLinks <- readRDS('Links_BackgroundTSS.rds')

tssLinksf <- FilterCoAccessibleLinks(tssLinks, threshold = 0.5)
allTSS_Network <- c(tssLinksf$Peak1, tssLinksf$Peak2, GRangesToString(specTSS)) %>%
                   unique() %>% StringsToGRanges(.) %>% 
                plyranges::filter_by_overlaps(daps, .)
backGround = filter_by_non_overlaps(allTSS_Network, altTSS_Network)

## Get motifs
posList <- getPositions(MonoDCE, 'cisbpMotif')

## Run enrichment: AltTSS Network vs all TSS Network, DAPs vs nonDAPs
enrich_df <- MotifEnrichment(foreGround,backGround, posList, numCores = 40)
    
write.csv(enrich_df, 'CD16_MotifEnrichment_AltTSS.csv')

enrich_df <- read.csv('CD16_MotifEnrichment_AltTSS.csv', row.names = TRUE)


enrich_df2 <- enrich_df  %>% 
                dplyr::mutate(TranscriptionFactor = gsub("_.*", "", rownames(.)),
                                                mlog10Padj = -log10(adjp_val)) %>%
                dplyr::mutate(label = ifelse(adjp_val < 0.05, TranscriptionFactor, NA))


pdf('Fig4F.pdf')

    ggplot(enrich_df2, aes(x = enrichment, y = -log10(adjp_val), label = label)) + 
        geom_point() + theme_bw() + 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) + 
        ggtitle('Motif Enrichment at Alternative TSS for CD16 in Early Infection') +
        ylab('-Log of Adjusted P-Value')+ 
        xlab('Enrichment')

    ggplot(enrichDAP_df2, aes(x = enrichment, y = -log10(adjp_val), label = label)) + 
        geom_point() + theme_bw() + 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) + 
        ggtitle('Motif Enrichment at DAPs for CD16 in Early Infection') +
        ylab('-Log of Adjusted P-Value')+ 
        xlab('Enrichment')

dev.off()


################################## Run Ligand-MSEA on the results


ligandtf <- readRDS('../Resubmission COVID19/ligand_tf_matrix.rds')   

filteredTF <- ligandtf[rownames(ligandtf) %in% unique(enrich_df2$TranscriptionFactor), ]             
                    

ligandMSEA <- MotifSetEnrichmentAnalysis(ligandtf, enrich_df2, 
                                       columnWithMotifs = "TranscriptionFactor", 
                                       ligands = colnames(filteredTF)[colSums(filteredTF) > 0],
                                       mlogPval_threshold = 2, 
                                       annotationName = 'CellType', annotation = "none", 
                                       numCores = 30, verbose = FALSE) %>%
                        mutate(Label = ifelse(adjp_val < 0.05, ligand, NA))

write.csv(ligandMSEA, 'CD16_AltTSS_LigandMSEA.csv')

pdf('Fig4G.pdf')
ggplot(ligandMSEA, aes(x = PercInNicheNet, y = -log10(adjp_val), label  = Label)) +
                    geom_point() + #geom_hline(yintercept = -log10(0.05), aes(color = 'red')) + 
                    theme_bw()+
                    geom_text_repel() + 
                    xlab('Percent of Related TFs in Ligand Pathway') + ylab('-Log of Adj P-value') +
                    ggtitle('Potential Ligands Influencing Alternative TSS Usage')
dev.off()


######### Pull out Motifs associated with each gene. 

Gene2MotifLinks <- Gene2Motif(TSS_Sites = filter(altTSS, FDR < 0.2), 
                              allPeaks = daps, 
                              TSS_Links = altTSSLinksf, 
                              motifPosList = posList, numCores = 20)

saveRDS(Gene2MotifLinks, 'CD16_AltTSS_Network_Gene2Motif.rds')

altTSSLinks <- readRDS('Links_AltTSS.rds')
altTSSLinksf <- FilterCoAccessibleLinks(altTSSLinks, threshold = 0.5)

tssLinks <- readRDS('Links_BackgroundTSS.rds')
tssLinksf <- FilterCoAccessibleLinks(tssLinks, threshold = 0.5)

allLinks <- rbind(tssLinksf, altTSSLinksf)
Gene2MotifLinks <- Gene2Motif(TSS_Sites = filter(allTSS, FDR < 0.2), 
                              allPeaks = daps, 
                              TSS_Links = allLinksf, 
                              motifPosList = posList, numCores = 20)
saveRDS(Gene2MotifLinks2, 'AllDAP_TSS_Links.RDS')

## Generate total Motif - Gene linkage
sigMotifs <- filter(enrichDAP_df,adjp_val < 0.05)  %>% mutate(Type = 'Motif') %>%
                dplyr::rename(AdjP = adjp_val, EffectSize = enrichment, name = TranscriptionFactor) %>%
                dplyr::select(name, AdjP, EffectSize, Type)
# Pull in EffectSize for strongest DAPs per gene
geneEffect <- filter(allTSS, FDR < 0.2) %>% group_by(name) %>% 
    summarize(EffectSize = max(abs(Log2FC_C)),
              AdjP = min(FDR)) %>%
    as.data.frame()  %>%
    mutate(Type = 'Gene') %>%
    dplyr::select(name, AdjP, EffectSize, Type)

GM_df <- melt(Gene2MotifLinks2) %>% dplyr::rename(Motif = value, Gene = L1) %>%
            dplyr::mutate(Motif = gsub("_.*","", Motif)) %>%
            dplyr::filter(Motif %in% sigMotifs$name) 

sigLigands <- dplyr::filter(ligandMSEA2, adjp_val <  0.05)\
filtedLigands <- filter(ligandMSEA2, adjp_val < 0.05) %>% mutate(Type = 'Ligand') %>%
                       rename(name = ligand, AdjP = adjp_val, EffectSize = PercInNicheNet)  %>%
                        dplyr::select(name, AdjP, EffectSize, Type)
LM_df <- ligandtf[rownames(ligandtf) %in% 
                        unique(sigMotifs$TranscriptionFactor), sigLigands$ligand] %>%
                melt() %>% dplyr::filter(value > 0) %>% dplyr::rename(Motif = Var1, Ligand = Var2) %>%
                dplyr::select(Ligand,Motif)
full_df <- full_join(LM_df, GM_df, by = 'Motif')
part1 <- dplyr::select(full_df, Ligand, Motif) %>% 
            rename(From = Ligand, To = Motif)
part2 <- dplyr::select(full_df, Motif, Gene) %>% 
            rename(From = Motif, To = Gene)
both1 <- distinct(rbind(part1, part2))

write.table(both1, 'DAPS_Network_Edges.tsv', sep='\t')

sigMotifs$WeightPVal <- rescale(sigMotifs$AdjP)
sigMotifs$WeightEffect <- rescale(sigMotifs$EffectSize)
geneEffect$WeightPVal <- rescale(geneEffect$AdjP)
geneEffect$WeightEffect <- rescale(geneEffect$EffectSize)
filtedLigands$WeightPVal <- rescale(filtedLigands$AdjP)
filtedLigands$WeightEffect <- rescale(filtedLigands$EffectSize)

NodeTable <- do.call('rbind', list(sigMotifs, geneEffect, filtedLigands))
write.csv(NodeTable, 'DAPs_NodeTable.csv')



############################################################################################

###### Functions used. 

############################################################################################
                        
simplifiedORA <- function(database, foreground, background){
    
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                            referenceGeneType= "genesymbol")
    
    }