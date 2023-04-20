
############### Figure 4 Code:

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


STM <- readRDS('CD16_SampleTileObj.rds')

################## Add peakset to ArchR Project

MonoDCE <- addPeakSet(MonoDCE,rowRanges(STM), force = TRUE)
MonoDCE <- addMotifAnnotations(MonoDCE,  motifSet = "cisbp", name = "cisbpMotif", force = TRUE)
saveArchRProject(MonoDCE)


################## Find TSSs that are differential
daps <- makeGRangesFromDataFrame(read.csv('Fig4_AllDAPs_CD16s_EarlyInfection.csv'), 
                                 keep.extra.columns = TRUE)

allTSS <- getAltTSS(sort(daps), returnAllTSS = 'TRUE')    
write.csv(allTSS, 'Fig4_AllTSS_CD16s_EarlyInfection.csv')

allTSS <- makeGRangesFromDataFrame(read.csv('Fig4_AllTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)


### Pathway enrichment 
enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")
                           

allTSS <- makeGRangesFromDataFrame(read.csv('Fig4_AllTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)


                           
allTSSReactome <- simplifiedORA(enrichDataBaseList$name[9],
                               as.character(allTSS_geneSet), refList)
                        
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
PathAnnot = lapply(allTSSReactome$description, function(y)
                    unlist(lapply(y, 
                             function(x) findTrunk(x, Hierarch, ReactID))))  %>% unlist()
allTSSReactome$TopLevel = PathAnnot
allTSSReactome$SecondLevel = dplyr::case_when(
                        grepl('TLR|MAP Kinase|MyD88|MAP|Toll', allTSSReactome$description) ~ 'Innate Immune System',
                        grepl('euro|Cancer', allTSSReactome$description) ~ 'Disease',
                        grepl('Interleukin', allTSSReactome$description) ~ 'Signaling by Interleukins',
                        grepl('B cell|B Cell', allTSSReactome$description) ~ 'Adaptive Immune System',
                        grepl('Cell Cycle|kinetochores|Mitotic', allTSSReactome$description) ~ 'Cell Cycle',
                        grepl('UPR|GTPas|WNT|TP53|MAPK1/3', allTSSReactome$description) ~ 'Other')
                                  
allTSSReactome <- arrange(allTSSReactome, SecondLevel, enrichmentRatio)
allTSSReactome$description = factor(allTSSReactome$description, levels = unique(allTSSReactome$description))

write.csv(allTSSReactome, 'SuppDataFile2_allTSS_Differential_ReactomePathways.csv')                           
                        
                        
pdf('SuppFig6_A.pdf')

ggplot(allTSSReactome, aes(x = description,
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
            
saveRDS(STM, 'CD16_SampleTileObj.rds')
                           
posList <- STM@metadata$CISBP

enrichDAP_df <- MotifEnrichment(filter(daps, FDR < 0.2),
                                filter(daps, FDR >= 0.2), posList, numCores = 5)
    
write.csv(enrichDAP_df, 'Fig4_CD16s_DAPs_MotifEnrichment.csv')

tmp <- read.csv('Fig4_CD16s_DAPs_MotifEnrichment.csv')
colnames(tmp)[1] = 'TranscriptionFactor'      
                           
enrichDAP_df2 <- enrichDAP_df  %>% 
                dplyr::mutate(TranscriptionFactor = gsub("_.*", "", rownames(.))) %>%
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

pdf('Fig4_A.pdf')

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
                           
length(unique(altTSS$name))
                           
##Annotate by type:
                           
altTSS <- altTSS %>% group_by(name) %>% 
    mutate(Type = dplyr::case_when(
        any(Log2FC_C > 0 & FDR <= 0.2) & any(Log2FC_C < 0 & FDR <= 0.2) ~ 'Type II',
                                   TRUE ~ 'Type I'))

write.csv(altTSS, 'AlternativeTSS_CD16s_EarlyInfection.csv')
altTSS <- makeGRangesFromDataFrame(read.csv('AlternativeTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)
                                  
## Annotate with previous analyzed differential gene expression from COVID19 infection

alterDEGs <- read.csv('covid_alteredDEGs.csv') %>% 
            filter(P.value.6h < 0.05 | P.value.2h < 0.05 | P.value.10h < 0.05 |
                   P.value.24h < 0.05)
dim(alterDEGs)
sum(unique(altTSS$name) %in% alterDEGs$Gene.Symbol01)
                                  
alterProtein <- read.csv('covid_alteredProtein.csv') %>% 
            filter(P.value.6h < 0.05 | P.value < 0.05 | P.value.10h < 0.05 |
                   P.value.24h < 0.05)
dim(alterProtein)
sum(unique(altTSS$name) %in% alterProtein$Gene.Symbol)
 
sum(unique(altTSS$name) %in% c(alterProtein$Gene.Symbol,alterDEGs$Gene.Symbol01))
#56 out of 283             

####### Pathway enrichment at altTSS genes

##DAP pathway analysis: Reactome pathway

enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

altTSS_geneSet <- unique(altTSS$name) 
altTSSReactome<- simplifiedORA(enrichDataBaseList$name[9], as.character(altTSS_geneSet), refList)

                           
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                           
PathAnnot = lapply(altTSSReactome$description, function(y)
                    unlist(lapply(y, 
                             function(x) findTrunk(x, Hierarch, ReactID))))  %>% unlist()
altTSSReactome$TopLevel = PathAnnot   
altTSSReactome <- altTSSReactome %>% 
                    mutate(SecondLevel = dplyr::case_when(
                        grepl('TLR|MAP|MyD88|MAP|Toll', description) ~ 'TLR Cascades',
                        grepl('euro|Cancer', description) ~ 'Disease',
                        grepl('Interleukin', description) ~ 'Interleukin',
                        grepl('UPR|GTPas|WNT', description) ~ 'Other'))
                        
altTSSReactome <- arrange(altTSSReactome, SecondLevel, enrichmentRatio)
altTSSReactome$description = factor(altTSSReactome$description, levels = unique(altTSSReactome$description))
                                  
write.csv(altTSSReactome, 'AltTSS_Pathways.csv')
          
pdf('Fig4_B.pdf')


ggplot(altTSSReactome, aes(x = description,40,
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

pdf('ExamplePlot_SOCS3_Links.pdf')

plotRegion(SOCS3Count, collapseGenes = 'longestTx', 
           linkdf = loopsf, #plotType = 'line',  
          relativeHeights = c(`Chr` = 0.9, `Normalized Counts` = 7, 
                              `Links` =1.5, `Genes` = 2))
dev.off()


############### Let's run co-accessibility for just the Alternative TSS sites
altTSS <- makeGRangesFromDataFrame(read.csv('AlternativeTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)
daps <- makeGRangesFromDataFrame(read.csv('Fig4_AllDAPs_CD16s_EarlyInfection.csv'), 
                                 keep.extra.columns = TRUE)
allTSS <- makeGRangesFromDataFrame(read.csv('Fig4_AllTSS_CD16s_EarlyInfection.csv'), keep.extra.columns = TRUE)
                                  
altTSSLinks <- getCoAccessibleLinks(STM, cellPopulation = 'CD16 Mono',
                                     regions = plyranges::filter(altTSS, FDR < 0.2), 
                                     chrChunks = 8,
                                    approximateTile = TRUE,
                                     numCores = 55)
STM <- readRDS('CD16_SampleTileObj.rds')                                  

saveRDS(altTSSLinks , 'Links_AltTSS.rds') 
altTSSLinks <- readRDS('Links_AltTSS.rds')
                                  
## older results saved as '_Final'

altTSSLinksf <- filterCoAccessibleLinks(altTSSLinks, threshold = 0.5)

altTSS_Network <- c(altTSSLinksf$Tile1, altTSSLinksf$Tile2, 
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
tssLinks <- getCoAccessibleLinks(STM, cellPopulation = 'CD16 Mono',
                                     regions = specTSS, 
                                     chrChunks = 3,
                                    approximateTile = TRUE,
                                     numCores = 30)
endT <- Sys.time() - startT    

saveRDS(tssLinks, 'Links_BackgroundTSS_Unfiltered.rds')
tssLinks <- readRDS('Links_BackgroundTSS_Unfiltered.rds')

tssLinksf <- filterCoAccessibleLinks(tssLinks, threshold = 0.5)
allTSS_Network <- c(tssLinksf$Tile1, tssLinksf$Tile2) %>%
                   unique() %>% StringsToGRanges(.) %>% c(.,specTSS) %>%
                plyranges::filter_by_overlaps(daps, .)

backGround = filter_by_non_overlaps(allTSS_Network, altTSS_Network)

## Get motifs
posList <- metadata(STM)$CISBP

## Run enrichment: AltTSS Network vs all TSS Network, DAPs vs nonDAPs

enrich_df <- MotifEnrichment(altTSS_Network,backGround, posList)
                    
write.csv(enrich_df, 'CD16_MotifEnrichment_AltTSS_v2.csv')

enrich_df <- read.csv('CD16_MotifEnrichment_AltTSS_v2.csv', row.names= 1)

enrich_df2 <- enrich_df  %>% 
                dplyr::mutate(TranscriptionFactor = gsub("_.*", "", rownames(.))) %>%
                dplyr::mutate(label = ifelse(adjp_val < 0.05, TranscriptionFactor, NA))

enrichDAP_df <- MotifEnrichment(filter(daps, FDR <= 0.2),filter(daps, FDR >= 0.2),
                                posList) %>% 
                dplyr::mutate(TranscriptionFactor = gsub("_.*", "", rownames(.))) %>%
                dplyr::mutate(label = ifelse(adjp_val < 0.05, TranscriptionFactor, NA))
write.csv(enrichDAP_df, 'CD16_MotifEnrichment_All_DAPs.csv')
                                  
pdf('Fig5F.pdf')

    ggplot(enrich_df2, aes(x = enrichment, y = mlog10Padj, label = label)) + 
        geom_point() + theme_bw() + 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) + 
        ggtitle('Motif Enrichment at Alternative TSS for CD16 in Early Infection') +
        ylab('-Log of Adjusted P-Value')+ 
        xlab('Enrichment')

    ggplot(enrichDAP_df, aes(x = enrichment, y = mlog10Padj, label = label)) + 
        geom_point() + theme_bw() + 
    geom_text_repel(max.overlaps = 80,show.legend = FALSE) + 
        ggtitle('Motif Enrichment at DAPs for CD16 in Early Infection') +
        ylab('-Log of Adjusted P-Value')+ 
        xlab('Enrichment')

dev.off()


################################## Run Ligand-MSEA on the results


ligandtf <- readRDS('ligand_tf_matrix.rds')   

filteredTF <- ligandtf[rownames(ligandtf) %in% unique(enrich_df2$TranscriptionFactor), ]                  

ligandMSEA <- MotifSetEnrichmentAnalysis(ligandtf, enrich_df2, 
                                       motifColumn = "TranscriptionFactor", 
                                       ligands = colnames(filteredTF)[colSums(filteredTF) > 0],
                                        stat_column = 'mlog10Padj',
                                       stat_threshold = 2, 
                                       annotationName = 'CellType', annotation = "none", 
                                       numCores = 30, verbose = FALSE) %>%
                        mutate(Label = ifelse(adjp_val < 0.05, ligand, NA)) 
                                  
ligandMSEA[ligandMSEA$adjp_val < 0.05,]
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
                              allTiles = daps, 
                              TSS_Links = altTSSLinksf, 
                              motifPosList = posList, numCores =28)

saveRDS(Gene2MotifLinks, 'CD16_AltTSS_Network_Gene2Motif_final.rds')
##final is the latest version
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
sigMotifs <- filter(enrich_df2,adjp_val < 0.05)  %>% mutate(Type = 'Motif') %>%
                dplyr::rename(AdjP = adjp_val, EffectSize = enrichment, name = TranscriptionFactor) %>%
                dplyr::select(name, AdjP, EffectSize, Type)
# Pull in EffectSize for strongest DATs per gene
geneEffect <- filter(altTSS, FDR < 0.2) %>% group_by(name) %>% 
    summarize(EffectSize = max(abs(Log2FC_C)),
              AdjP = min(FDR)) %>%
    as.data.frame()  %>%
    mutate(Type = 'Gene') %>%
    dplyr::select(name, AdjP, EffectSize, Type)

GM_df <- melt(Gene2MotifLinks) %>% dplyr::rename(Motif = value, Gene = L1) %>%
            dplyr::mutate(Motif = gsub("_.*","", Motif)) %>%
            dplyr::filter(Motif %in% sigMotifs$name) 

sigLigands <- dplyr::filter(ligandMSEA, adjp_val <  0.05)
filtedLigands <- filter(ligandMSEA, adjp_val < 0.05) %>% mutate(Type = 'Ligand') %>%
                       rename(name = ligand, AdjP = adjp_val, EffectSize = PercInNicheNet)  %>%
                        dplyr::select(name, AdjP, EffectSize, Type)
LM_df <- ligandtf[rownames(ligandtf) %in% 
                        unique(sigMotifs$name), sigLigands$ligand] %>%
                melt() %>% dplyr::filter(value > 0) %>% dplyr::rename(Motif = Var1, Ligand = Var2) %>%
                dplyr::select(Ligand,Motif)
colnames(LM_df) <- c('From', 'To')
colnames(GM_df) <- c('From', 'To')                                  

write.table(rbind(LM_df, GM_df), 'DAPS_Network_Edges_final.tsv', sep='\t')
                                  
NodeTable <- do.call('rbind', list(sigMotifs, geneEffect, filtedLigands))
                                  
notes <- read.csv('SupplementalTable_X_CEBPA_NodeTable.csv') %>% 
                                  dplyr::select(name, COVID, Detail, Reference)
NodeTable <- left_join(NodeTable, notes, by = 'name')                                  
write.csv(NodeTable, 'DAPs_NodeTable_final.csv')
                                  
                        
