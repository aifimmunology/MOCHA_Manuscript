

library(MOCHA)
library(GenomicRanges)

#Load ArchR DATs
ArchRDAT <- makeGRangesFromDataFrame(read.csv('Fig3/generateData/cd16_ArchR.csv')[,-1], keep.extra.columns = TRUE)

#Load Signac DATs
SignacTiles <- read.csv('Fig3/generateData/cd16_signac.csv')
SignacDAT <- MOCHA::StringsToGRanges(sub("-",":",SignacTiles[,1]))
mcols(SignacDAT) <- SignacTiles


#Now run getAltTSS - which requires specific columns
SignacDAT$FDR = SignacDAT$p_val_adj
SignacDAT$Log2FC_C = SignacDAT$avg_log2FC

ArchRDAT$Log2FC_C = ArchRDAT$Log2FC

### getAltTSS needs an FDR

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

archrAlt <- getAltTSS(ArchRDAT, 
                      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene, 
                      OrgDb = org.Hs.eg.db)
                      
signacAlt <- getAltTSS(SignacDAT, 
                      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene, 
                      OrgDb = org.Hs.eg.db)

## Now let's load MOCHA for comparison. 
altTSS <- makeGRangesFromDataFrame(read.csv('../AlternativeTSS_CD16s_EarlyInfection.csv'), 
                                   keep.extra.columns = TRUE)


write.csv(archrAlt, 'AltTSS_ArchR.csv')
write.csv(signacAlt, 'AltTSS_Signac.csv')

#######################  Now run pathway enrichment

####### Pathway enrichment at altTSS genes

##DAP pathway analysis: Reactome pathway

enrichDataBaseList <- WebGestaltR::listGeneSet()

allAltGenes <- list('Signac' = as.character(unique(signacAlt$name)), 
                    'ArchR' = as.character(unique(archrAlt$name)), 
                    'MOCHA' = as.character(unique(altTSS$name)))

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

methodPathways <- lapply(allAltGenes[c('Signac','ArchR')], function(x) {
    
                    simplifiedORA(enrichDataBaseList$name[9], x, refList)
    
    })
                           
Hierarch <- read.csv('../../ReactomePathwaysRelation.csv', col.names = c('V1', 'V2'))
ReactID <- read.csv('../../ReactomePathways.csv', col.names = c('V1', 'V2','V3'))
ReactID$V2 <- sub(" $","",ReactID$V2)
                    

PathAnnot_A = lapply(methodPathways[[2]]$description, function(y)
                    unlist(lapply(y, 
                             function(x) findTrunk(x, Hierarch, ReactID))))  %>% unlist()
methodPathways[[2]]$TopLevel = PathAnnot_A  
## No immune pathways
                                  
write.csv(methodPathways[[2]], 'AltTSS_ArchR_Pathways.csv')
                        
mocha_altTSSPathways <- read.csv('../../AltTSS_Pathways.csv')
PathAnnot_M = lapply(mocha_altTSSPathways$description, function(y)
                    unlist(lapply(y, 
                             function(x) findTrunk(x, Hierarch, ReactID))))  %>% unlist()
mocha_altTSSPathways$TopLevel = PathAnnot_M
                                  
allPaths <- rbindlist(
                                  
##Summary over methods
pathwaySum1 <- data.frame( Method = c('Mocha', 'ArchR', 'Signac'),
                          PathwayNumber = c(dim(mocha_altTSSPathways)[1], length(PathAnnot_A), 0),
                          immunePathwayNumber = c(sum(PathAnnot_M == 'Immune System'), 
                                                  sum(PathAnnot_A == 'Immune System'), 0))
                                  
### Let's compare by the number of genes with alternative TSS sites. 

altGenes <- data.frame(AltGeneNumber = c(length(unique(signacAlt$name)), length(unique(archrAlt$name)), length(unique(altTSS$name))),
                        Method = c('Signac', 'ArchR', 'Mocha'))

allAltGenes <- list('Signac' = as.character(unique(signacAlt$name)), 
                    'ArchR' = as.character(unique(archrAlt$name)), 
                    'Mocha' = as.character(unique(altTSS$name)))
d <- ggVennDiagram::process_data(Venn(allAltGenes))
                                  
## Compare with pathway overlap too.                                   
allAltPathways <- list('Signac' = c(), 
                    'ArchR' = as.character(unique(methodPathways[[2]]$description)), 
                    'Mocha' = as.character(unique(mocha_altTSSPathways$description)))
d2 <- ggVennDiagram::process_data(Venn(allAltPathways))
                                  

### Plot all of it
pdf('SuppFig14_AltTSS_otherMethods.pdf')
ggplot(altGenes, aes(x = factor(Method, levels = c('Mocha','ArchR', 'Signac')), y = AltGeneNumber, fill = Method)) + geom_col() + scale_fill_MOCHA() + xlab('Method') + ylab('Number of Genes with Alternative TSS Regulation') +
theme_bw() + theme(legend.position = 'none') + ggtitle('Number of Genes with alternative TSS Regulation by Method')

colorList <- c('#89ACDA', '#F26E65', '#31B34A','white', 'white', 'grey','grey' )
names(colorList) <- c('Mocha','ArchR', 'Signac', 'Signac..ArchR','Signac..ArchR..Mocha', "Signac..Mocha", "ArchR..Mocha")

ggVennDiagram::ggVennDiagram(allAltGenes) + ggtitle('Overlap in Genes Identified with Alternatively Regulated TSSs') +
                  theme(legend.position = 'none') + 
        scale_x_continuous(expand = expansion(mult = .2)) 

ggplot() +
  geom_sf(aes(fill = name), data = venn_region(d)) +
  geom_sf(aes(color = name), data = venn_setedge(d)) +
  geom_sf_text(aes(label = name), data = venn_setlabel(d)) +
  geom_sf_text(aes(label = count), data = venn_region(d)) +
  scale_fill_manual(values = alpha(colorList, .2)) +
  scale_color_manual(values = colorList) + theme_void() + theme(legend.position = 'none')
                                  
ggplot(pathwaySum1, aes(x = factor(Method, levels = c('Mocha','ArchR', 'Signac')), y = PathwayNumber, fill = Method)) + geom_col() + scale_fill_MOCHA() + xlab('Method') + ylab('Number of enriched Pathways from Geneset with Alternative TSS Regulation') +
theme_bw() + theme(legend.position = 'none') + ggtitle("MOCHA's alt TSS's have more enriched Pathways.")
                                  
ggplot(pathwaySum1, aes(x = factor(Method, levels = c('Mocha','ArchR', 'Signac')), y = immunePathwayNumber, fill = Method)) + geom_col() + scale_fill_MOCHA() + xlab('Method') + ylab('Number of enriched Immune Pathways from Geneset with Alternative TSS Regulation') +
theme_bw() + theme(legend.position = 'none')  + ggtitle("Only MOCHA's alt TSS's have enriched Immune System Pathways.")
                                  
ggplot(methodPathways[[2]], aes(x = stringr::str_wrap(description,30), y = enrichmentRatio, fill = -log(FDR))) + coord_flip() + geom_col() + 
       ylab('Number of enriched Immune Pathways from Geneset with Alternative TSS Regulation') +
        theme_bw() + theme(legend.position = 'none') + ggtitle("Enriched Pathways for alternative TSS's from ArchR's DATs")
                                  
ggplot() +
  geom_sf(aes(fill = name), data = venn_region(d2)) +
  geom_sf(aes(color = name), data = venn_setedge(d2)) +
  geom_sf_text(aes(label = name), data = venn_setlabel(d2)) +
  geom_sf_text(aes(label = count), data = venn_region(d2)) +
  scale_fill_manual(values = alpha(colorList, .2)) +
  scale_color_manual(values = colorList) + theme_void() + theme(legend.position = 'none')

dev.off()
                        
####################################################



mochaTmp <-  read.csv('Fig3/generateData/cd16_mocha.csv')
mochaDAT <- MOCHA::StringsToGRanges(mochaTmp$Tile)
mcols(mochaDAT) <- mochaTmp
allArchR <- read.csv('Fig3/generateData/cd16_ArchR.csv')
archrTiles <- paste(
                paste(read.csv('Fig3/generateData/cd16_ArchR.csv')[,2], 
                    read.csv('Fig3/generateData/cd16_ArchR.csv')[,4],
                      sep= ':'),
                    read.csv('Fig3/generateData/cd16_ArchR.csv')[,5],
                    sep = '-')
archrTiles <- archrTiles[allArchR$FDR < 0.05]
allSignac <- read.csv('Fig3/generateData/cd16_signac.csv')
SignacTiles <- sub("-", ":", allSignac$X[allSignac$p_val_adj < 0.05])

common = intersect(intersect(archrTiles, mochaTmp$Tile),SignacTiles)

MACommon = setdiff(intersect(archrTiles, mochaTmp$Tile), SignacTiles)
SACommon = setdiff(intersect(archrTiles, SignacTiles), mochaTmp$Tile)
SMCommon = setdiff(intersect(mochaTmp$Tile, SignacTiles),archrTiles )
unique_archr = setdiff(setdiff(archrTiles, SignacTiles), mochaTmp$Tile)
unique_scmacs = setdiff(mochaTmp$Tile, union(SignacTiles,archrTiles))
unique_signac = setdiff(SignacTiles, union(mochaTmp$Tile, archrTiles))
                                  
annotations1 <- rowRanges(STM)

tmp1 <- lapply(list(common, MACommon, SACommon, SMCommon, unique_archr, unique_scmacs, unique_signac, SignacTiles, mochaTmp$Tile, archrTiles), function(XX){
    
                tmp <- as.data.frame(table(annotations1$tileType[names(annotations1) %in% XX]))
                tmp$Perc = tmp$Freq/sum(tmp$Freq)
                tmp
    
    })
names(tmp1) <- c('Common', 'M-A', 'S-A', 'S-M', 'A only', 'M only', 'S only', 'S', 'M', 'A')
                                  
                                  


#### Support functions

scale_fill_MOCHA<- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('#89ACDA', 
                            '#F26E65', 
                           '#31B34A'
                            ), 
                          c('Mocha','ArchR',
                            'Signac'
                            )),
    )
}