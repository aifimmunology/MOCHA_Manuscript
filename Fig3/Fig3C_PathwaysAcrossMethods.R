## I've removed all subset analysis

# ##

setwd('scMACS_Analysis')


library(ArchR)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(tidyverse)

Signac <- read.csv('cd16_signac.csv')
Signac <- Signac[Signac$p_val_adj < 0.05,]
Mocha <- read.csv('cd16_mocha.csv')
ArchR <- read.csv('cd16_ArchR.csv')

Mocha_features <- MOCHA::StringsToGRanges(Mocha$Tile)
ArchR_features <- MOCHA::StringsToGRanges(paste(ArchR$seqnames,":",
                                                 ArchR$start,"-",ArchR$end, sep = ""))
Signac_features <- sub("-", ":", Signac$X) %>% 
                    MOCHA::StringsToGRanges(.) %>% sort()

allFeatures_GR <- list(Mocha_features, ArchR_features, Signac_features)
allFeatures_GR <- lapply(allFeatures_GR, function(x) { annotateTiles(x, TxDb.Hsapiens.UCSC.hg38.refGene, org.Hs.eg.db)})
names(allFeatures_GR) <- c('Mocha_Mat','ArchR_Mat', 'Signac_Mat')


################################################

### MOCHA analysis

################################################

STM <- readRDS('SampleTileObject_CD16s_EarlyComp.rds')
saveRDS(STM, 'SampleTileObject_CD16s_EarlyComp.rds')

Signac <- read.csv('cd16_signac.csv')
Signac <- Signac[Signac$p_val_adj < 0.05,]
Mocha <- read.csv('cd16_mocha.csv')
ArchR <- read.csv('cd16_ArchR.csv')

Mocha_features <- Mocha$Tile
ArchR_features <- paste(ArchR$seqnames,":",
                                                 ArchR$start,"-",ArchR$end, sep = "")
Signac_features <- sub("-", ":", Signac$X) 

allFeatures <- list(Mocha_features, ArchR_features, Signac_features)
names(allFeatures) <- c('Mocha_Mat','ArchR_Mat', 'Signac_Mat')

allFeatures_GR <- lapply(allFeatures, StringsToGRanges)
allFeatures_GR <- lapply(allFeatures_GR, function(x) { annotateTiles(x, TxDb.Hsapiens.UCSC.hg38.refGene, org.Hs.eg.db)})
names(allFeatures_GR) <- c('Mocha','ArchR', 'Signac')

library(WebGestaltR)

enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

geneLists <- lapply(allFeatures_GR, function(x){
    
                filter(x, tileType == 'Promoter') %>% mcols(.) %>% as.data.frame() %>% dplyr::select(Gene) %>% unlist() %>%
                paste(., collapse = ", ") %>% stringr::str_split(., pattern = ", ") %>% unlist() %>% unique()
    
    })

allORA <- lapply(enrichDataBaseList$name[c(2,7,9,10)], function(y){
    tmp <- lapply(geneLists, function(x){
    
        WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = y,
                         interestGene = x, 
                        interestGeneType ="genesymbol",
                        referenceGene =  refList, 
                            referenceGeneType= "genesymbol")
    
    })
    names(tmp) <- names(allFeatures)
    tmp
})

names(allORA) <- enrichDataBaseList$name[c(2,7,9,10)]

saveRDS(allORA, 'Figure3_DATs_Pathways.rds')

MOCHA_Pathways <- lapply(enrichDataBaseList$name[c(2,7,9,10)], function(y){
                                mocha_path <- allORA[[y]][[1]]
                                mocha_path$Database = y
                                mocha_path}) %>% do.call('rbind',.)

write.csv(MOCHA_Pathways, 'Figure3_MOCHA_Pathways.csv')

ArchR_Pathways <- lapply(enrichDataBaseList$name[c(2,7,9,10)], function(y){
                                archr_path <- allORA[[y]][[2]]
                                archr_path$Database = y
                                archr_path}) %>% do.call('rbind',.)

write.csv(ArchR_Pathways, 'Figure3_ArchR_Pathways.csv')

Signac_Pathways <- lapply(enrichDataBaseList$name[c(2,7,9,10)], function(y){
                                if(!is.null(allORA[[y]][[3]])){
                                    Signac_path <- allORA[[y]][[3]]
                                    Signac_path$Database = y
                                   Signac_path
                                }else{NULL} }) %>% do.call('rbind',.)

write.csv(Signac_Pathways, 'Figure3_Signac_Pathways.csv')

pdf('allEnrichments.pdf')
lapply(seq_along(allORA), function(x){
    
    lapply(seq_along(allORA[[x]]), function(y){

        if(!is.null(allORA[[x]][[y]])){ 
        
            ggplot(allORA[[x]][[y]], aes(x = FDR, y = str_wrap(description, 20))) + geom_col() + 
            ggtitle(paste(gsub("_Mat","",names(allFeatures)[y]), 
                       gsub('pathway_|_noRedundant','', 
                            names(allORA)[x]), sep = ' '))
        }
        
        })
  
    })

lapply(seq_along(allORA), function(x){
 
    names(allORA[[x]]) <- names(geneLists)
    tmp_enriched <- lapply(allORA[[x]][!unlist(lapply(allORA[[x]], is.null))], function(x){
        
                        x$description
        
        })
    tryCatch(
        {
             ggVennDiagram::ggVennDiagram(tmp_enriched) + ggtitle(names(allORA)[x])
        },  error=function(cond) { return(NULL)}
            )


})

dev.off()

#### Analyze pathway and gene sets across methods
names(geneLists) <- gsub("_Mat", "", names(allFeatures))

geneDF <- as.data.frame(lengths(geneLists)) %>% dplyr::mutate(Method = rownames(.))  %>%
                               dplyr::rename('Genes' = 'lengths(geneLists)')

geneDF$Method <- factor(geneDF$Method, levels = c('Mocha', 'ArchR', 'Signac'))           


names(allORA[[4]]) <- names(geneLists)
pathDF_List <- lapply(allORA[[4]], function(x) dim(x)[1])
pathDF_List[unlist(lapply(pathDF_List, is.null))] = 0 
pathDF <-  data.frame(Pathways = unlist(pathDF_List), Method = names(pathDF_List))

pathDF$Method <- factor(pathDF$Method, levels = c('Mocha', 'ArchR', 'Signac'))           

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

scale_fill_MOCHA2 <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(c('#89ACDA', 
                            '#F26E65', 
                            '#31B34A',
                            '#2BB892',
                            '#11B5EA',
                            '#9589C1',
                            '#D66DAB'
                            ), 
                          c('Mocha_Unique','ArchR_Unique',
                            'Signac_Unique', 'Mocha_ArchR',
                            'Mocha_Signac', 'ArchR_Signac',
                            'All_Methods'
                            )),
    )
}



pdf('DATs_GeneLists & Pathways.pdf')

ggplot(geneDF, 
       aes(x = factor(Method, levels = orderList),
           y = Genes, fill = Method)) + geom_col() +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                               ylab('Number of Unique Genes (Promoters)') + scale_fill_MOCHA()

ggplot(pathDF, 
       aes(x = factor(Method, levels = orderList),
           y = Pathways, fill = Method)) + geom_col() +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                               ylab('Number of Enriched Reactome Pathways') +
                           scale_fill_MOCHA() 

dev.off()