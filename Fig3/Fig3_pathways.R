############################################


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

library(WebGestaltR)

enrichDataBaseList <- WebGestaltR::listGeneSet()

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

geneLists <- lapply(seq_along(allFeatures), function(x){
    
        geneSet <- plyranges::filter_by_overlaps(rowRanges(STM),
                                         MOCHA::StringsToGRanges(allFeatures[[x]])) %>%
                    plyranges::filter(tileType == 'Promoter')
        unique(unlist(lapply(geneSet$Gene, function(y){str_split(y, pattern = ', ')})))
    
    })
names(geneLists) <- gsub("_Mat", "", names(allFeatures))

allORA <- lapply(enrichDataBaseList$name[c(2,7,8,9,10)], function(y){
    tmp <- lapply(geneLists, function(x){
    
        WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = y,
                         interestGene = x, 
                        interestGeneType ="genesymbol",
                        referenceGene =  refList, 
                            referenceGeneType= "genesymbol")
    
    })
    names(tmp) <-  gsub("_Mat", "", names(allFeatures))
    tmp
})

names(allORA) <- enrichDataBaseList$name[c(2,7,8,9,10)]

saveRDS(allORA, 'Figure3_DATs_Pathways.rds')


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
                                