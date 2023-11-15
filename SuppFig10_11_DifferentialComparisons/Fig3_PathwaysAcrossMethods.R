    ]]]]]]]library(ArchR)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(tidyverse)
library(WebGestaltR)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(pbapply)

enrichDataBaseList <- WebGestaltR::listGeneSet()
allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- AnnotationDbi::mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

deseq <- read.csv("DATs_DESeq2_EarlyInfection.csv")
deseq <- (na.omit(deseq))
deseq <- deseq[deseq$padj < 0.05, ]
diffchipl <- read.csv("DATS_DiffChipL_Differentials.csv")
diffchipl <- diffchipl[diffchipl$adj.P.Val < 0.05, ]
Signac <- read.csv("cd16_signac.csv")
Mocha <- read.csv("cd16_mocha.csv")
ArchR <- read.csv("cd16_ArchR.csv")

Mocha_features <- paste0(Mocha$seqnames, ":", Mocha$start, "-", Mocha$end)
ArchR_features <- paste0(ArchR$seqnames, ":", ArchR$start, "-", ArchR$end)
Signac_features <- paste0(Signac$seqnames, ":", Signac$start, "-", Signac$end)
deseq_features <- paste0(deseq$seqnames, ":", deseq$start, "-", deseq$end)
diffchipl_features <- diffchipl$X

allFeatures <- list(Mocha_features, ArchR_features, Signac_features, deseq_features, diffchipl_features)

allFeatures_GR <- lapply(allFeatures, MOCHA::StringsToGRanges)
names(allFeatures_GR) <- c("MOCHA", "ArchR", "Signac", "DESeq2", "DiffChipL")

allFeatures_GR <- lapply(allFeatures_GR, function(x) {
  MOCHA::annotateTiles(x, TxDb.Hsapiens.UCSC.hg38.refGene, org.Hs.eg.db)
})


# ###############################################

# ## MOCHA analysis

# ###############################################

geneLists <- lapply(allFeatures_GR, function(x) {
  plyranges::filter(x, tileType == "Promoter") %>%
    mcols(.) %>%
    as.data.frame() %>%
    dplyr::select(Gene) %>%
    unlist() %>%
    paste(., collapse = ", ") %>%
    stringr::str_split(., pattern = ", ") %>%
    unlist() %>%
    unique()
})

# dblist <- c(2, 7, 9, 10)
dblist <- c(9) # only pathway_Reactome
# allORA <- lapply(enrichDataBaseList$name[c(2, 7, 9, 10)], function(y) {
allORA <- lapply(enrichDataBaseList$name[dblist], function(y) {
  print(y)
  tmp <- lapply(geneLists, function(x) {
    WebGestaltR::WebGestaltR(
      enrichMethods = "ORA", organism = "hsapiens",
      enrichDatabase = y,
      interestGene = x,
      interestGeneType = "genesymbol",
      referenceGene = refList,
      referenceGeneType = "genesymbol"
    )
  })
  names(tmp) <- names(allFeatures)
  tmp
})

names(allORA) <- enrichDataBaseList$name[dblist]

saveRDS(allORA, "Figure3_DATs_Pathways.rds")

# MOCHA_Pathways <- lapply(enrichDataBaseList$name[dblist], function(y) {
#   mocha_path <- allORA[[y]][[1]]
#   mocha_path$Database <- y
#   mocha_path
# }) %>% do.call("rbind", .)

# write.csv(MOCHA_Pathways, "Figure3_MOCHA_Pathways.csv")

# ArchR_Pathways <- lapply(enrichDataBaseList$name[dblist], function(y) {
#   archr_path <- allORA[[y]][[2]]
#   archr_path$Database <- y
#   archr_path
# }) %>% do.call("rbind", .)

# write.csv(ArchR_Pathways, "Figure3_ArchR_Pathways.csv")

# Signac_Pathways <- lapply(enrichDataBaseList$name[dblist], function(y) {
#   if (!is.null(allORA[[y]][[3]])) {
#     Signac_path <- allORA[[y]][[3]]
#     Signac_path$Database <- y
#     Signac_path
#   } else {
#     NULL
#   }
# }) %>% do.call("rbind", .)

# write.csv(Signac_Pathways, "Figure3_Signac_Pathways.csv")

# pdf("allEnrichments.pdf")
# lapply(seq_along(allORA), function(x) {
#   lapply(seq_along(allORA[[x]]), function(y) {
#     if (!is.null(allORA[[x]][[y]])) {
#       ggplot(allORA[[x]][[y]], aes(x = FDR, y = str_wrap(description, 20))) +
#         geom_col() +
#         ggtitle(paste(gsub("_Mat", "", names(allFeatures)[y]),
#           gsub(
#             "pathway_|_noRedundant", "",
#             names(allORA)[x]
#           ),
#           sep = " "
#         ))
#     }
#   })
# })

# lapply(seq_along(allORA), function(x) {
#   names(allORA[[x]]) <- names(geneLists)
#   tmp_enriched <- lapply(allORA[[x]][!unlist(lapply(allORA[[x]], is.null))], function(x) {
#     x$description
#   })
#   tryCatch(
#     {
#       ggVennDiagram::ggVennDiagram(tmp_enriched) + ggtitle(names(allORA)[x])
#     },
#     error = function(cond) {
#       return(NULL)
#     }
#   )
# })

# dev.off()

# #### Analyze pathway and gene sets across methods
# names(geneLists) <- gsub("_Mat", "", names(allFeatures))

# geneDF <- as.data.frame(lengths(geneLists)) %>%
#   dplyr::mutate(Method = rownames(.)) %>%
#   dplyr::rename("Genes" = "lengths(geneLists)")

# geneDF$Method <- factor(geneDF$Method, levels = c("Mocha", "ArchR", "Signac"))


# names(allORA[[4]]) <- names(geneLists)
# pathDF_List <- lapply(allORA[[4]], function(x) dim(x)[1])
# pathDF_List[unlist(lapply(pathDF_List, is.null))] <- 0
# pathDF <- data.frame(Pathways = unlist(pathDF_List), Method = names(pathDF_List))

# pathDF$Method <- factor(pathDF$Method, levels = c("Mocha", "ArchR", "Signac"))

# scale_fill_MOCHA <- function(...) {
#   ggplot2:::manual_scale(
#     "fill",
#     values = setNames(
#       c(
#         "#89ACDA",
#         "#F26E65",
#         "#31B34A"
#       ),
#       c(
#         "Mocha", "ArchR",
#         "Signac"
#       )
#     ),
#   )
# }

# scale_fill_MOCHA2 <- function(...) {
#   ggplot2:::manual_scale(
#     "fill",
#     values = setNames(
#       c(
#         "#89ACDA",
#         "#F26E65",
#         "#31B34A",
#         "#2BB892",
#         "#11B5EA",
#         "#9589C1",
#         "#D66DAB"
#       ),
#       c(
#         "Mocha_Unique", "ArchR_Unique",
#         "Signac_Unique", "Mocha_ArchR",
#         "Mocha_Signac", "ArchR_Signac",
#         "All_Methods"
#       )
#     ),9
#   )
# }



# pdf("DATs_GeneLists & Pathways.pdf")

# ggplot(
#   geneDF,
#   aes(
#     x = factor(Method, levels = orderList),
#     y = Genes, fill = Method
#   )
# ) +
#   geom_col() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   ylab("Number of Unique Genes (Promoters)") +
#   scale_fill_MOCHA()

# ggplot(
#   pathDF,
#   aes(
#     x = factor(Method, levels = orderList),
#     y = Pathways, fill = Method
#   )
# ) +
#   geom_col() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   ylab("Number of Enriched Reactome Pathways") +
#   scale_fill_MOCHA()

# dev.off()
