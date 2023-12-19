# #####################################################
library(parallel)
library(dplyr)
ncores <- 4
STM <- readRDS('CD16_SampleTileObj.rds')

# Randomly select 8 DATs, do kmeans clustering across samples within STAM, measure accuracy. 
# Repeat 1000x, then generate violin plots for each method

Signac <- read.csv('cd16_signac.csv')
Mocha <- read.csv('cd16_mocha.csv')
ArchR <- read.csv('cd16_ArchR.csv')

deseq <- read.csv("DATs_DESeq2_EarlyInfection.csv")
dim(deseq)
deseq <- (na.omit(deseq))
dim(deseq)
deseq <- deseq[deseq$padj < 0.05,]
dim(deseq)

diffchipl <- read.csv("DATS_DiffChipL_Differentials.csv")
dim(diffchipl)
diffchipl <- diffchipl[diffchipl$adj.P.Val < 0.05,]
dim(diffchipl)


Mocha_features <- paste0(Mocha$seqnames, ":", Mocha$start, "-", Mocha$end)
ArchR_features <- paste0(ArchR$seqnames, ":", ArchR$start, "-", ArchR$end)
Signac_features <- paste0(Signac$seqnames, ":", Signac$start, "-", Signac$end)

deseq_features <- paste0(deseq$seqnames, ":", deseq$start, "-", deseq$end)
diffchipl_features <- diffchipl$X

allFeatures <- list(Mocha_features, ArchR_features, Signac_features, deseq_features, diffchipl_features)
names(allFeatures) <- c('MOCHA','ArchR', 'Signac', 'DESeq2','DiffChipL')

somefeatures <- list(Mocha_features, deseq_features, diffchipl_features)
names(somefeatures) <- c('Mocha_Mat', 'DESeq2_Mat','DiffChipL_Mat')

# 2D Venn diagram
library(ggVennDiagram)
pdf("venn.pdf")
ggVennDiagram(somefeatures)
dev.off()

# Upset plot
library(UpSetR)
pdf("upset_daps.pdf")
upset(fromList(allFeatures))
dev.off()

# Different Upset plots 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("krassowski/complex-upset")
library(ggplot2)
library(ComplexUpset)
df <- data.frame(tileID = unique(unlist(allFeatures)))
rownames(df) <- df$tileID
df$MOCHA <- df$tileID %in% allFeatures$MOCHA
df$ArchR <- df$tileID %in% allFeatures$ArchR
df$Signac <- df$tileID %in% allFeatures$Signac
df$DESeq2 <- df$tileID %in% allFeatures$DESeq2
df$DiffChipL <- df$tileID %in% allFeatures$DiffChipL
df <- df[c('MOCHA','ArchR', 'Signac', 'DESeq2','DiffChipL')]

pdf("upset_daps_nocolor.pdf")
upset(
    df, 
    c('MOCHA','ArchR', 'Signac', 'DESeq2','DiffChipL'), 
    name='Method',
    base_annotations=list(
        'Intersection size'=intersection_size(
            text=list(
                size = 2.7
            )
        )
    ),
    set_sizes=(
        upset_set_size()
        + theme(axis.text.x=element_text(angle=90))
        + ylab("# DAPs")
    ),
    # queries=list(
    #     upset_query(set='MOCHA', color='blue', fill='blue'),
    #     upset_query(set='ArchR', color='red', fill='red'),
    #     upset_query(set='Signac', color='pink', fill='pink'),
    #     upset_query(set='DESeq2', color='yellow', fill='yellow'),
    #     upset_query(set='DiffChipL', color='brown', fill='brown'),
    #     upset_query(
    #         intersect=c('MOCHA', 'ArchR'),
    #         color='purple',
    #         fill='purple',
    #         only_components=c('intersections_matrix', 'Intersection size')
    #     ),
    #     upset_query(
    #         intersect=c('DESeq2', 'DiffChipL'),
    #         color='orange',
    #         fill='orange',
    #         only_components=c('intersections_matrix', 'Intersection size')
    #     ),
    #     upset_query(
    #         intersect=c('MOCHA', 'DESeq2'),
    #         color='green',
    #         fill='green',
    #         only_components=c('intersections_matrix', 'Intersection size')
    #     )
    # ),
    min_size=9,
    width_ratio=0.1
)
dev.off()


## Extract tile matrix.
tileMat <- assays(STM)[[1]]
tileMat[is.na(tileMat)] = 0
meta <- colData(STM)

## Calculate Holly's G
findConcordance <- function(meta, cluster){
    
     cM <- ArchR::confusionMatrix(meta$COVID_status[
                 match(names(cluster), rownames(meta))], cluster) 

     s1 = cM[1,1] + cM[2,2]
     s2 = cM[1,2] + cM[2,1]
     s = abs(s1-s2)
     s/(s1+s2)
    return(s/(s1+s2))
}

k <- 50
ncores <- 4
for (k in c(25, 50, 100)) {
    message("k : ", k)
    subClustListsMocha <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[1]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)
    
    
    
    subClustListsArchR <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[2]],k),])
                    kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta,kClust)
        }, mc.cores = ncores)
    
    
    subClustListsSignac <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[3]],k),])
                    kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta,kClust)
        }, mc.cores = ncores)  
    
    subClustListsDESeq <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[4]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)
    
    subClustListsDiffChipL <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[5]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)
    
    
    
    
    # MUnique <- setdiff(allFeatures[[1]], c(allFeatures[[2]],allFeatures[[3]], allFeatures[[4]], allFeatures[[5]]))
    # subClustMocha_Unique <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(MUnique,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores)  
    
    # AUnique <- setdiff(allFeatures[[2]], c(allFeatures[[1]],allFeatures[[3]], allFeatures[[4]], allFeatures[[5]]))
    # subClustArchR_Unique <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(AUnique,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores)
    
    # DESeqUnique <- setdiff(allFeatures[[4]], c(allFeatures[[1]],allFeatures[[2]], allFeatures[[3]], allFeatures[[5]]))
    # print(length(DESeqUnique))
    # subClustDESeq_Unique <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(DESeqUnique,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores)
    
    # DiffChipLUnique <- setdiff(allFeatures[[5]], c(allFeatures[[1]],allFeatures[[2]], allFeatures[[3]], allFeatures[[4]]))
    # print(length(DiffChipLUnique))
    # subClustDiffChipL_Unique <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(DiffChipLUnique,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores)  
      
    # All_u <- intersect(allFeatures[[1]],
    #                    intersect(allFeatures[[2]],
    #                              intersect(allFeatures[[3]],
    #                                        intersect(allFeatures[[4]], allFeatures[[5]]))))
    
    # All_u <- intersect(allFeatures[[1]], intersect(allFeatures[[2]],intersect(allFeatures[[3]])  )
    # subClustAll <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(All_u,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores)
    
    # MA_u <- setdiff(intersect(allFeatures[[1]],allFeatures[[2]]), All_u)   
    # subClustMochaArchR <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(MA_u,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores) 
    
    # SA_u <-  setdiff(intersect(allFeatures[[2]],allFeatures[[3]]),   All_u)
    # subClustSignacArchR <- mclapply(1:1000, function(x){
    #                 subMat <- t(tileMat[rownames(tileMat) %in% sample(SA_u,k),])
    #                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
    #                 findConcordance(meta, kClust)
    #     }, mc.cores = ncores) 
    
    # methods <- c(
    #             'Mocha', 'ArchR', 'Signac',
    #             'DESeq2', 'DiffChipL', 
    #             'Mocha_Unique','ArchR_Unique',
    #             'Mocha_ArchR', 'Signac_ArchR', 
    #             'DESeq2_Unique', 'DiffChipL_Unique', 
    #             'All_Methods')
                 
    # allGs <- data.frame(
    #     Mocha = unlist(subClustListsMocha),
    #     ArchR = unlist(subClustListsArchR),
    #     Signac = unlist(subClustListsSignac),
    #     DESeq2 = unlist(subClustListsDESeq),
    #     DiffChipL = unlist(subClustListsDiffChipL),
    #     Mocha_Unique = unlist(subClustMocha_Unique),
    #     ArchR_Unique = unlist(subClustArchR_Unique), 
    #     Mocha_ArchR = unlist(subClustMochaArchR),
    #     Signac_ArchR = unlist(subClustSignacArchR),
    #     DESeq2_Unique = unlist(subClustDESeq_Unique),
    #     DiffChipL_Unique = unlist(subClustDiffChipL_Unique),
    #     All_Methods = unlist(subClustAll)) %>% 
    #         tidyr::pivot_longer(cols = methods,
    #                             names_to = 'Method', values_to = 'G')
    
    methods <- c('Mocha', 'ArchR', 'Signac',
                'DESeq2', 'DiffChipL')
    
                 
    allGs <- data.frame(
        Mocha = unlist(subClustListsMocha),
        ArchR = unlist(subClustListsArchR),
        Signac = unlist(subClustListsSignac),
        DESeq2 = unlist(subClustListsDESeq),
        DiffChipL = unlist(subClustListsDiffChipL)
    ) %>% tidyr::pivot_longer(cols = methods, names_to = 'Method', values_to = 'G')
        
    
    allGs$Method <- factor(allGs$Method, levels = methods)
    
    write.csv(allGs, stringr::str_interp("Fig3_HolleysG_${k}.csv"))
}

# wilcox.test(filter(allGs, Method == 'ArchR')$G, filter(allGs, Method == 'Mocha')$G) # p-value < 2.2e-16
# wilcox.test(filter(allGs, Method == 'Signac')$G, filter(allGs, Method == 'Mocha')$G) # p-value < 2.2e-16
# wilcox.test(filter(allGs, Method == 'Signac')$G, filter(allGs, Method == 'ArchR')$G) # p-value = 0.01809                        

# wilcox.test(filter(allGs, Method == 'DESeq2')$G, filter(allGs, Method == 'DiffChipL')$G)
# W = 517004, p-value = 0.1825
# alternative hypothesis: true location shift is not equal to 0

allG_unclean <- read.csv('Fig3_HolleysG.csv')
allGs <- allG_unclean %>% dplyr::mutate(Method=recode(Method, All_Methods="Signac_MOCHA_ArchR", Mocha="MOCHA", Mocha_Unique="MOCHA_Unique", Mocha_ArchR="MOCHA_ArchR"))

median_ranking <- allGs %>% 
  group_by(Method) %>% 
  summarise(median = median(G, na.rm = TRUE)) %>% 
  arrange(desc(median))

allGs <- allGs %>% mutate(Method=factor(Method,levels=median_ranking$Method))

library(ggpubr)
pdf("PredictiveValue_ByMethod_all_v2.pdf")
ggplot(allGs, aes(x = Method, y = G, fill = Method)) + 
  geom_violin() + 
  theme_bw() +
  ggtitle("Performance of 50 randomly sub-sampled Differential Regions by Method") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  ylab("Absolute of Holley's G") +
  scale_fill_manual(values = list(
    'ArchR' = '#F26E65', 'MOCHA' =  '#89ACDA', 'Signac' = '#31B34A',
    'MOCHA_Unique' = '#89ACDA', 'ArchR_Unique' = '#F26E65',
    'MOCHA_ArchR' = '#C14282', 'Signac_ArchR' = '#C0813F',
    'DESeq2' = '#8E44AD', 'DiffChipL' = '#F4D03F',
    'DESeq2_Unique' = '#8E44AD', 'DiffChipL_Unique' = '#F4D03F', 
    'Signac_MOCHA_ArchR' = '#9E5F3F'))

dev.off()                              



k <- 50
ncores <- 4
allks <- c(25, 50, 100)
allks <- c(75)
for (k in allks) {
    message("k : ", k)
    subClustListsMocha <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[1]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)
    
    
    
    subClustListsArchR <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[2]],k),])
                    kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta,kClust)
        }, mc.cores = ncores)
    
    
    subClustListsSignac <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[3]],k),])
                    kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta,kClust)
        }, mc.cores = ncores)  
    
    subClustListsDESeq <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[4]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)
    
    subClustListsDiffChipL <- mclapply(1:1000, function(x){
                    subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[5]],k),])
                     kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                    findConcordance(meta, kClust)
        }, mc.cores = ncores)

    
    methods <- c('Mocha', 'ArchR', 'Signac',
                'DESeq2', 'DiffChipL')
    
                 
    allGs <- data.frame(
        Mocha = unlist(subClustListsMocha),
        ArchR = unlist(subClustListsArchR),
        Signac = unlist(subClustListsSignac),
        DESeq2 = unlist(subClustListsDESeq),
        DiffChipL = unlist(subClustListsDiffChipL)
    ) %>% tidyr::pivot_longer(cols = methods, names_to = 'Method', values_to = 'G')
        
    
    allGs$Method <- factor(allGs$Method, levels = methods)
    
    write.csv(allGs, stringr::str_interp("Fig3_HolleysG_${k}.csv"))
}

allks <- c(25, 50, 75, 100)
allks <- c(75)
for (k in allks){
    allGs_k <- read.csv(stringr::str_interp("Fig3_HolleysG_${k}.csv"))

    allGs_k <- allGs_k %>% dplyr::mutate(Method=recode(Method, Mocha="MOCHA"))
  
    median_ranking <- allGs_k %>% 
      group_by(Method) %>% 
      summarise(median = median(G, na.rm = TRUE)) %>% 
      arrange(desc(median))
    
    allGs_k <- allGs_k %>% mutate(Method=factor(Method,levels=median_ranking$Method))
  
    library(ggpubr)
    pdf(stringr::str_interp("PredictiveValue_ByMethod_all_${k}.pdf"))
    ggplot(allGs_k, aes(x = Method, y = G, fill = Method)) + geom_violin() + theme_bw() +
                ggtitle(stringr::str_interp("Performance of ${k} randomly sub-sampled Differential Regions by Method")) +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
                ylab("Absolute of Holley's G") +
                scale_fill_manual(
                    values = list(
                        'ArchR' = '#F26E65', 'MOCHA' =  '#89ACDA', 'Signac' = '#31B34A',
                        'DESeq2' = '#8E44AD', 'DiffChipL' = '#F4D03F'))
    
    dev.off()
}
