# #####################################################

STM <- readRDS('CD16_SampleTileObj.rds')

# Randomly select 8 DATs, do kmeans clustering across samples within STAM, measure accuracy. 
# Repeat 1000x, then generate violin plots for each method

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

## Extract tile matrix.
tileMat <- assays(STM)[[1]]
tileMat[is.na(tileMat)] = 0
meta <- colData(STM)

##Calculate Holly's G
findConcordance <- function(meta, cluster){
    
     cM <- ArchR::confusionMatrix(meta$COVID_status[
                 match(names(cluster), rownames(meta))], cluster) 

     s1 = cM[1,1] + cM[2,2]
     s2 = cM[1,2] + cM[2,1]
     s = abs(s1-s2)
     s/(s1+s2)
    return(s/(s1+s2))
}


subClustListsMocha <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[1]],50),])
                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30)



subClustListsArchR <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[2]],50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta,kClust)
    }, mc.cores = 30)


subClustListsSignac <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[3]],50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta,kClust)
    }, mc.cores = 30)  

MUnique <- setdiff(allFeatures[[1]], c(allFeatures[[2]],allFeatures[[3]]))

subClustMocha_Unique <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(MUnique,50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30)  

AUnique <- setdiff(allFeatures[[2]], c(allFeatures[[1]],allFeatures[[3]]))                               

subClustArchR_Unique <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(AUnique,50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30)

MA_u <- intersect(allFeatures[[1]],allFeatures[[2]])    

subClustMochaArchR <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(MA_u,50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30) 

SA_u <- intersect(allFeatures[[2]],allFeatures[[3]])    

subClustSignacArchR <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(SA_u,50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30) 

All_u <- intersect(allFeatures[[1]], intersect(allFeatures[[2]],allFeatures[[3]])  )

subClustAll <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(All_u,50),])
                kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 30) 

allGs <- data.frame(Mocha = unlist(subClustListsMocha), ArchR = unlist(subClustListsArchR), 
                          Signac = unlist(subClustListsSignac), 
                      Mocha_Unique = unlist(subClustMocha_Unique),
                     ArchR_Unique = unlist(subClustArchR_Unique), 
                      Mocha_ArchR = unlist(subClustMochaArchR),
                     Signac_ArchR = unlist(subClustSignacArchR),
                     All_Methods = unlist(subClustAll)) %>%
                tidyr::pivot_longer(cols = c('Mocha', 'ArchR', 'Signac', 
                                            'Mocha_Unique','ArchR_Unique', 'Mocha_ArchR', 'Signac_ArchR',
                                            'All_Methods'),
                                    names_to = 'Method', values_to = 'G')

allGs$Method <- factor(allGs$Method, levels = c('Mocha', 'ArchR', 'Signac',
                                                    'Mocha_Unique', 'ArchR_Unique',
                                                    'Mocha_ArchR', 'Signac_ArchR',
                                                    'All_Methods'))

write.csv(allGs, 'Fig3_HolleysG.csv')


wilcox.test(filter(allGs, Method == 'ArchR')$G, filter(allGs, Method == 'Mocha')$G) # p-value < 2.2e-16
wilcox.test(filter(allGs, Method == 'Signac')$G, filter(allGs, Method == 'Mocha')$G) # p-value < 2.2e-16
wilcox.test(filter(allGs, Method == 'Signac')$G, filter(allGs, Method == 'ArchR')$G) # p-value = 0.01809

allDiffs <- append(allFeatures, list(MUnique, AUnique, MA_u, SA_u, All_u))
names(allDiffs) <- c('Mocha', 'ArchR', 'Signac', 'Mocha_Unique', 'ArchR_Unique', 'Mocha_ArchR',
                          'Signac_ArchR', 'All_Methods')
                   
mColsAll_Perc <- lapply(seq_along(allDiffs), function(x){
    
                subList <- mcols(plyranges::filter_by_overlaps(rowRanges(STM), 
                                                               MOCHA::StringsToGRanges(allDiffs[[x]]))) %>% 
                    as.data.frame()
                group_by(subList, tileType) %>% summarize(Percentage = dplyr::n()/dim(subList)[1]) %>%
                mutate(type = names(allDiffs)[x])
    })
names(mColsAll_Perc) <- c('Mocha', 'ArchR', 'Signac', 'Mocha_Unique', 'ArchR_Unique', 'Mocha_ArchR',
                          'Signac_ArchR', 'All_Methods')

percentBreakDown <- do.call('rbind', mColsAll_Perc)  %>% 
                        pivot_wider(id_cols = 'type', names_from = 'tileType', values_from = 'Percentage') %>%
                             dplyr::rename('Method' = 'type') 

## Code used from here:
# https://gist.github.com/tomsing1/6df5e3a9d98617eac85892b3d8056a59       
clr <- function(m) {
  apply(m, 2, function(x) {
    log1p(x = x/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x))))
  })
}

percentBreakDown2 <- as.data.frame(clr(percentBreakDown[,c(2:4)]))
percentBreakDown2$Method = percentBreakDown$Method                       

sumAllGs <- allGs %>% group_by(Method) %>% summarize(MedG = median(G)) %>%  
                             left_join(percentBreakDown2, by = 'Method') %>% 
                             filter(!is.na(Distal))

library(ggpubr)


pdf('PredictiveValue_ByMethod.pdf')
ggplot(allGs, aes(x = Method, y = G, fill = Method)) + geom_violin() + theme_bw() +
            ggtitle('Performance of 50 randomly sub-sampled Differential Regions by Method') +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
            ylab("Absolute of Holley's G") +
            scale_fill_manual(values = 
                              list('ArchR' = '#F26E65', 'Mocha' =  '#89ACDA', 'Signac' = '#31B34A',
                                   'Mocha_Unique' = '#89ACDA', 'ArchR_Unique' = '#F26E65', 
                                  'Mocha_ArchR' = '#C14282', 'Signac_ArchR' = '#C0813F',
                                  'All_Methods' = '#9E5F3F'))

ggplot(sumAllGs, aes(x = Promoter, y = MedG)) + geom_point() + theme_bw() + 
                              geom_smooth(method='lm',  show.legend = TRUE) +
                             stat_cor(method = 'spearman') +
                             ylab("Median value of Holley's G across iterations") +
                             xlab("Percent Promoter Tiles Detected") +
                             ggtitle('Percent Promoter within DATs correlated to Sample Discrimination')
ggplot(sumAllGs, aes(x = Intragenic, y = MedG)) + geom_point() + theme_bw()+ 
                              geom_smooth(method='lm') + 
                             stat_cor(method = 'spearman') +
                             ylab("Median value of Holley's G across iterations") +
                             xlab("Percent Intragenic Tiles Detected") +
                             ggtitle('Percent Intragenic Tiles within DATs anti-correlated to Sample Discrimination')
ggplot(sumAllGs, aes(x = Distal, y = MedG)) + geom_point() + theme_bw()+ 
                              geom_smooth(method='lm') +
                             stat_cor(method = 'spearman')+
                              ylab("Median value of Holley's G across iterations") +
                             xlab("Percent Distal Tiles Detected") +
                             ggtitle('Percent Promoter within DATs trends with Sample Discrimination')

dev.off()                              
