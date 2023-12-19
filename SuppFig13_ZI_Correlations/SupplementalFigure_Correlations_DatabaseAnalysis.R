library(MOCHA)
library(SummarizedExperiment)
library(tidyverse)
library(plyranges)

## Load in database links and TSAM object
STM <- readRDS('SampleTileMatrix_AllCellTypes.rds')
databaseLinks <- read.csv('ActivePromoterEnhancerLinks_Hg38.csv')

##Filter to just look at CD8 and CD4 Naive links. 
matchedLinks <- databaseLinks[grepl('nCD',databaseLinks$cellTypes),]
## Filter to look for only peaks detected in CD4 and CD8 Naives
matchedTiles <- rowRanges(STM)[rowRanges(STM)$'CD8 Naive' | rowRanges(STM)$'CD4 Naive']

## Tile the links so that we look at all the links 
tiledLinks <- matchInteractions(NULL, matchedLinks, returnTileLinks = TRUE)
fTiledLinks <- tiledLinks[which(tiledLinks$Tile1 %in% GRangesToString(matchedTiles) &
      tiledLinks$Tile2 %in% GRangesToString(matchedTiles)),]


fullObj <- MOCHA:::combineSampleTileMatrix(STM)
  cl <- parallel::makeCluster(35)
accMat <- SummarizedExperiment::assays(fullObj)[[1]]
combPairs <- data.frame(fTiledLinks$Tile1, fTiledLinks$Tile2)

combinedMat <- cbind(accMat[combPairs[c(1:100), 1],], accMat[combPairs[c(1:100), 2],])

matList <- as.list(as.data.frame(t(combinedMat)))


testVar <- runCoAccessibility(accMat, combPairs[1:100,], ZI = TRUE, numCores = cl)

       
  sampleNum <- dim(subAccMat)[2]

parallel::stopCluster(cl)
cl <- parallel::makeCluster(35)
#
       # getSpearman <- function(x){}
testVar <- runCoAccessibility(accMat, combPairs[1:100,], numCores = cl)

getSpearman = function(x){}

runCoAccessibility <- function(accMat, pairs, ZI = TRUE, verbose = TRUE, numCores = 1) {

  # Generate matrix for just Tile1, and just Tile2, then combined by column. 
  combinedMat <- cbind(accMat[pairs[, 1],], accMat[pairs[, 2],])
  matList <- as.list(as.data.frame(t(combinedMat)))
  #rm(combinedMat)
  #rm(accMat)
  message('Now exporting...')

  parallel::clusterExport(numCores, varlist = c('matList'), envir = environment())

  #Remember how many samples you have, so you can split the list later
  sampleNum <- dim(accMat)[2]
  if(verbose){
    message('Data.frame wrangled for co-accessibility')
  }

  ## split combinedMat into blocks that represent the number of cores we are using. 
  ## Use split command for that. 
  ## Then test each pair, but only within each block. 

 
 # parallel::clusterExport(numCores, varlist = c('getSpearman'), envir = environment())
  #message('SpearmanMade')
     #newEnv <- new.env(parent = emptyenv())

  #parallel::clusterExport(cl, 
   #   varlist = c('getSpearman', 'ZI', 'sampleNum','verbose'), 
   #    envir = environment())
  newFunctioNow <- function(x){}
  message('VariablesExported')

  tmpFun <- function(x, cl){
      
    zero_inflated_spearman <- unlist(pbapply::pblapply(X=x,
        FUN = newZI,
        cl = numCores
      ))
  }
 zero_inflated_spearman <- tmpFun(matList, cl)
 return(zero_inflated_spearman)
    
}

newZI <- function(x, verbose = FALSE, ZI = TRUE){
    
    tryCatch(
        {
            numLen <- length(x)/2
         MOCHA:::weightedZISpearman(x = x[1:numLen] ,y=x[(numLen+1):(2*numLen)], 
                                    verbose , ZI)
           },
            error = function(e) {
              NA
            }
          )
    
}

newZI(x=matList[[1]])
numLen = 1911
MOCHA:::weightedZISpearman(x = matList[[1]][1:numLen] ,y=matList[[1]][(numLen+1):(2*numLen)], verbose =TRUE, ZI = TRUE)


 assign("getSpearman", function(x) {
          tryCatch(
            {
              MOCHA:::weightedZISpearman(
                x = x[1:1911],
                y = x[(1911+1):(2*1911)],
                verbose = TRUE,
                ZI = TRUE
              )
            },
            error = function(e) {
              NA
            }
          )
        },  envir = globalenv())
  

runSpearman <- function(matList, cl){
    
      newMat <- lapply(matList, function(z) z*10)
                        
      getSpearman <- function(z){}
                        
      parallel::clusterExport(cl,varlist = c('getSpearman'), env = environment())

      unlist(pbapply::pblapply(X = matList,
                    FUN = getSpearman,
                        cl = cl))

  }
  runSpearman(matList, cl)

## Run ZI-Spearman correlations (will return both the ZI-Spearman results, and background peakset. 
ZI_Corr <- MOCHA::testCoAccessibilityRandom(STM,fTiledLinks$Tile1, fTiledLinks$Tile2, 
                    numCores = 25,
                    ZI = TRUE, backNumber = 1000, returnBackGround = TRUE)
saveRDS(ZI_Corr, 'ZI_Correlations_w_backgroundCorrelations.rds')
write.csv(ZI_Corr[[1]], 'Correlations_PromoterEnhancers_ZISpearman.csv')
write.csv(ZI_Corr[[2]], 'Correlations_PromoterEnhancers_ZISpearman_Background.csv') 
## Run a Spearman-based CoAccessibility using the same backgrounds
startTime <- Sys.time()
Norm_Corr <- MOCHA::testCoAccessibilityRandom(STM,fTiledLinks$Tile1, fTiledLinks$Tile2, 
                                       numCores = 45,
                    ZI = FALSE, backNumber = ZI_Corr[['Background']], returnBackGround = TRUE)
stopTime <- Sys.time() # 8.8646 minutes
saveRDS(Norm_Corr, 'Spearman_Correlations.rds')
write.csv(Norm_Corr[[1]], 'Correlations_PromoterEnhancers_Spearman.csv')
write.csv(Norm_Corr[[2]], 'Correlations_PromoterEnhancers_Spearman_Background.csv')       

par(mfrow=c(3,2))
pdf('Coaccc_Foreground_Background.pdf')
hist(ZI_Corr[[1]]$pValues)
hist(Norm_Corr[[1]]$pValues)
hist(ZI_Corr[[1]]$Correlation)
hist(Norm_Corr[[1]]$Correlation)
hist(ZI_Corr[[2]]$Correlation)
hist(Norm_Corr[[2]]$Correlation)
dev.off()


par(mfrow=c(3,2))
pdf('Coaccc_Examples.pdf')
p1 <- ggplot(joinedBack[[1]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround$Correlation[1], col='blue') + ggtitle('ZI Spearman')
p1.2 <- ggplot(joinedBack_N[[1]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround_Norm$Correlation[1], col='blue') +
ggtitle('Spearman')
x=200
p2.1 <- ggplot(joinedBack[[x]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround$Correlation[x], col='blue')
p2.2 <- ggplot(joinedBack_N[[x]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround_Norm$Correlation[x], col='blue')
x=700
p3.1 <- ggplot(joinedBack[[x]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround$Correlation[x], col='blue')
p3.2 <-ggplot(joinedBack_N[[x]], aes(x = Correlation)) + geom_histogram() +
     geom_vline(xintercept = foreGround_Norm$Correlation[x], col='blue')
gridExtra::grid.arrange(p1,p1.2, p2.1,p2.2, p3.1, p3.2, ncol = 2, nrow = 3)
dev.off()


pdf('CoAccessibility_Figure_Plots.pdf')
ggplot(ZI_Spearman, aes(x = Correlation, Group = Type, fill = Type)) + geom_density(alpha = 0.5) +  ggtitle('ZI Spearman: Distance = 0.34')  + theme_bw()
ggplot(Spearman, aes(x = Correlation, Group = Type, fill = Type)) + 
geom_density(alpha = 0.5) +  ggtitle('Spearman: Distance = 0.19') + theme_bw()
ggplot(data.frame(ZI = foreGround$FDR, Norm = foreGround_Norm$FDR)[sample(length(foreGround$FDR),2500),], 
       aes(x = Norm, y = ZI)) +
    geom_point() + theme_bw() + xlab('Normal Spearman FDR') + 
    ylab('Zero-Inflated Spearman FDR ') + 
    xlim(0,1) + ylim(0,1) + 
    ggtitle('Subsampled Spearman vs ZI Spearman FDRs') + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red')

ggplot(data.frame(ZI = foreGround$FDR, Norm = foreGround_Norm$FDR), 
       aes(x = Norm, y = ZI)) +
    stat_density_2d_filled() + theme_bw() + 
    xlab('Normal Spearman FDR') + 
    ylab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red')
dev.off()


df_FDRs <- data.frame(ZI = foreGround$FDR, Norm = foreGround_Norm$FDR) %>%
                mutate(Group = dplyr::case_when(ZI < 0.05 & Norm < 0.05 ~ 'Both',
                                         ZI < 0.05 & Norm >= 0.05  ~ 'ZI Spearman',
                                         ZI >= 0.05 & Norm < 0.05  ~ 'Spearman',
                                          ZI >= 0.05 & Norm >= 0.05  ~ 'Neither'))
correlationDF <- data.frame(Tile1 = foreGround$Tile1, Tile2 = foreGround$Tile2, 
                            ZI_Correlation = foreGround$Correlation, 
                            ZI_pValues = foreGround$pValues,
                            ZI_FDR = foreGround$FDR,
                            Norm_Correlation = foreGround_Norm$Correlation, 
                            Norm_pValues = foreGround_Norm$pValues,
                            Norm_FDR = foreGround_Norm$FDR)
write.csv(correlationDF, 'Correlation_Dataframe.csv')

pdf('CoAccessibility_MoreFigurePlots.pdf')
library(scattermore)

ggplot(df_FDRs, 
       aes(y = Norm, x = ZI, color = Group)) +
    geom_scattermore() + theme_bw() + 
    ylab('Spearman FDR') + 
    xlab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red') +
    color_scale_manual(values=c("red", "blue", "green"))


ggplot(df_FDRs, 
       aes(y = Norm, x = ZI, color = Group)) +
    geom_scattermore() + theme_bw() + 
    ylab('Spearman FDR') + 
    xlab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red') +
    scale_color_manual(values=c("red", "black", "purple", "green"))

ggplot(data.frame(ZI = foreGround$Correlation, Norm = foreGround_Norm$Correlation), 
       aes(x = Norm, y = ZI)) +
    stat_density_2d_filled() + theme_bw() + 
    xlab('Normal Spearman FDR') + 
    ylab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') + 
 geom_hline(yintercept = 0.45, color = 'red') 
dev.off()

##### Plot specific examples
ZI_Only <- which(correlationDF$ZI_FDR < 0.05 & correlationDF$Norm_FDR > 0.05)
Spear_Only <- which(correlationDF$ZI_FDR > 0.05 & correlationDF$Norm_FDR < 0.05)
Both_Only <- which(correlationDF$ZI_FDR < 0.05 & correlationDF$Norm_FDR < 0.05)

idx = abs(correlationDF$Norm_Correlation - correlationDF$ZI_Correlation)>0.2 & correlationDF$ZI_FDR < 0.05 & correlationDF$Norm_FDR > 0.05

idx2 = abs(correlationDF$Norm_Correlation - correlationDF$ZI_Correlation)>0.4 & correlationDF$ZI_FDR > 0.05 & correlationDF$Norm_FDR < 0.05


pdf('CoAccessibility_SpecificExamples.pdf', width = 3, height = 3)

plot_coaccessibility(correlationDF, index = which(idx)[1], peak_mat = accMat, addLabel = T)
plot_coaccessibility(correlationDF, index = which(idx)[2], peak_mat = accMat, addLabel = T)
plot_coaccessibility(correlationDF, index = which(idx)[3], peak_mat = accMat, addLabel = T)

for(i in 1:sum(idx2)){
plot_coaccessibility(correlationDF, index = which(idx2)[i], peak_mat = accMat, addLabel = T)
}
#plot_coaccessibility(correlationDF, index = which(idx2)[2], peak_mat = accMat, addLabel = T)
#plot_coaccessibility(correlationDF, index = which(idx2)[3], peak_mat = accMat, addLabel = T)
dev.off()

### Subset foreground and background to the largest correlation per database interaction. 
mergedSet_ZI <- dplyr::full_join(fTiledLinks, 
                                 correlationDF, by = c('Tile1' = 'Tile1', 'Tile2' = 'Tile2')) %>%
                dplyr::distinct() %>% 
                dplyr::group_by(Partition) %>% dplyr::arrange(ZI_FDR) %>% dplyr::slice_head(n=1) 
                dplyr::mutate(Set = 'ZI_Spearman')
mergedSet_N <- dplyr::full_join(fTiledLinks, 
                                 correlationDF, by = c('Tile1' = 'Tile1', 'Tile2' = 'Tile2')) %>%
                dplyr::distinct() %>% 
                dplyr::group_by(Partition) %>% dplyr::arrange(Norm_FDR) %>% dplyr::slice_head(n=1)# %>%
                dplyr::mutate(Set = 'Spearman')



subCorrelations <- rbind(filter(, mergedSet_N) %>% 
                    tidyr::pivot_wider(id_cols = c('Tile1', 'Tile2', 'Partition'),
                                                    names_from = 'Set', 
                                                    values_from = c('Correlation', 'pValues', 'FDR'))
filt_Core <- filter(correlationDF, 
                    (paste(Tile1, Tile2) %in% paste(mergedSet_N$Tile1,mergedSet_N$Tile2)) |
                     (paste(Tile1, Tile2) %in% paste(mergedSet_ZI$Tile1,mergedSet_ZI$Tile2))) %>%
                dplyr::inner_join(.,fTiledLinks, by = c('Tile1' = 'Tile1', 'Tile2' = 'Tile2')) %>%
                            dplyr::distinct()  %>%
                mutate(Group = dplyr::case_when(ZI_FDR < 0.05 & Norm_FDR < 0.05 ~ 'Both',
                                         ZI_FDR < 0.05 & Norm_FDR >= 0.05  ~ 'ZI Spearman',
                                         ZI_FDR >= 0.05 & Norm_FDR < 0.05  ~ 'Spearman',
                                          ZI_FDR >= 0.05 & Norm_FDR >= 0.05  ~ 'Neither'))

pdf('CoAccessibility_MoreFigures_UniqueInteractions.pdf')
ggplot(filt_Core, 
       aes(y = Norm_FDR, x = ZI_FDR, color = Group)) +
    geom_scattermore() + theme_bw() + 
    ylab('Spearman FDR') + 
    xlab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red') +
    scale_color_manual(values=c("red", "black", "purple", "green"))

ggplot(filt_Core, 
       aes(y = Norm_FDR, x = ZI_FDR, color = Group)) +
    geom_scattermore() + theme_bw() + geom_density_2d_filled() +
    ylab('Spearman FDR') + 
    xlab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red') +
    scale_color_manual(values=c("red", "black", "purple", "green"))

ggplot(filt_Core, 
       aes(y = Norm_FDR, x = ZI_FDR, color = Group)) +
    stat_density_2d_filled() + theme_bw() + 
    ylab('Spearman FDR') + 
    xlab('Zero-Inflated Spearman FDR ') + xlim(0,1) + 
    ylim(0,1) + geom_vline(xintercept = 0.05, color = 'red') +
     geom_hline(yintercept = 0.05, color = 'red') +
    scale_color_manual(values=c("red", "black", "purple", "green"))

dev.off()

# ######################## functions


calculatePValues <- function(foreground, background, cl){
    
    pValues <- pbapply::pblapply(1:length(foreground$Correlation), function(x){

        cor1 <- foreground$Correlation[x]

        if(cor1 >= 0){

            sum(cor1 < background$Correlation, na.rm = TRUE)/length(background$Correlation)

        }else if(cor1 < 0){

            sum(cor1 > background$Correlation, na.rm = TRUE)/length(background$Correlation)

        }

    }, cl = cl) %>% unlist()
    
  return(pValues)

}

matchInteractions <- function(filteredCorrelations, enhancerLinks, returnTileLinks = FALSE){
    
    ## Let's tile the interactions from the database, so that they match the correlations. 
    Tile1Links <- MOCHA::StringsToGRanges(enhancerLinks$Tile1)
    Tile1Links <-  plyranges::stretch(plyranges::anchor_end(Tile1Links), 
                                   extend = start(Tile1Links) %% 500)
    Tile1Links <-  plyranges::stretch(plyranges::anchor_start(Tile1Links), 
                                   extend = 499 - (end(Tile1Links) %% 500))     
    Tile2Links <- MOCHA::StringsToGRanges(enhancerLinks$Tile2)
    Tile2Links <-  plyranges::stretch(plyranges::anchor_end(Tile2Links), 
                                       extend = start(Tile2Links) %% 500)
    Tile2Links <-  plyranges::stretch(plyranges::anchor_start(Tile2Links), 
                                       extend = 499 - (end(Tile2Links) %% 500))    

    Tile1Links2 <- plyranges::tile_ranges(Tile1Links, 500) 
    Tile2Links2 <- plyranges::tile_ranges(Tile2Links, 500)
    
    #Now let's put the dataframe back together of all tile-tile relationship found within our database. 
    tiledLinks <- dplyr::full_join(
        data.frame(Tile1 = MOCHA::GRangesToString(Tile1Links2), Partition = Tile1Links2$partition),
        data.frame(Tile2 = MOCHA::GRangesToString(Tile2Links2), Partition = Tile2Links2$partition),
        by = 'Partition')
    
    if(returnTileLinks){
        
        return(tiledLinks)
        
    }
    ## join this tiledLinks list with the correlation dataframe, so that we only have validated interactions with correlations.
    filterMatchedCorr <- rbind(dplyr::inner_join(tiledLinks,filteredCorrelations, 
                                           by = c('Tile1' = 'Tile1', 'Tile2' = 'Tile2')),
                               dplyr::inner_join(tiledLinks,filteredCorrelations, 
                                           by = c('Tile1' = 'Tile2', 'Tile2' = 'Tile1')))
    return(filterMatchedCorr)

}

