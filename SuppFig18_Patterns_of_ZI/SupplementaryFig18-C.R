###############################################################################
#install.packages(c('lme4','glmmTMB','ggbreak','lmerTest','DHARMa'))
require(data.table)
require(MOCHA)
require(ggplot2)
require(parallel)
require(ggpubr)
library(ArchR)
library(GenomicRanges)
library(plyranges)
require(parallel)
require(pbapply)
require(glmmTMB)
require(DHARMa)


# Load the ArchR Project
ArchRProject <- loadArchRProject("/home/jupyter/FullCovid")
source('/home/jupyter/theme.R')

#######################################################################
metadf_dt <- as.data.table(getCellColData(ArchRProject)) 
numCores=10
metaColumn='new_sample_cellType'

### Select 10 cell types
### To Show normalization 
celltypes = unique(metadf_dt$predictedGroup_Co2)[1:10]
blackList = getBlacklist(ArchRProject)

# ### subset early visit

lookup_table <- unique(metadf_dt[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]


cell_counts <-  metadf_dt[new_cellType=='CD16 Mono'
                          ,.N, by=list(new_cellType, Sample)]

lookup_table <- dplyr::left_join(lookup_table, cell_counts)

xsec_lookup = lookup_table[Sample %in% samplesToKeep]
xsec_lookup = xsec_lookup[Sample %in% colnames(stm)]

## Load the ArchR Project
## and extract the metadata object
ArchRProj <- loadArchRProject("/home/jupyter/FullCovid")
metadata = as.data.table(getCellColData(ArchRProj))
studySignal = median(metadata$nFrags)

# Define your annotation package for TxDb object(s)
# and genome-wide annotation
# Here our samples are human using hg38 as a reference.
# For more info: https://bioconductor.org/packages/3.15/data/annotation/

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

# ##########################################################
# ##########################################################

## Get metadata information
## at the sample level
lookup_table <- unique(metadata[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'Age','Sex',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table$Sample[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]
xsec_lookup = lookup_table[Sample %in% samplesToKeep]
## subset ArchR Project
idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
cellsSample <- ArchRProj$cellNames[idxSample]
ArchRProj <- ArchRProj[cellsSample, ]

# ##########################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
# ##########################################################


tileResults <- callOpenTiles( 
    ArchRProj,
    cellPopLabel = "CellSubsets" ,
    cellPopulations = c("CD16 Mono",'CD4 Naive','B naive'),
    TxDb = 'TxDb.Hsapiens.UCSC.hg38.refGene',
    Org = 'org.Hs.eg.db',
    numCores = 5,
    outDir = 'tmp',
    studySignal = studySignal
)


###########################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
###########################################################
groupColumn <- "COVID_status" 
threshold <- 0.2

SampleTileMatrices <- lapply(c("CD16 Mono",'CD4 Naive','B naive'),
                             function(x)
                                 MOCHA::getSampleTileMatrix(
                                      tileResults,
                                      cellPopulations = x, 
                                      groupColumn = groupColumn,
                                      threshold = threshold,
                                      verbose = FALSE
                                    )
                             )

### Extract CD16 Sample tile matrix 
### and log transfrom 

test_zero_inflation <- function(cellIndex, cellname='CD16'){
    
    stm = as.data.frame(assays(SampleTileMatrices[[cellIndex]])[[1]])
    stm[is.na(stm)]<- 0
    #cd16_stm = log2(cd16_stm +1)
    stm$Tile = row.names(stm)
    stm = as.data.table(stm)

    ### Confirm Metadata and column
    ### names match to measure zero-inflation 
    xsec_lookup$Sample == colnames(stm)[1:39]

    ### create function to test zero inflation 
    test_zi <- function(accessibility, Tile){

        xsec_lookup$Accessibility <- round(as.numeric(accessibility))

        fittedModelNB1 <- glmmTMB(Accessibility ~  COVID_status + Age+ Sex, 
                               family=nbinom1,
                               data = xsec_lookup)
        fittedModelNB2 <- glmmTMB(Accessibility ~  COVID_status + Age+ Sex, 
                               family=nbinom2,
                               data = xsec_lookup)

        fittedModelPoisson <- glmmTMB(Accessibility ~ COVID_status + Age+ Sex , 
                               family=poisson,
                               data = xsec_lookup)

        test_NBModel1 <- DHARMa::testZeroInflation(fittedModelNB1)
        test_NBModel2 <- DHARMa::testZeroInflation(fittedModelNB2)

        res = data.frame(
            Tile = Tile,
            NB1_Pvalue = test_NBModel1$p.value,
            NB2_Pvalue = test_NBModel2$p.value
        )

        return(res)
    }

    cl <- makeCluster(60)
    clusterExport(cl, c("test_zi",'nbinom1','stm','nbinom2',
                        "glmmTMB","xsec_lookup","glmmTMBControl"),
                          envir=environment())

    start = Sys.time()
    zero_inflation_testing = parLapply(1:nrow(stm),
                 function(x)
                    try(test_zi(accessibility=stm[x,1:39], Tile=stm$Tile[x])),
                          cl=cl)

    Sys.time()-start

    stopCluster(cl)

    idx <- which(sapply(zero_inflation_testing, class) !='try-error')
    final_res <- rbindlist(zero_inflation_testing[idx])


    pdf(paste(cellname, 'ZI_Pval_NB1.pdf',sep=''))
    p<- ggplot(final_res,
           aes(x=NB1_Pvalue))+geom_histogram()+
            ggtitle('Negative Binom1 Distribution')
    print(p)
    dev.off()


    pdf(paste(cellname, 'ZI_Pval_NB2.pdf',sep=''))
    
    p = ggplot(final_res,
           aes(x=NB2_Pvalue))+geom_histogram()+
            ggtitle('Negative Binom2 Distribution')
    print(p)
    dev.off()
    
    return(final_res)
    
}


res1= test_zero_inflation(2, cellname='CD4 Naive')
res2= test_zero_inflation(3, cellname='B Naive')
res3= test_zero_inflation(1, cellname='CD16 Mono')

write.csv(res1, file='cd4 naive res.csv')
write.csv(res2, file='b naive res.csv')
write.csv(res3, file='cd16 mono res.csv')

### calculate statistics for 
### for quantifying number of 
### tiles with p < 0.05 
res1 = fread('cd4 naive res.csv')
res2 = fread('b naive res.csv')
res3 = fread('cd16 mono res.csv')

## CD4 Naive T Cell
round(c(mean(res1$NB1_Pvalue < 0.05),
mean(res1$NB2_Pvalue < 0.05)),2)

## Naive B cells
round(c(mean(res2$NB1_Pvalue < 0.05),
mean(res2$NB2_Pvalue < 0.05)),2)

## CD16 monocytes
round(c(mean(res3$NB1_Pvalue < 0.05),
  mean(res3$NB2_Pvalue < 0.05)),2)

### calculate statistics for 
### for quantifying number of 
### tiles with p < 0.05 
res1 = fread('cd4 naive res.csv')
res2 = fread('b naive res.csv')
res3 = fread('cd16 mono res.csv')

## CD4 Naive T Cell
round(c(mean(res1$NB1_Pvalue < 0.05 & res1$NB1_Statistic >=1),
mean(res1$NB2_Pvalue < 0.05 & res1$NB2_Statistic >=1)),2)

## Naive B cells
round(c(mean(res2$NB1_Pvalue < 0.05 & res2$NB1_Statistic >=1),
mean(res2$NB2_Pvalue < 0.05 & res2$NB2_Statistic >=1)),2)

## CD16 monocytes
round(c(mean(res3$NB1_Pvalue < 0.05 & res3$NB1_Statistic >=1),
mean(res3$NB2_Pvalue < 0.05 & res3$NB2_Statistic >=1)),2)

cd16 <- read.csv('cd16 mono res.csv')[,-1]
cd4 <- read.csv('cd4 naive res.csv')[,-1]
bnaive <- read.csv('b naive res.csv')[,-1]

pdf('SupplementalFigure18_ZeroInflation.pdf')

lapply(list(cd16, cd4, bnaive), function(XX){
    
    p1 <- ggplot(XX,
           aes(x=NB1_Pvalue))+geom_histogram() +
            ggtitle('Negative Binom1 Distribution')+
                geom_vline(xintercept = 0.05, color = 'red') + theme_bw() +
                ylab('Counts') + xlab('P-value')
    
    p2 <- ggplot(XX,
               aes(x=NB2_Pvalue))+geom_histogram()+
                ggtitle('Negative Binom2 Distribution') +
                geom_vline(xintercept = 0.05, color = 'red') + theme_bw()+
                ylab('Counts') + xlab('P-value')
    
    nb1 <- dplyr::mutate(XX, Group = 
        ifelse(NB1_Statistic > 1, 'Overinflated', 'Underinflated')) %>%
        dplyr::group_by(Group) %>%
        dplyr::filter(NB1_Pvalue < 0.05) %>%
        dplyr::summarize(TileNum = dplyr::n())
    nb2 <- dplyr::mutate(XX, Group = 
        ifelse(NB2_Statistic > 1, 'Overinflated', 'Underinflated')) %>%
        dplyr::group_by(Group) %>%
        dplyr::filter(NB2_Pvalue < 0.05) %>%
        dplyr::summarize(TileNum = dplyr::n())
    
    if(any(!grepl('Under', nb1$Group))){
        nb1 <- rbind(nb1, data.frame(Group = 'Underinflated', TileNum = 0))
    }
    if(any(!grepl('Under', nb2$Group))){
        nb2 <- rbind(nb2, data.frame(Group = 'Underinflated', TileNum = 0))
    }
    nb2$Group = factor(nb2$Group, levels = c('Underinflated', 'Overinflated'))
    nb1$Group = factor(nb1$Group, levels = c('Underinflated', 'Overinflated'))
    
    p3 <- ggplot(nb1, aes(x = Group, y = TileNum)) +
            geom_col() + theme_bw() + ylab('Number of Tiles') + xlab(NULL) +
        ggtitle('Under- vs Over-inflation of Zeroes')
    
    p4 <- ggplot(nb2, aes(x = Group, y = TileNum)) +
            geom_col() + theme_bw() + ylab('Number of Tiles') + xlab(NULL) +
        ggtitle('Under- vs Over-inflation of Zeroes')
    
    list(p1, p2, p3, p4)
    
    })
  
dev.off()

lapply(list(cd16, cd4, bnaive), function(XX){
    
       nb1 <- dplyr::mutate(XX, Group = 
        ifelse(NB1_Statistic > 1, 'Overinflated', 'Underinflated')) %>%
        dplyr::group_by(Group) %>%
        dplyr::filter(NB1_Pvalue < 0.05) %>%
        dplyr::summarize(TileNum = dplyr::n())
    nb2 <- dplyr::mutate(XX, Group = 
        ifelse(NB2_Statistic > 1, 'Overinflated', 'Underinflated')) %>%
        dplyr::group_by(Group) %>%
        dplyr::filter(NB2_Pvalue < 0.05) %>%
        dplyr::summarize(TileNum = dplyr::n())
    list(nb1, nb2)
    
    })