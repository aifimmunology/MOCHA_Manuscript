# ##############################################################

# ################ Let's identify invariant CTCF sites

# ##############################################################


# #Download all CTCF sites - 
# http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/assembled/Oth.ALL.05.CTCF.AllCell.bed

# Download all blood CTCF sites
# http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/assembled/Oth.Bld.05.CTCF.AllCell.bed

# ### Plan: Find invariant CTCF sites by finding the genome wide coverage of CTCF peaks.
# ### Then identify regions that have the max value. 
# ### We will then use the number of fragments overlapping with invariant CTCF sites as #### the control variable. 


source('ChIPseq_ProcessingFunctions.R')
toRemove = c("%20", "%","<br>", "2B|B3","B3")

allCTCF  <- read_bed("Oth.Bld.05.CTCF.AllCell.bed") %>%
                processChipMeta(., "name", ";", toRemove)

length(table(allCTCF$name))
#825 CTCF experiments
length(table(allCTCF$Name))
## 205 unique tissues

MonoDC <- loadArchRProject("MonoDC") 
blackList <- getBlacklist(MonoDC)

covCTCF <- allCTCF %>% group_by(Name) %>%
                plyranges::reduce_ranges() %>%
                plyranges::compute_coverage() 
saveRDS(covCTCF, file = "CTCF_coverage.RDS")

covCTCF <- readRDS("CTCF_coverage.RDS")

# ## Remove overlaps with blacklisted regions
# ## Filter these windows for the altius peakset. Looks like not all CTCF ChipSeq peaks are accessible.
# ## Choose all CTCF ChIP-seq peaks that show up in 95% of cell lines/primary cell types
# ## Expand width to 250 bp window

altius <- readRDS("hg38_altius_gr.rds")
altiusStretch <- plyranges::stretch(anchor_center(altius), extend = 500-width(altius))

invarCTCF <- covCTCF %>% plyranges::filter_by_non_overlaps(blackList, minoverlap = 50) %>%
               plyranges::filter_by_overlaps(altiusStretch) %>%
               plyranges::filter(score >= 0.99*204)
fixedWindows <- plyranges::stretch(anchor_center(invarCTCF), extend = 500-width(invarCTCF))

write_bed(fixedWindows, file ="InvariableCTCF_Domains_Blood_500bp.bed")

exportCTCF <- covCTCF %>% plyranges::filter_by_non_overlaps(blackList, minoverlap = 50) %>%
               plyranges::filter(width(.) <= 40000) %>%
                plyranges::filter(score >= 0.95*204) %>%
                plyranges::reduce_ranges(maxScore = max(score),
                                         minScore = min(score),
                                         meanScore = mean(score))

write_bed(exportCTCF, file ="InvariableCTCF_Domains_Blood.bed")

# ##############################################################

# ################ ### Run through tests of variation.

# ##############################################################

FullCovid <- loadArchRProject('FullCovidf')

metadf <- getCellColData(FullCovid) %>% as.data.frame()

saveArchRProject(FullCovid)
blackList <- getBlacklist(FullCovid)

CellSubset <- grep("CD8 Naive|CD4 Naive|CD14 Mono|CD16 Mono|NKs", unique(metadf$CellSubsets), value = TRUE)

frags1 <- getPopFrags(FullCovid, "CellSubsets", 
                      cellSubsets = CellSubset,
                     numCores = 35, blackList = blackList)

#CTCF regions that show up in 95% of tissues (ChIP-seq data) & in Altius peakset
invarCTCF <- read_bed("InvariableCTCF_Domains_Blood_500bp.bed")

CTCFCounts <- lapply(frags1, function(x) {
    
        mclapply(x, function(y) { count_overlaps(invarCTCF, y)}, mc.cores = 30)

})
names(CTCFCounts) <- NULL
CTCFCounts <- unlist(CTCFCounts, recursive = FALSE)

#Make sure CTCF sites are all detected
any(unlist(lapply(CTCFCounts, function(x){all(sum(unlist(x)) ==0)})))

### Let's create the list of peaks, celltypes, and samples to match CTCFCounts, and invarCTCF\
peaks <- MOCHA::GRangesToString(invarCTCF)

SampleList <- names(CTCFCounts) %>% gsub(".*#","",.) %>% gsub("__.*", "",.) 

CellTypeList <- names(CTCFCounts) %>% gsub(".FS.*","", .) %>% gsub(".B0.*","", .)

########  Generate normalization factors
## Total Sample Fragments
allArrows <- getArrowFiles(FullCovid)
fragsList<-  mclapply(seq_along(allArrows), function(x){
                    getFragmentsFromArrow(allArrows[x])
            }, mc.cores = 10)
names(fragsList) <- allArrows

totalsampleFrags = lapply(as(fragsList, "GRangesList"), length)
    
    lapply(fragsList, length)

names(totalsampleFrags) = allArrows
saveRDS(totalsampleFrags, file ="TotalFragsments_perSample.RDS")
totalsampleFrags <- readRDS("TotalFragsments_perSample.RDS")
rm(fragsList)

## Cell Number per sample/cell type
cellNumber <- getCellColData(FullCovid, c('Sample', 'CellSubsets')) %>% as.data.frame() %>%
                group_by(Sample, CellSubsets) %>% summarize(CellNumber = dplyr::n())

cellNumberList <- unlist(lapply(seq_along(CTCFCounts) , function(x){
    
    normF <- filter(cellNumber, Sample == SampleList[x],
                    CellSubsets == gsub("_", " ", CellTypeList[x])) %>%
                    ungroup() %>%
                    select(CellNumber) %>% unlist()
    normF/1000
    }))
names(cellNumberList) <- names(CTCFCounts)

## Fragment number per sample/cell type
fragNumber <- lapply(frags1, lengths)
names(fragNumber ) <- NULL
fragNumber  <- unlist(fragNumber , recursive = FALSE)

#Normalize by the number of cells per 1000 in the group
nCellsCount <- lapply(seq_along(CTCFCounts), function(x) {
        CTCFCounts[[x]]/cellNumberList[[x]]
})
names(nCellsCount) <- names(CTCFCounts)

#Normalized by total nFrags by population
nFragsCount <-lapply(seq_along(CTCFCounts), function(x) {
        normF <- fragNumber[[x]]/10^6
        CTCFCounts[[x]]/normF
})
names(nFragsCount) <- names(CTCFCounts)

#normalize by total nFrags per sample 
sampleTotalCount <- lapply(seq_along(CTCFCounts), function(x) {
        normF <- as.numeric(totalsampleFrags[grepl(pattern = SampleList[x], names(totalsampleFrags))])/10^6
        CTCFCounts[[x]]/normF
})
names(sampleTotalCount) <- names(CTCFCounts)

    

### Function for creating matrix
createMatrix <- function(CountList, peaks){
    #### SAMIR code
    tmp <- as.data.table(CountList)
    tmp$Peaks <- peaks                  
    melt_tmp <- melt(tmp, id.var='Peaks')    
    melt_tmp$variable <- as.character(melt_tmp$variable)
    melt_tmp$cellType <-gsub('#.*','',melt_tmp$variable)
    melt_tmp$Sample <- gsub("__.*","", gsub('.*#','',melt_tmp$variable))
    return(melt_tmp)
}


allNorms <- list(CTCFCounts, nCellsCount, nFragsCount, sampleTotalCount)
names(allNorms) <-  c('Raw', 'NCells', 'NFrags', 'SampleFrags') 

allPeaks <- mclapply(seq_along(allNorms), function(x) { 
                peakMat <- createMatrix(allNorms[[x]], peaks)
                peakMat$NormMethod = rep(names(allNorms)[x], dim(peakMat)[1])
                peakMat
}, mc.cores = 45)
names(allPeaks) <- names(allNorms)

peakMat <- rbindlist(allPeaks)
saveRDS(peakMat, file = "PeakMatrix.RDS")

peakMat <- readRDS("PeakMatrix.RDS")

#Calculate CV and mean counts
CVdt <- peakMat[,list(CV = sd(value)/mean(value), MeanCount = mean(value)), by=list(Peaks, cellType, NormMethod)]

library(ggridges)
library(ggforce)

CVdt$NormMethod[grepl('NFrags', CVdt$NormMethod)] = 
                                        'Number of Fragments (MOCHA)'   
CVdt$NormMethod[grepl('NCells', CVdt$NormMethod)] = 
                                        'Cell Number'
CVdt$NormMethod[grepl('SampleFrags', CVdt$NormMethod)] = 
                                       'Total Fragment #'                                          
colnames(CVdt) <- c('CTCF_Site', 'Cell_Type', 'Normalization_Method', 'CV', 'MeanCount')
write.csv(CVdt, 'Normalization_CV_CTCF_Counts.csv')                                   

# #Plots for Figure 1


pdf('Fig1_CVPlots.pdf')
                               
    ggplot(CVdt, aes(x = CV, fill = Normalization_Method)) + geom_density_line(alpha = 0.5) +
                     theme_bw() +
                     xlim(0,2) + facet_wrap(~Cell_Type, ncol = 1, scales = "free_y") +
                     ggtitle("CV per Peak Across All Samples - Cropped") + 
                      xlab("Coefficient of Variation") +
                    ylab("Distribution")

dev.off()

# ################ Code for calculating CV across cell types by Sample


#Calculate CV and mean counts
CVdt3 <- peakMat[,list(CV = sd(value)/mean(value), MeanCount = mean(value)), by=list(Peaks, NormMethod)]


library(ggridges)
library(ggforce)

CVdt3$NormMethod[grepl('NFrags', CVdt3$NormMethod)] = 
                                        'MOCHA'   
CVdt3$NormMethod[grepl('NCells', CVdt3$NormMethod)] = 
                                        'Number of Cells'
CVdt3$NormMethod[grepl('SampleFrags', CVdt3$NormMethod)] =  'Total Fragments'            
CVdt3$NormMethod[grepl('Raw', CVdt3$NormMethod)] = 'Raw'                              
colnames(CVdt3) <- c('CTCF_Site', 'Normalization_Method', 'CV', 'MeanCount')
write.csv(CVdt3, 'Normalization_CVAcrossCellTypes_CTCF_Counts.csv')                                   

pdf("Fig1_CVPlots_Part2.pdf") 
                               
    ggplot(CVdt3, aes(x = CV, fill = Normalization_Method)) + geom_density_line(alpha = 0.5) +
                     theme_bw() +
                     xlim(0,2) + facet_wrap(~Normalization_Method, ncol = 1, scales = "free_y") +
                     ggtitle("CV per Peak Across All Samples - Cropped") + 
                      xlab("Coefficient of Variation") +
                    ylab("Distribution")
                               
    ggplot(CVdt3, aes(x = CV, fill =Normalization_Method)) + geom_density_line(alpha = 0.5) +
                     theme_bw() +
                     xlim(0,2) + 
                     ggtitle("CV per Peak Across All Samples - Cropped") + 
                      xlab("Coefficient of Variation") +
                    ylab("Distribution")

dev.off()


# #### 

