## altius pre-factor analysis
require(ArchR)
require(MOCHA)
require(plyranges)
require(parallel)
require(data.table)
source('/home/jupyter/theme.R')

setwd('/home/jupyter/MOCHA_Manuscript2/SuppFig3_Downsampling')
hg38_altius <- readRDS('hg38_altius_gr.rds')
hg19_altius <- readRDS('hg19_altius_gr.rds')
lp <- ArchR::loadArchRProject('/home/jupyter/longPilot')
bm <- ArchR::loadArchRProject('/home/jupyter/DoubletFreeBoneMarrow/')
cov <- ArchR::loadArchRProject('/home/jupyter/FullCovid')



### Get Covid Altius peakset Overlap
### clear memory first
cl = parallel::makeCluster(50)
parallel::clusterEvalQ(cl, {library(ArchR)})
parallel::clusterEvalQ(cl, {library(rhdf5)})
cov_meta = as.data.table(cov@cellColData)
ctypes_names = unique(cov_meta$new_cellType)


### Run cell type by cell type 
### to avoid memory issues with 
### loading all 1.4 million cells 
cov_altius <- lapply(ctypes_names,
       function(x){
           
           tmp_frags =  MOCHA::getPopFrags(cov, 
                                 cellPopLabel = "predictedGroup_Co2", 
                                 cellSubsets = x, 
                                 numCores = cl,
                                 poolSamples = TRUE,  
                                 verbose = TRUE)

            sum(countOverlaps(tmp_frags[[1]], hg38_altius)>0) / length(tmp_frags[[1]])}
                      )
stopCluster(cl)
covid_altius_summary <- unlist(cov_altius)
covid_mat = data.frame(
    CellType = ctypes_names,
    Altius =covid_altius_summary,
    Dataset ='COVID19'
)


### Get BM Altius peakset Overlap
cl = parallel::makeCluster(50)
parallel::clusterEvalQ(cl, {library(ArchR)})
parallel::clusterEvalQ(cl, {library(rhdf5)})
bm_meta = as.data.table(bm@cellColData)
ctypes_names = unique(bm_meta$new_cellType)

frags_bm <-  MOCHA:::getPopFrags(bm, 
                                 cellPopLabel = "predictedGroup", 
                                 cellSubsets = NULL, 
                                 numCores = cl,
                                 poolSamples = TRUE,
                                 verbose = TRUE)
 

bm_altius <- mclapply(ctypes_names,
       function(x){
           idx = grep(x, names(frags_bm))
           tmpFrags = IRanges::stack(frags_bm[idx])
           sum(countOverlaps(tmpFrags, hg19_altius)>0) / length(tmpFrags)
    },
                      mc.cores=10
                      )
rm(frags_bm)
gc()
pryr::mem_used()
bm_altius_summary <- unlist(bm_altius)

bm_mat = data.frame(
    CellType = ctypes_names,
    Altius =bm_altius_summary,
    Dataset ='Hematopoiesis'
)

### Get Long Pilot Altius peakset Overlap
cl = parallel::makeCluster(50)
parallel::clusterEvalQ(cl, {library(ArchR)})
parallel::clusterEvalQ(cl, {library(rhdf5)})


frags_lp <- MOCHA:::getPopFrags( ArchRProj = lp,
                                 cellPopLabel='predictedGroup_Col2.5',
                                 cellSubsets = NULL, 
                                 numCores = cl,
                                 poolSamples = TRUE,
                                 verbose = TRUE)
stopCluster(cl)

lp_meta = as.data.table(lp@cellColData)
ctypes_names = unique(lp_meta$predictedGroup_Col2.5)

lp_altius <- mclapply(ctypes_names,
       function(x){
           idx = grep(x, names(frags_lp))
           tmpFrags = IRanges::stack(frags_lp[idx])
           sum(countOverlaps(tmpFrags, hg38_altius)>0) / length(tmpFrags)
    },
                      mc.cores=10
                      )

rm(frags_lp)

lp_altius_summary <- unlist(lp_altius)
pryr::mem_used()


lp_mat = data.frame(
    CellType = ctypes_names,
    Altius =lp_altius_summary,
    Dataset ='HealthyDonors'
)


combined_df = rbind(
    covid_mat,
    bm_mat,
    lp_mat)

write.csv(combined_df,
          file='/home/jupyter/MOCHA_Manuscript2/SuppFig3_Downsampling/FRIA.csv'
          )

pdf('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling/altius.pdf',
    width=8,
    height=6)
ggplot(combined_df,
       aes(x=Dataset,
           y=FRIA))+
        geom_point()+geom_violin(scale='width')+ThemeMain+
        ggtitle('Fraction Reads in Altius Peakset by Cell Type')+ylim(c(0.3,1))+
        theme_minimal()+
        theme(text=element_text(size=20))+
        ggpubr::stat_compare_means()
dev.off()
                              

metadf_covid <-  as.data.table(getCellColData(cov))
metadf_bm <-  as.data.table(getCellColData(bm))
metadf_lp <-  as.data.table(getCellColData(lp))

df_frags <- data.table(
    nFrags = c(metadf_covid$nFrags,
               metadf_bm$nFrags,
               metadf_lp$nFrags),
    Sample = c(metadf_covid$Sample,
               metadf_bm$Sample,
               metadf_lp$Sample),
    Dataset= c(rep('COVID-19',nrow(metadf_covid)),
               rep('Hematopoiesis',nrow(metadf_bm)),
               rep('Healthy Donors',nrow(metadf_lp)))
    )

df_frags$Dataset <- factor(df_frags$Dataset, levels=c('COVID-19',
                                          'Healthy Donors',
                                          'Hematopoiesis'))

pdf('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling/nfrags.pdf',
    width=8,
    height=6)
ggplot(df_frags,
       aes(x=Dataset,
           y=nFrags))+
        geom_boxplot(scale='width')+ThemeMain+
        ggtitle('# Frags per Cell Across Datasets')+
        theme_minimal()+
        theme(text=element_text(size=20))+
ggpubr::stat_compare_means()
dev.off()

pdf('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling/nfrags2.pdf',
    width=8,
    height=6)
ggplot(df_frags,
       aes(x=Dataset,
           y=nFrags))+
        geom_violin(scale='width')+ThemeMain+
        ggtitle('# Frags per Cell Across Datasets')+
        theme_minimal()+
        theme(text=element_text(size=20))+
    ggpubr::stat_compare_means()
dev.off()
              
            
write.csv(df_frags,
          file='frags_by_dataset.csv')





metadf_covid <-  as.data.table(getCellColData(cov))
metadf_bm <-  as.data.table(getCellColData(bm))
metadf_lp <-  as.data.table(getCellColData(lp))


covid_cellCounts = metadf_covid[, list(Count=.N), 
                                by=list(Sample, CellSubsets)]

lp_cellCounts = metadf_lp[, list(Count=.N), 
                                by=list(Sample, predictedGroup_Col2.5)]

bm_cellCounts = metadf_bm[, list(Count=.N),
                          by=new_cellType]
bm_cellCounts$Sample = NA


covid_cellCounts$Dataset='COVID-19'
colnames(covid_cellCounts) = c('Sample','Cell','Count','Dataset')
lp_cellCounts$Dataset='HealthyDonors'
colnames(lp_cellCounts) = c('Sample','Cell','Count','Dataset')
bm_cellCounts$Dataset='Hematopoiesis'
colnames(bm_cellCounts) = c('Cell','Count','Sample','Dataset')


df = rbind(covid_cellCounts,
      lp_cellCounts,
      bm_cellCounts)

write.csv(df, 
          file='cell_counts.csv')