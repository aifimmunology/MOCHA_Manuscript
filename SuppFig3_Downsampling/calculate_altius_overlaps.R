## altius pre-factor analysis
require(ArchR)
require(MOCHA)
require(plyranges)
require(parallel)
source('/home/jupyter/theme.R')

setwd('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling')
hg38_altius <- readRDS('hg38_altius_gr.rds')
hg19_altius <- readRDS('hg19_altius_gr.rds')
lp <- ArchR::loadArchRProject('/home/jupyter/longPilot')
bm <- ArchR::loadArchRProject('/home/jupyter/DoubletFreeBoneMarrow/')
cov <- ArchR::loadArchRProject('/home/jupyter/FullCovid')

### Get BM Altius peakset Overlap

cl = parallel::makeCluster(50)
parallel::clusterEvalQ(cl, {library(ArchR)})
parallel::clusterEvalQ(cl, {library(rhdf5)})

frags_bm <-  MOCHA:::getPopFrags(bm, 
                                 metaColumn = "predictedGroup", 
                         cellSubsets = NULL, 
                         region = NULL,
                         numCores = cl,
                         blackList = NULL,
                         verbose = TRUE,
                         overlapList = 50
      )

bm_altius <- mclapply(frags_bm,
       function(x)
    sum(countOverlaps(x, hg19_altius)>0) / length(x),
       mc.cores=10
                      )
rm(frags_bm)
gc()
pryr::mem_used()
bm_altius_summary <- unlist(bm_altius)

### Get Long Pilot Altius peakset Overlap
cl = parallel::makeCluster(20)
parallel::clusterEvalQ(cl, {library(ArchR)})
parallel::clusterEvalQ(cl, {library(rhdf5)})

frags_lp <- MOCHA:::getPopFrags(
          ArchRProj = lp,
          metaColumn = 'predictedCell_Col2.5',
          cellSubsets = NULL,
          region = NULL,
          numCores = cl,
          blackList = NULL,
          verbose = TRUE,
          overlapList = 50
      )
stopCluster(cl)

frags_lp <-  MOCHA::getPopFrags(lp, 
                                 metaColumn = "predictedCell_Col2.5", 
                         numCores= 30)


lp_altius <- mclapply(frags_lp,
       function(x)
    sum(countOverlaps(x, hg38_altius)>0) / length(x),
       mc.cores=10
                      )
rm(frags_lp)

lp_altius_summary <- unlist(lp_altius)
pryr::mem_used()


### Get Covid Altius peakset Overlap
### clear memory first
cl = parallel::makeCluster(20)
parallel::clusterEvalQ(cl, {library(ArchR)
                            library(rhdf5) })


frags_covid <-  MOCHA::getPopFrags(cov, 
                                 metaColumn = "predictedGroup_Co2", 
                                 cellSubsets = NULL, 
                                 region = NULL,
                                 numCores = cl,
                                 sampleSpecific = T,
                                 NormMethod = "nfrags",
                                 blackList = NULL,
                                 verbose = TRUE,
                                 overlapList = 50
      )

cov_altius <- mclapply(frags_covid,
       function(x)
    sum(countOverlaps(x, hg38_altius)>0) / length(x),
       mc.cores=10
                      )


covid_altius_summary <- unlist(cov_altius)

df = data.frame(FRIA = c(bm_altius_summary, 
                  lp_altius_summary,
                  covid_altius_summary),
                Dataset=c(rep('Bone Marrow', length(bm_altius_summary)),
                          rep('Healthy Donors', length(lp_altius_summary)),
                          rep('COVID-19', length(covid_altius_summary)))
                )

df$Dataset <- factor(df$Dataset, levels=c('COVID-19',
                                          'Healthy Donors',
                                          'Bone Marrow'))

pdf('/home/jupyter/MOCHA_Manuscript/Fig2/altius.pdf',
    width=8,
    height=6)
ggplot(df,
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
    Dataset= c(rep('COVID-19',nrow(metadf_covid)),
               rep('Bone Marrow',nrow(metadf_bm)),
               rep('Healthy Donors',nrow(metadf_lp)))
    )

df_frags$Dataset <- factor(df_frags$Dataset, levels=c('COVID-19',
                                          'Healthy Donors',
                                          'Bone Marrow'))

pdf('/home/jupyter/MOCHA_Manuscript/Fig2/nfrags.pdf',
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

pdf('/home/jupyter/MOCHA_Manuscript/Fig2/nfrags2.pdf',
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
              
                              