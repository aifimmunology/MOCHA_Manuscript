### Get cell counts all datasets 
require(MOCHA)
require(ArchR)
require(ggplot2)
require(data.table)

setwd('/home/jupyter/MOCHA_Manuscript/Fig2')

#### CovidDatasets
covid = loadArchRProject('/home/jupyter/FullCovid')
covid_meta <- as.data.table(getCellColData(covid))
covid_cell_counts <- covid_meta[, .N, by=c("Sample",'CellSubsets')]
cell_factors =covid_cell_counts[, median(N), by=CellSubsets]
cell_factors = cell_factors %>% arrange(-V1)
covid_cell_counts$CellSubsets= factor(covid_cell_counts$CellSubsets,
                                      levels=cell_factors$CellSubsets)

setwd('/home/jupyter/MOCHA_Manuscript/SuppFig3_Downsampling')

pdf('cell_counts_covid.pdf')
ggplot(covid_cell_counts,
       aes(x=CellSubsets,
           y=N))+
        geom_boxplot()+
   scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                 axis.text.y=element_text(size=16),
                 axis.text.x=element_text(size=16, angle=90),
                 title=element_blank())
dev.off()


#### Long Pilot
lp = loadArchRProject('/home/jupyter/longPilot')
lp_meta <- as.data.table(getCellColData(lp))
lp_cell_counts <- lp_meta[, .N, by=c("Sample",'predictedGroup_Col2.5')]
cell_factors =lp_cell_counts[, median(N), by=predictedGroup_Col2.5]
cell_factors = cell_factors %>% arrange(-V1)
lp_cell_counts$CellSubsets= factor(lp_cell_counts$predictedGroup_Col2.5,
                                      levels=cell_factors$predictedGroup_Col2.5)

pdf('cell_counts_lp.pdf')
ggplot(lp_cell_counts,
       aes(x=CellSubsets,
           y=N))+
        geom_boxplot()+
   scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                 axis.text.y=element_text(size=16),
                 axis.text.x=element_text(size=16, angle=90),
                 title=element_blank())
dev.off()


#### Bone Marrow
BM = loadArchRProject('/home/jupyter/DoubletFreeBoneMarrow')
BM_meta <- as.data.table(getCellColData(BM))
BM_cell_counts <- BM_meta[, .N, by=c("Sample",'new_cellType')]
cell_factors =BM_cell_counts[, sum(N), by=new_cellType]
cell_factors = cell_factors %>% arrange(-V1)
cell_factors$CellSubsets= factor(cell_factors$new_cellType,
                                      levels=cell_factors$new_cellType)


pdf('cell_counts_BM.pdf')
ggplot(cell_factors,
       aes(x=CellSubsets,
           y=V1))+
        geom_bar(stat='identity')+
   scale_y_continuous(trans='log10')+  theme_minimal()+
            theme(legend.position='none',
                 axis.text.y=element_text(size=16),
                 axis.text.x=element_text(size=16, angle=90),
                 title=element_blank())
dev.off()



### 
lp_cell_counts$Data = 'HealthyDonors'
lp_cell_counts$Cells= lp_cell_counts$CellSubsets
lp_cell_counts$Counts = lp_cell_counts$N


cell_factors$Data='BoneMarrow'
cell_factors$Cells=cell_factors$CellSubsets
cell_factors$Counts=cell_factors$V1

covid_cell_counts$Data ='COVID19'
covid_cell_counts$Cells =covid_cell_counts$CellSubsets
covid_cell_counts$Counts = covid_cell_counts$N


all_cell_counts = rbind(
    covid_cell_counts[,c('Data','Cells','Counts')],
    lp_cell_counts[,c('Data','Cells','Counts')],
    cell_factors[,c('Data','Cells','Counts')]   
    )


write.csv(all_cell_counts,
          file='cellCounts.csv')


### 
frags_df = data.table(
    Data = c(rep('COVID19', nrow(covid_meta)),
             rep('HealthyDonors', nrow(lp_meta)),
             rep('BoneMarrow', nrow(BM_meta))),    
    nFrags= c(covid_meta$nFrags,
              lp_meta$nFrags,
              BM_meta$nFrags)
    )

frags_df$Data = factor(frags_df$Data,
                       levels=c('COVID19','HealthyDonors','BoneMarrow'))

pdf('nFrags.pdf')
ggplot(frags_df,
       aes(x=Data,
           y=nFrags))+geom_violin()+
    theme_minimal()
dev.off()

write.csv(frags_df,
          file='frags.csv')

save(frags_df,
          file='frags_per_cell.RDS')