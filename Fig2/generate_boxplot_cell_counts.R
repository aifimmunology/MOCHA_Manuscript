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