# ##################################################################################
# ##################################################################################

require(ggplot2)
require(dplyr)
require(data.table)



cellCounts = data.table(readxl::read_xlsx('MOCHA_Manuscript2/Fig2/SourceData_Figure2-1.xlsx'))


pdf('MOCHA_Manuscript2/Fig2/cell_counts.pdf',
    width=20, height=4)
p1= ggplot(cellCounts[Dataset=='COVID19X'],
       aes(x=CellSubsets,
           y=Count,
           fill=CellGroup))+
    geom_jitter(alpha = 0.5, width = 0.2)+ geom_violin(alpha=0.5)+
    theme_minimal()+ 
    scale_y_continuous(trans='log10', limits=c(1,5000))+
    theme(axis.text.x=element_text(size=10, angle=90),
          text = element_text(size=10),
         legend.position='none')  

p2= ggplot(cellCounts[Dataset=='HealthyDonor'],
       aes(x=CellSubsets,
           y=Count,
           fill=CellGroup))+
    geom_jitter(alpha = 0.5, width = 0.2)+ geom_violin(alpha=0.5)+
    theme_minimal()+
    scale_y_continuous(trans='log10', limits=c(1,5000))+
    theme(axis.text.x=element_text(size=10, angle=90),
          text = element_text(size=10),
         legend.position='none')

p3= ggplot(cellCounts[Dataset=='Hematopoiesis'],
       aes(x=reorder(CellSubsets, Count, mean),
           y=Count,
           fill=CellGroup))+
    geom_bar(stat='identity', position='dodge')+
    theme_minimal()+
    scale_y_continuous(trans='log10', limits=c(1,30000))+
    theme(axis.text.x=element_text(size=10, angle=90),
          text = element_text(size=10),
         legend.position='none')

ggpubr::ggarrange(p1,p2,p3, ncol=3)

dev.off()





tileCounts = data.table(readxl::read_xlsx('MOCHA_Manuscript2/Fig2/SourceData_Figure2-1.xlsx', sheet=2))

pdf('MOCHA_Manuscript2/Fig2/tile_counts.pdf',
    width=24, height=4)

p1= ggplot(tileCounts[Dataset=='COVID19X'],
       aes(x=Model,
           y=Tiles,
           fill=Model))+
    geom_jitter(alpha = 0.5, width = 0.2)+ geom_violin(alpha=0.5)+

    theme_minimal()+ facet_wrap(~CellPop)+
    scale_y_continuous(limits=c(1,1300000))+
    theme(axis.text.x=element_text(size=10),
          text = element_text(size=10),
         legend.position='none')  

p2= ggplot(tileCounts[Dataset=='HealthyDonor'],
       aes(x=Model,
           y=Tiles,
           fill=Model))+
   geom_jitter(alpha = 0.5, width = 0.2)+ geom_violin(alpha=0.5)+
    theme_minimal()+ facet_wrap(~CellPop)+
    scale_y_continuous(limits=c(1,1300000))+
    theme(axis.text.x=element_text(size=10),
          text = element_text(size=10),
         legend.position='none')  



p3= ggplot(tileCounts[Dataset=='Hematopoiesis'],
       aes(x=reorder(CellPop, Tiles, max), 
           y=Tiles,
           fill=Model))+
    geom_bar(stat='identity', position='dodge')+
    theme_minimal()+
    scale_y_continuous(limits=c(0,1300000))+
    theme(axis.text.x=element_text(size=10),
          text = element_text(size=10),
         legend.position='none')

ggpubr::ggarrange(p1,p2,p3, ncol=3)

dev.off()

### Get Exact P-values 
tiles = tileCounts[Dataset !='Hematopoiesis']
tiles$variable = tiles$Model

### HYPOTHESIS TEST
pairwise_wilcox <- function(groupA, groupB, cell,dataset){

    valsA = tiles[variable==groupA & CellPop==cell & Dataset==dataset]$Tiles
    valsB = tiles[variable==groupB & CellPop==cell & Dataset==dataset]$Tiles
    
    res = data.frame(
        A = groupA,
        medianA= median(valsA),
        B = groupB,
        medianB=median(valsB),
        P_value=wilcox.test(valsA, valsB)$p.value,
        Dataset = dataset,
        Cell = cell
    )
    res
}

tile_pvals <- rbind(
    pairwise_wilcox('MOCHA','MACS2','B naive','COVID19X'),
    pairwise_wilcox('MOCHA','MACS2','CD16 Mono','COVID19X'),
    pairwise_wilcox('MOCHA','MACS2','CD4 CTL TEM','COVID19X'),
    pairwise_wilcox('MOCHA','HOMER','B naive','COVID19X'),
    pairwise_wilcox('MOCHA','HOMER','CD16 Mono','COVID19X'),
    pairwise_wilcox('MOCHA','HOMER','CD4 CTL TEM','COVID19X'),
    pairwise_wilcox('MOCHA','MACS2','B naive','HealthyDonor'),
    pairwise_wilcox('MOCHA','MACS2','CD16 Mono','HealthyDonor'),
    pairwise_wilcox('MOCHA','MACS2','CD8 TEM','HealthyDonor'),
    pairwise_wilcox('MOCHA','HOMER','B naive','HealthyDonor'),
    pairwise_wilcox('MOCHA','HOMER','CD16 Mono','HealthyDonor'),
    pairwise_wilcox('MOCHA','HOMER','CD8 TEM','HealthyDonor')
)

write.csv(tile_pvals, file='MOCHA_Manuscript2/Fig2/Fig2B_pvals.csv')

summarized_tiles = tiles[, median(value), by=list(variable, CellPop)]


intensity_covid = fread('MOCHA_Manuscript2/Fig2/Covid-19/results/intensity_covid.csv')
intensity_lp= fread('MOCHA_Manuscript2/Fig2/LongPilot/results/intensity_LP.csv')
intensity_bm= fread('MOCHA_Manuscript2/Fig2/BoneMarrow/results/intensity_bm.csv')

intensities_df = rbind(intensity_covid,
                       intensity_lp,
                       intensity_bm)


### HYPOTHESIS TEST
pairwise_wilcox_intensity <- function(cell,dataset){

    valsA = intensities_df[tile=='MOCHA Unique' & Cell==cell & Dataset==Dataset]$values
    valsB = intensities_df[tile!='MOCHA Unique' & Cell==cell & Dataset==Dataset]$values
    
    res = data.frame(
        A = 'MOCHA Unique',
        medianA= median(valsA),
        B = 'Missed by MOCHA',
        medianB=median(valsB),
        P_value=wilcox.test(valsA, valsB)$p.value,
        Dataset = Dataset,
        Cell = cell
    )
    res
}

intensities_pvals <- rbind(
    pairwise_wilcox_intensity('B naive','COVID19X'),
    pairwise_wilcox_intensity('CD16 Mono','COVID19X'),
    pairwise_wilcox_intensity('CD4 CTL TEM','COVID19X'),
    
    pairwise_wilcox_intensity('B naive','HealthyDonor'),
    pairwise_wilcox_intensity('CD16 Mono','HealthyDonor'),
    pairwise_wilcox_intensity('CD8 TEM','HealthyDonor'),
    
    pairwise_wilcox_intensity('10_cDC','BoneMarrow'),
    pairwise_wilcox_intensity('20_CD4.N1','BoneMarrow'),
    pairwise_wilcox_intensity('12_CD14.Mono.2','BoneMarrow'))

write.csv(intensities_pvals, file='MOCHA_Manuscript2/Fig2/Fig2C_pvals.csv')
