#######################################################################################
#######################################################################################

### Generates graph
### for Figure 1 panel e
### Consensus reproducibility
#######################################################################################
#######################################################################################

## load libraries 
require(MOCHA)
require(ggpubr)
require(data.table)
require(ggplot2)
require(ArchR)
require(MultiAssayExperiment)
require(RaggedExperiment)

#######################################################################################
#######################################################################################
### set directory 
homeDir = '/home/jupyter/MOCHA_Manuscript/Fig2/'
setwd(homeDir)
source('../theme.R')
source('helper_granges.R')
source('utils.R')

## load ArchR project
ArchRProj = loadArchRProject('/home/jupyter/FullCovid')
metadata = as.data.table(ArchRProj@cellColData)
metadf_dt <- as.data.table(metadata) 
numCores=10

#######################################################################################
#######################################################################################

### 
celltypes = unique(metadf_dt$predictedGroup_Co2)
blackList = getBlacklist(ArchRProj)

#### subset early visit
lookup_table <- unique(metadf_dt[,c('Sample',
                                 'COVID_status',
                                 'Visit',
                                 'days_since_symptoms'),       
                              with=F])

## Subset to visit 1 and extract samples
samplesToKeep <- lookup_table[lookup_table$Visit =='FH3 COVID-19 Visit 1' & lookup_table$days_since_symptoms <= 15 | is.na(lookup_table$days_since_symptoms)]

#### and extract only 
#### 10 from each group 
covid_neg <- samplesToKeep[COVID_status =='Negative']$Sample
covid_pos <- samplesToKeep[COVID_status =='Positive']$Sample

#######################################################################################
#######################################################################################
mocha = readRDS("/home/jupyter/MOCHA_Manuscript/Fig2/Covid-19/MOCHA.RDS")
cd16 = mocha[['CD16 Mono']]

### identify positive & negative samples 
positive = cd16[,colnames(cd16) %in% covid_pos] 
negative = cd16[,colnames(cd16) %in% covid_neg] 

covidPeaks = compactAssay(positive, 'peak')
controlPeaks = compactAssay(negative, 'peak')

calculate_consensus <- function(covidPeaks){
        grid = c(0.0001,seq(0.05, 0.5, by=0.05))
        numPeaks = rowSums(covidPeaks, na.rm=T)/ncol(covidPeaks)
        
        consensus = sapply(grid, 
               function(x)
               sum(numPeaks >= x, na.rm=T
                  )
               )
    
        res = data.table(
            grid=grid,
            peak=consensus)
    
    return(list(Res=res,
                Tiles= names(numPeaks)[which(numPeaks >=0.2)]
               )
           )
}
           
           

ctl_repro = calculate_consensus(controlPeaks)
cov_repro = calculate_consensus(covidPeaks)

union = union(ctl_repro$Tiles, cov_repro$Tiles)
length(union)


tsam = compactAssay(cd16, 'TotalIntensity')
tsam = tsam[row.names(tsam) %in% union,]
mean(is.na(tsam))

df = data.frame(
    Grid = ctl_repro$Res$grid,
    'covid'= cov_repro$Res$peak,
    'control'=ctl_repro$Res$peak
    )
colnames(df) = c('Cutoff','COVID(+)','COVID(-)')
    
df = melt(df, id.var='Cutoff')
df$Cohort=df$variable
df = data.table(df)

setwd('../Fig1/')
pdf("Consensus_1e.pdf", width=8, height=4)
ggplot(df[Cutoff<=0.5], 
       aes(x=Cutoff,
           y=value,
           col=Cohort,
          fill=Cohort))+geom_point()+
    geom_vline(xintercept=0.2, lty=2, col='black')+
    geom_smooth()+theme_minimal()+

            theme(strip.text = element_text(size=25),
              axis.text.x = element_text(size=14, angle=90),
              axis.text.y =  element_text(size=14),
                 legend.position='right')
dev.off()

write.csv(df,
          file='../Fig1/supplementalFiles/Fig1e_consensus_plot.csv',
          row.names=F)