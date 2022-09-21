##################################################################
##################################################################

### This script generates Panel B results for the 
### MOCHA manuscript. 

##################################################################
##################################################################

rm(list=ls())

########################################################################
########################################################################

### load libraries 
require(data.table)
require(ggplot2)
require(ggpubr)

### load datasets
archr = fread('n-1.daps.ArchR.csv')
scMACS= fread('n-1.daps.MOCHA.csv')
seurat= fread('n-1.daps.Signac.csv')

########################################################################
########################################################################
### Global Comparisons
scMACS$Conserved_Recall[1] <- scMACS$Common[1]
archr$Conserved = archr$Recall * archr$PeakNumber[1]

### Extract Global
### Precision & Recall 

df1= data.table(
  Iteration=rep(c(1:40),3),
  ConservedPeaks= c(archr$Conserved, scMACS$Conserved_Recall, seurat$Conserved),
  NewPeaks= c(archr$Unique+archr$PeakNumber[1], 
              scMACS$CountUnique + scMACS$Common[1], 
              seurat$Unique + seurat$Conserved[1]),
  Model = c(rep('Archr',nrow(archr)), 
            rep('scMACS',nrow(scMACS)),
            rep('Seurat',nrow(seurat)))
)
   

### Relative table
df_relative = data.frame(
  Model=c('ArchR','ArchR','Seurat','Seurat','scMACS','scMACS'),
  Value=c(2.54, 0.2, .58, .36, 0.85, 0.6),
  Evaluation=rep(c('New','Conserved'),3)
)
df_relative$Model <- factor(df_relative$Model,
                               levels=c('ArchR','Seurat','scMACS'))

### Add +1 to make 
### size relative to 
### full peakset

df1$Model <- factor(df1$Model, levels=c('Archr','Seurat','scMACS'))
df2 = melt(df1, id.vars = c('Iteration','Model'))         
df2$ModelEvaluation = paste(df2$Model,df2$variable,sep=':')


pdf('panelC.pdf')

p <-  ggplot(df2,
       aes(x=Iteration,
           y=value,
           col=Model,
           linetype=variable
       ))+geom_line(size=2)+#facet_wrap(~variable, ncol=1, scales='free_y')+
  theme_minimal() + ylab('# DAPS') + xlab("Sample Removed")+
  theme(legend.position = 'none')+
  theme(legend.position = 'none', 
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        strip.text.x = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90),
        axis.text.y = element_text(size=14),
  )

print(p)
dev.off()
