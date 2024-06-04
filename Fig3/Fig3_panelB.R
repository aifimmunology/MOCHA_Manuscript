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
MOCHA= fread('n-1.daps.MOCHA.csv')
signac= fread('n-1.daps.Signac.csv')

########################################################################
########################################################################
### Global Comparisons
MOCHA$Conserved_Recall[1] <- MOCHA$Common[1]
archr$Conserved = archr$Recall * archr$PeakNumber[1]

### Extract Global
### Precision & Recall 

df1= data.table(
  Iteration=rep(c(1:40),3),
  ConservedPeaks= c(archr$Conserved, MOCHA$Conserved_Recall, signac$Conserved),
  NewPeaks= c(archr$Unique+archr$PeakNumber[1], 
              MOCHA$CountUnique + MOCHA$Common[1], 
              signac$Unique + signac$Conserved[1]),
  Model = c(rep('Archr',nrow(archr)), 
            rep('MOCHA',nrow(MOCHA)),
            rep('Signac',nrow(signac)))
)
   

### Relative table
df_relative = data.frame(
  Model=c('ArchR','ArchR','Signac','Signac','MOCHA','MOCHA'),
  Value=c(2.54, 0.2, .58, .36, 0.85, 0.6),
  Evaluation=rep(c('New','Conserved'),3)
)
df_relative$Model <- factor(df_relative$Model,
                               levels=c('ArchR','Signac','MOCHA'))

### Add +1 to make 
### size relative to 
### full peakset

df1$Model <- factor(df1$Model, levels=c('Archr','Signac','MOCHA'))
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
