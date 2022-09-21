##################################################################
##################################################################

### This script generates Panel C results for the 
### MOCHA manuscript. 

##################################################################
##################################################################

require(data.table)
require(ggplot2)

signac = fread('Seurat_RunTimeAnalysis.csv')
archr = fread('ArchR_RunTime.csv')
scmacs = fread('MOCHA_runtime.csv')

### format 
signac = signac[, c('TotalPeaks','RunTime')]
signac$Model = 'Signac'

### 
archr$RunTime = archr$RunTime * 60
archr$TotalPeaks = archr$PeakNumber
archr$Model = 'ArchR'
archr = archr[, c('TotalPeaks','RunTime','Model')]

### 
scmacs$TotalPeaks = scmacs$PeaksTested
scmacs = scmacs[, c('TotalPeaks','RunTime','Model')]

### combined data
combined_df = rbind(archr, signac, scmacs)
combined_df = combined_df[TotalPeaks > 4000]

combined_df$Model = factor(combined_df$Model,
                           levels=c('ArchR', 'Signac','scMACS'))

pdf('panelC_runtime.pdf')
p <- ggplot(combined_df,
       aes(x=TotalPeaks,
           y=RunTime,
           col=Model,
           fill=Model)) + geom_point() +geom_line(size=2)+
  scale_y_continuous(trans='log10')+ylab('Runtime\n(Log10 Seconds)')+
  xlab('Peaks Tested')+theme_bw()+
  theme(legend.position = 'none', 
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        strip.text.x = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90),
        axis.text.y = element_text(size=14),
  )
print(p)

dev.off()

get_ratios <- function(combined_df, ModelComparison){
  
  df = dplyr::left_join(combined_df, 
                        combined_df[Model==ModelComparison], 
                                    by='TotalPeaks')
  
  df$Ratio = df$RunTime.y/df$RunTime.x
  ratios= df$Ratio[df$Model.x=='scMACS']
  return(ratios)
}

print(get_ratios(combined_df, 'Signac'))
print(get_ratios(combined_df, 'ArchR') )