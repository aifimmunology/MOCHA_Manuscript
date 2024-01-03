### Plot Results
library(tidyverse)
library(ggplot2)

herit <- read.table('sumstats/AllResults.tsv', header = TRUE)

head(herit)

## There was some oddities related to hise. 
herit <- dplyr::filter(herit, !grepl('_Results-checkpoint.results',file))

herit <- dplyr::mutate(herit, 
        file = gsub("_Results.results", "", file))


## Check which results are consistent across studies for the same sample. 
dplyr::filter(herit, file == 'MOCHA_CD8_TEM_intensities_lp' & Category == 'Coding_UCSCL2_1')

dim(herit)
#10584 12
summary(herit$Enrichment_p)
hist(herit$Enrichment_p)
#P-values show real signal. 

##Adjust p-value
herit$Enrichment_padj <- p.adjust(herit$Enrichment_p, method = 'fdr')

## Let's split out by cell type, cohort, and method

herit <- dplyr::mutate(herit, 
    Method = sub('_.*', '', file),
    CellType = sub("HOMER_|MOCHA_|MACS2_", "", sub('_intensities_.*', '', file)), 
   Cohort = sub(".*intensities_", "", file))
        
heritf <- dplyr::filter(herit, Enrichment_padj < 0.05)

table(heritf$Cohort)
#   bm covid    lp 
#  230   229   226 
table(heritf$Method)
#HOMER MACS2 MOCHA 
#  228   226   231

#Look at reproducibility
heritf_sum1 <- dplyr::group_by(heritf, Category, GWAS_Study, CellType, Cohort) %>%
dplyr::summarize(Number = dplyr::n(),
                 Methods = paste0(Method, collapse =', ')) 

### Some of the models are mispecified. 
### Removed negative enrichment values. 
tmp <- dplyr::filter(heritf, Enrichment < 0) %>% dplyr::select(Category, CellType, Cohort, GWAS_Study) %>% table() 
sum(!as.data.frame(tmp)$Freq %in% c(0,3))
##### The methods all fail for a given celltype/cohort/study combination
## Let's filter them out. 
heritf2 <- dplyr::filter(heritf, 
                    Enrichment > 0)

heritf_sum2 <- dplyr::group_by(heritf2, 
                    Category, GWAS_Study, CellType, Cohort) %>%
    dplyr::summarize(Number = dplyr::n(),
                 Methods = paste0(Method, collapse =', ')) 

## Percent reproducible
sum(heritf_sum2$Number == 3)/length(heritf_sum2$Number)
#95.43% reproducible

uniqueHits <- dplyr::filter(heritf_sum2, Number < 3) # Only 9 groups
uniqueHits %>% ungroup() %>%
    dplyr::select(Methods) %>% unlist() %>% table()
## When it was not detected across all three methods, 
#5/9 were only detected in MOCHA. 
#7/9 were MOCHA onlu, or MOCHA & HOMER
#The remaining 2 are MACS2 and HOMER

## Now create plots. 
## First show degree of agreement across methods
## Then highlight disagreement


pdf('Heritability_Summary.pdf')

ggplot(heritf2, aes(x = Method, fill = CellType)) +
    geom_bar() + theme_bw()

ggplot(heritf_sum2, aes(x = GWAS_Study, 
                        group = as.character(Number), fill = as.character(Number)))  +
    geom_bar() + theme_bw() + 
    ylab('Number of Significant Hits for LDSC Enrichment') +
    xlab('GWAS Study') +
    ggtitle("Peaksets for All Three Methods agree for 96% of Significant Categories")

### log10 of Enrichment_p-value (adjusted)
### vs x-axis of trait, facet_wrap by study

dplyr::mutate(herit, log10Padj = -log10(Enrichment_padj)) %>%
    dplyr::filter(paste(Category,GWAS_Study,CellType,Cohort, sep = '_') %in%
                  paste(uniqueHits$Category, uniqueHits$GWAS_Study,
                        uniqueHits$CellType, uniqueHits$Cohort, sep = '_')) %>%

    ggplot(., aes(x = Category, y= log10Padj, color = Method)) +
        geom_jitter() + theme_bw() + theme( axis.text.x = element_text(angle=90,hjust=1)) +
        geom_hline(yintercept = -log10(0.05), color = 'red') + facet_wrap(~GWAS_Study, scales = 'free_x')

dplyr::filter(heritf2, paste(Category,GWAS_Study,CellType,Cohort, sep = '_') %in%
                  paste(uniqueHits$Category, uniqueHits$GWAS_Study,
                        uniqueHits$CellType, uniqueHits$Cohort, sep = '_')) %>%
    dplyr::group_by(Method) %>% dplyr::summarize(Number = dplyr::n()) %>%
    ggplot(., aes(x = "", y = Number, fill = Method)) +
    geom_bar(stat = 'identity',width=1, color="white") + coord_polar("y", start =0) +
    theme_bw() + 
    ggtitle("Significantly Enriched Categories that do not replicate across Methods.")

uniqueHits %>%
    dplyr::group_by(Methods) %>% dplyr::summarize(Number = sum(Number)) %>%
    dplyr::mutate(Methods = ifelse('HOMER, MOCHA' == Methods, 'MOCHA, HOMER', Methods)) %>%
    dplyr::ungroup() %>%
    arrange(desc(Methods)) %>% 
    mutate(prop = Number / sum(.$Number) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
    ggplot(., aes(x = "", y = prop, fill = Methods)) +
    geom_bar(stat = 'identity',width=1, color="white") + coord_polar("y", start =0) +
    theme_void() + 
    ggtitle("Significantly Enriched Categories that do not replicate across Methods.") +
    scale_fill_manual(values = c('MOCHA' = 'deepskyblue3', 
                                 'MOCHA, HOMER' = 'cornflowerblue', 
                                 'MACS2, HOMER' = '#484848')) +
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = Methods), color = "white", size=6) + 
     geom_text(aes(y = ypos, label = paste(as.character(round(prop,0)), '%',sep='')), 
            color = "white", size=6, vjust = 3) 

uniqueHits %>%
    dplyr::group_by(Methods) %>% dplyr::summarize(Number = sum(Number)) %>%
    dplyr::mutate(Methods = ifelse('HOMER, MOCHA' == Methods, 'MOCHA & HOMER Peaksets', Methods)) %>%
 dplyr::mutate(Methods = ifelse('MACS2, HOMER' == Methods, 'MACS2 & HOMER Peaksets', Methods)) %>%
 dplyr::mutate(Methods = ifelse('MOCHA' == Methods, 'MOCHA Peakset', Methods)) %>%
    dplyr::mutate(MethodType = ifelse(grepl("MOCHA", Methods), 'Significant with MOCHA Peaksets', 'Significant by Other Peaksets')) %>% 
    ggplot(., aes(x =MethodType, y = Number, fill = Methods)) + 
     geom_bar(stat = 'identity') + theme_bw() + 
     scale_fill_manual(values = c('MOCHA Peakset' = 'deepskyblue1', 
                                 'MOCHA & HOMER Peaksets' = 'deepskyblue3', 
                                 'MACS2 & HOMER Peaksets' = '#484848')) +
    xlab(NULL)+
    ylab('Number of Significantly Enriched Categories')

### Plot top categories across cell types

dplyr::group_by(heritf_sum2, Category, GWAS_Study) %>% 
    dplyr::select(Category, GWAS_Study, CellType) %>% 
    dplyr::filter(GWAS_Study != 'Abnormal_Immune_Study') %>%
    distinct() %>% dplyr::summarize(NumberOfCellType = round(dplyr::n()/7,2)) %>% 
    dplyr::ungroup() %>% 
    arrange(NumberOfCellType) %>% 
    dplyr::mutate(Category = factor(Category, levels = unique(.$Category))) %>% as.data.frame() %>%
    ggplot(., aes(y = Category, x = NumberOfCellType, fill = GWAS_Study)) +
     geom_col(position = 'dodge') + theme_bw()  + theme(legend.position = 'bottom') + 
    xlab('Percentage of Cell Types Test') +ylab(NULL)


### Plot top categories across cell types
dplyr::group_by(heritf_sum2, Category, GWAS_Study) %>% 
    dplyr::select(Category, GWAS_Study, CellType) %>% 
    dplyr::filter(GWAS_Study != 'Abnormal_Immune_Study') %>%
    distinct() %>% dplyr::summarize(NumberOfCellType = round(dplyr::n()/7,2)) %>% 
    dplyr::ungroup() %>% 
    arrange(NumberOfCellType) %>% 
    dplyr::mutate(Category = factor(Category, levels = unique(.$Category))) %>% as.data.frame() %>%
    ggplot(., aes(y = Category, x = NumberOfCellType, fill = GWAS_Study)) +
     geom_col(position = 'dodge') + theme_bw()  + theme(legend.position = 'bottom') + 
    xlab('Percentage of Cell Types Test') +ylab(NULL)
dev.off()
