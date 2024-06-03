#######################################################################################
#######################################################################################

### Generates graph
### for Figure 3 cutoff 
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
require(mixtools)
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
mocha = readRDS("/home/jupyter/MOCHA_Manuscript2/Fig2/Covid-19/MOCHA.RDS")
cd16 = mocha[['CD16 Mono']]

### identify positive & negative samples 
positive = cd16[,colnames(cd16) %in% covid_pos] 
negative = cd16[,colnames(cd16) %in% covid_neg] 

covidPeaks = compactAssay(positive, 'TotalIntensity')
controlPeaks = compactAssay(negative, 'TotalIntensity')

SampleTileMatrices <- getSampleTileMatrix( 
    mocha,
    cellPopulations = "CD16 Mono",
    groupColumn = "COVID_status",
    threshold = 0.2,
    verbose=T
)
tsam  = SampleTileMatrices@assays@data$`CD16 Mono`
tsam[is.na(tsam)] <- 0

tsam = log2(tsam+1)

mean_intensity =rowMeans(tsam)

require(mixtools)
mixture_model = normalmixEM(mean_intensity,
            lambda=0.5, mu=c(4.5,12),
            sigma=1)


pdf('/home/jupyter/MOCHA_Manuscript2/SuppFig4_ZI_differentials/differentials_cutoff.pdf')
plot(mixture_model, density=T, breaks=200)
dev.off()

zeroes =rowMeans(tsam==0)

summary(zeroes[mean_intensity > 12])
summary(zeroes[mean_intensity < 12])

df = data.table(
    zeroes=zeroes,
    intensity=mean_intensity)

require(scattermore)
png('/home/jupyter/MOCHA_Manuscript/Fig3/zeroes.png')
ggplot(df,
       aes(x=intensity,
           y=zeroes))+geom_point()+geom_smooth()+
        theme_minimal()+
        theme(text=element_text(size=14))+
        geom_vline(xintercept=12)
dev.off()


df$Group = ifelse(df$intensity < 12, 'Below cutoff', 'Above Cutoff')
pdf('/home/jupyter/MOCHA_Manuscript/Fig3/zeroes_violin.pdf')
ggplot(df,
       aes(x=Group,
           y=zeroes))+geom_violin(scale='width')+
        theme_minimal()+
        theme(text=element_text(size=14))+
        geom_vline(xintercept=12)+
        xlab('Tiles') + ylab('% of Samples without Fragments Across Tiles')
dev.off()


pdf('/home/jupyter/MOCHA_Manuscript/Fig3/zeroes_violin.pdf')
ggplot(df,
       aes(x=zeroes,
          fill=Group))+geom_histogram(position='identity', alpha=0.4)+
        theme_minimal()+
        theme(text=element_text(size=14))+
        ylab('# Tiles') + xlab('% of Samples without Fragments Across Tiles')
dev.off()


#######################################################################################
#######################################################################################
chr1 = tsam[grep('chr1:', row.names(tsam)),]

avg_zero_chr1 = rowMeans(chr1==0)
avg_zero = rowMeans(chr1==0)
avg_int = rowMeans(chr1)

zero_df = data.table(
    Perc_zero = avg_zero_chr1,
    Name = names(avg_zero_chr1),
    Location = gsub('chr1:', '', names(avg_zero_chr1))
    )

zero_df$Location = as.numeric(gsub('-[0-9]*','',zero_df$Location))

zero_df = zero_df %>% arrange(Location)

pdf('/home/jupyter/MOCHA_Manuscript2/SuppFig4_ZI_differentials/zeros_chr1.pdf')

p= ggplot(zero_df,
       aes(x=Location,
           y=Perc_zero))+
            geom_point(size=0.3)+
        theme(axis.text.x=element_text(size=5, angle=90))#+
        #scale_x_continuous(labels=y, breaks=0, 10000000, 20000000)
print(p)
dev.off()

chr1 = data.frame(chr1)
chr1$Location = row.names(chr1)
chr1 = data.table(chr1)
chr1_melted = melt(chr1)

#chr1_melted[chr1_melted==0] <- NA
pdf('/home/jupyter/MOCHA_Manuscript2/SuppFig4_ZI_differentials/chr1_hm.pdf')
chr1_melted$Sample = as.numeric(factor(chr1_melted$variable))
ggplot(chr1_melted,
       aes(x=reorder(Location,  value, mean),
           y=Sample,
          fill=value))+
        geom_tile()+scale_fill_gradient2(low='light blue',                                         high='dark blue', na.value='white')+
        theme(axis.text.x=element_blank())+
        xlab('Chr1 Tiles') + ylab('Sample')
dev.off()


zero_df_bins  = data.table(
    Perc_zero = avg_zero,
    avg_int = round(rowMeans(chr1[,colnames(chr1) !='Location', with=F]))
    )

zero_df_bins_means = zero_df_bins[, list(
                                        Mu=mean(Perc_zero)*100, 
                                        Sigma = sd(Perc_zero)*100),
                                  by=avg_int]

zero_df_bins_means = zero_df_bins_means %>% arrange(avg_int)

pdf('/home/jupyter/MOCHA_Manuscript2/SuppFig4_ZI_differentials/zeros_bin.pdf')
zero_df_bins_means$upper = zero_df_bins_means$Mu + zero_df_bins_means$Sigma
zero_df_bins_means$lower = zero_df_bins_means$Mu - zero_df_bins_means$Sigma

zero_df_bins_means$lower[zero_df_bins_means$lower < 0] <- 0

p = ggplot(zero_df_bins_means,
       aes(x=avg_int,
           y=Mu))+
           geom_bar(stat='identity')+
                geom_errorbar(aes(ymin=lower, ymax=upper))+
            theme(text=element_text(size=12))+
        theme_minimal()+
        xlab('Average Log2 Intensity')+
        ylab('Percent 0s in tile')+
        geom_vline(xintercept=12)
        
print(p)
dev.off()

pdf('/home/jupyter/MOCHA_Manuscript2/SuppFig4_ZI_differentials/suppFig2F.pdf')

ggplot(zero_df_bins,
       aes(x=avg_int, 
           y = Perc_zero))+
    geom_jitter(alpha = 0.5, width = 0.2)+
    geom_smooth()+
   xlab('Log2(Accessibility+1)')+
        ylab('percent of zeroes')+
    theme_minimal()

dev.off()

#######################################################################################
#######################################################################################

mocha_dats = fread('../Fig3/generateData/cd16_mocha.csv')
mocha_dats = mocha_dats %>% arrange(FDR)

mocha_dats$diffs0s = abs(mocha_dats$Pct0_Case - mocha_dats$Pct0_Control)
zero_diffs_dats = which(mocha_dats$diffs0s > 0.5)

zero_diffs_dats = mocha_dats[zero_diffs_dats[1:5],]
high_log2FC = mocha_dats[1:5,]

example_dats = melt(tsam[row.names(tsam) %in% c(zero_diffs_dats$Tile[1:5], high_log2FC$Tile[1:5]),])

colnames(example_dats) = c('Tile','Sample','Log2_l1')

example_dats$Group = 'COVID+'
example_dats$Group[example_dats$Sample %in% covid_neg] = 'COVID-'

pdf('../SuppFig4_ZI_differentials/example_dat1.pdf', width=6, height=8)
ggplot(example_dats[example_dats$Tile == high_log2FC$Tile[1],],
      aes(x=Log2_l1,
          fill=Group))+geom_histogram()+
        facet_wrap(~Group, scales='free_y', ncol=1 )+
        theme_minimal()+
        xlab('Log2 Normalized Frag Counts')+
        ylab('Distribution')+
        theme(text=element_text(size=14))+
        ggtitle(high_log2FC$Tile[1])
dev.off()

pdf('../SuppFig4_ZI_differentials/example_dat2.pdf', width=6, height=8)
ggplot(example_dats[example_dats$Tile == high_log2FC$Tile[2],],
      aes(x=Log2_l1,
          fill=Group))+geom_histogram()+
        facet_wrap(~Group, scales='free_y', ncol=1 )+
        theme_minimal()+
        xlab('Log2 Normalized Frag Counts')+
        ylab('Distribution')+
        theme(text=element_text(size=14))+
        ggtitle(high_log2FC$Tile[2])
dev.off()



pdf('../SuppFig4_ZI_differentials/example_dat3.pdf', width=6, height=8)
ggplot(example_dats[example_dats$Tile == zero_diffs_dats$Tile[1],],
      aes(x=Log2_l1,
          fill=Group))+geom_histogram()+
        facet_wrap(~Group, scales='free_y', ncol=1 )+
        theme_minimal()+
        xlab('Log2 Normalized Frag Counts')+
        ylab('Distribution')+
        theme(text=element_text(size=14))+
                ggtitle(zero_diffs_dats$Tile[1])

dev.off()


pdf('../SuppFig4_ZI_differentials/example_dat4.pdf', width=6, height=8)
ggplot(example_dats[example_dats$Tile == zero_diffs_dats$Tile[2],],
      aes(x=Log2_l1,
          fill=Group))+geom_histogram()+
        facet_wrap(~Group, scales='free_y', ncol=1 )+
        theme_minimal()+
        xlab('Log2 Normalized Frag Counts')+
        ylab('Distribution')+
        theme(text=element_text(size=14))+
                ggtitle(zero_diffs_dats$Tile[2])

dev.off()




#### Write results to file for 
#### creating 'source data' 
setwd('../SuppFig4_ZI_differentials/')

write.csv(example_dats, 
          file='exemplar_differentials.csv'
          )

write.csv(zero_df_bins_means,
          file='zeros_against_intensity.csv')

write.csv(chr1,
          file='chr1_heatmap.csv')

write.csv(mean_intensity,
          file='avg_log2_intensity_for_MixModels.csv')