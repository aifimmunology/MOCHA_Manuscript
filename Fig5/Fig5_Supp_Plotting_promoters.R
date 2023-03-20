###############################################################################
###############################################################################
### GLMMTMBs
###
###
###############################################################################

###############################################################################
###############################################################################
rm(list=ls())
require(glmmTMB)
require(data.table)
require(MOCHA)
require(ggplot2)
require(parallel)
require(lmerTest)
require(dplyr)

### load data 
stm = readRDS('~/Downloads/SampleTileObject-2.rds')
model_results = fread('~/Downloads/gene_hits (5).csv')
### extract CD16 Promoters 
tmp <- plyranges::filter(SummarizedExperiment::rowRanges(stm), tileType == 'Promoter')
cd16_promoters = tmp[tmp$`CD16 Mono`==T]
cd16_tsam = MOCHA::getCellPopMatrix(stm, cellPopulation = 'CD16 Mono')
cd16_promoter_tsam <- cd16_tsam[row.names(cd16_tsam) %in% MOCHA::GRangesToString(cd16_promoters),]

### Extract Metadata 
meta = as.data.table(stm@colData)

###############################################################################
###############################################################################

### 
meta$time <- meta$days_since_symptoms-mean(meta$days_since_symptoms, na.rm=T)
meta$time_sqrd <- (meta$time)^2
numCores = 1

### Establish group1 
### and group 2
Covid_subjects <- c(32209, 31207, 32140, 32245, 32054, 32255, 32131, 31945, 32416, 42409,
								 32220, 32038, 31874, 32124, 31924, 32251, 32196, 32415)   

### Extract samples
### from group 1 and 2
covid_samples <- meta$Sample[meta$PTID %in% Covid_subjects]

### create separate 
### TSAMs for each 

tsam = cd16_promoter_tsam[, colnames(cd16_promoter_tsam) %in% covid_samples]
tsam = as.data.frame(t(tsam)); tsam$Sample = row.names(tsam)

meta <- meta[Sample %in% tsam$Sample, c('Sex','time','time_sqrd','Age','PTID','Sample')]

cellCounts = data.table(stm@metadata$CellCounts)
cellCounts = cellCounts[cellTypeLabelList=='CD16 Mono' & V2 %in% covid_samples]
cellCounts$Sample = cellCounts$V2

###############################################################################
###############################################################################

###############################################################################
###############################################################################

plot_examples <- function(i, meta, tsam, gene){
	X = meta
	y = tsam[,i]
	X = left_join(X, cellCounts)
	X$y = y
	X$days = X$time + abs(min(X$time))
	

	if(sum(y==0)==0){
		
		### Specify no zero-inflation
		### if Y is fully continuous
		fit = try(glmmTMB(y ~  Age + Sex + 
												time +  (1|PTID) ,
											ziformula = ~ 0,
											data=X,
											family = gaussian(),
											REML=T))
		
	} else{
		fit = try(glmmTMB(y ~  Age + Sex +
												time + (1|PTID) ,                     
											ziformula = ~ 0+N,
											data=X,
											family = gaussian(),
											REML=T))
	}
	
	sum_fit = summary(fit)
	print(sum_fit)
	coefs = sum_fit$coefficients$cond[,1]
	X$Prediction = coefs[1] + 
								 coefs[2]* mean(X$Age) + 
								 coefs[4]* X$time
	
	p = ggplot(X,
				 aes(x=days,
				 		y=y, 
				 		col=as.character(PTID)))+geom_point(size=0.8,
				 																				alpha=0.4)+
		geom_line(data=X[X$y>0],
								aes(x=days,
										y=y,
										alpha=0.3),
			linewidth=0.5)+ggtitle(paste(gene, i, sep='\n'))+
		geom_line(data=X,
							aes(x=days,
									y=Prediction),
							linewidth=1,
							col='black')+
		theme_minimal()+
		xlim(0,100)+
		xlab('')+
		ylab('')+
		theme(legend.position = 'none')
	
	X$Gene = gene
	X$Promoter = i

	return(list(Plot=p,
							DF=X))

}

model_results = model_results %>% arrange(cont_pval)

pdf('~/Desktop/Figure 6/top6_hits.pdf', height=5, width=10)
ggpubr::ggarrange(plotlist = lapply(1:6, function(x)
	plot_examples(model_results$Promoter[x],meta, tsam, model_results$Gene[x])$Plot),
	nrow=2, ncol=3)
dev.off()

example_data =	lapply(1:6, function(x)
	plot_examples(model_results$Promoter[x],meta, tsam, model_results$Gene[x])$DF
	)

example_data = rbindlist(example_data)
write.csv(example_data,
					'~/Desktop/Figure 6/top_gene_promoter_hits.csv')
################################################################################
################################################################################

