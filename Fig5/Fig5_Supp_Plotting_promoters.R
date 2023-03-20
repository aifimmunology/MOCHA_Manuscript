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

### train LMMs for group 1
install.packages('pbapply')
require(pbapply)

################################################################################
################################################################################
yMat = g1_tsam[,which(zeroes1==0)]



res = rbindlist(res_list[which(sapply(res_list, class)!='try-error')])
res[is.nan(res$zi_pval)]$zi_pval <- 1

ggplot(res,
		 aes(x=-log(lmm_pval),
		 		 y=-log(zi_pval))) + geom_point()+
theme_minimal()+
theme(text=element_text(size=14))+
xlab('Linear Mixed Model')+
ylab('Zero Inflated LMM')+
ggtitle('Time -log(P values) Per Model')+
geom_abline()


ggplot(res,
		 aes(x=lmm_coef,
		 		y= zi_coef)) + geom_point()+
theme_minimal()+
theme(text=element_text(size=14))+
xlab('Linear Mixed Model')+
ylab('Zero Inflated LMM')+
ggtitle('Time Coefficients Per Model')+
geom_abline()



################################################################################

res_list = lapply(1:10,
									function(x)
									{
										print(x)
										X = g1_meta
										y = yMat[,x]
										X$y = y
										
										fit = glmmTMB(y ~ time +(1|PTID) ,
																	ziformula = ~0,
																	data=X,
																	family = gaussian(),
																	REML=F)

										plot(hist(y))
										print(summary(fit))
									})

################################################################################
################################################################################
X = rbind(g1_meta, g2_meta)
idx=sample(ncol(g1_tsam),120)
y = rbind(g1_tsam[,idx], g2_tsam[,idx])
y$Sample = row.names(y)
X= dplyr::left_join(X, cellCounts)

X= dplyr::left_join(X, y)
X=melt(X, id.vars = c('Sex','time','time_sqrd','Age','PTID','Sample','V2','cellTypeLabelList','N'))
X$Group = ifelse(X$Sample %in% group1_samples,
								 'group1','group2')
################################################################################
################################################################################

### Plot Cell Counts 
ggplot(X,
		 aes(x=N,
		 		y=value,
		 		#col=Group,
		 		alpha=0.5))+geom_point()+
xlab('Cell Counts')+ylab('Log2 Intensity')+
	theme_bw()+
	theme(title=element_text(size=20, hjust=0.5),
		text=element_text(size=12),
		legend.position = 'none',
		axis.text.x=element_text(size=12, angle=90),
		strip.text.x = element_blank())+
	facet_wrap(~variable,ncol=10)+
	geom_vline(xintercept = 300, lty=2, col='black')+
	ggtitle('Cell Counts vs. Log2 Intensity')+
	scale_x_log10(
	) 
tab1=table(X$N > 300, X$value>0)
tab1
fisher.test(tab1)


### Plot Age 
ggplot(X,
			 aes(x=Age,
			 		y=value))+geom_point()+
	xlab('Age')+ylab('Log2 Intensity')+
	theme_bw()+
	theme(text=element_text(size=12),
			 strip.text.x = element_blank())+
	geom_vline(xintercept = 40, lty=2, col='black')+
	facet_wrap(~variable,ncol=10)+
	ggtitle('Age vs. Log2 Intensity')
tab1 = table(X$Age > 40, X$value>0)
tab1
fisher.test(tab1)


ggplot(X,
			 aes(x=Sex,
			 		y=value))+geom_violin()+
	xlab('Age')+ylab('Log2 Intensity')+
	theme_bw()+
	theme(text=element_text(size=12),
				strip.text.x = element_blank())+
	facet_wrap(~variable, ncol=10)+
	ggtitle('Sex vs. Log2 Intensity')

tab1 = table(X$Sex, X$value>0)
tab1
fisher.test(tab1)

ggplot(X,
			 aes(x=time,
			 		y=value))+geom_point()+
	xlab('Age')+ylab('Log2 Intensity')+
	theme_bw()+
	theme(text=element_text(size=12),
				strip.text.x = element_blank())+
	facet_wrap(~variable, ncol=10)+
	ggtitle('Time vs. Log2 Intensity')
tab1= table(X$time>0, X$value>0)
tab1
fisher.test(tab1)



########################################################################
########################################################################
### Test fisher test / ODDs ratio 
### of 0s vs. covariate 

test_features <- function(i){
	
	X = rbind(g1_meta, g2_meta)
	idx=sample(ncol(g1_tsam),1)
	y = data.frame(c(g1_tsam[,idx], (g2_tsam[,idx])))
	y$Sample = c(row.names(g1_tsam), 
							 row.names(g2_tsam))
	X= dplyr::left_join(X, cellCounts)
	
	X= dplyr::left_join(X, y)
	X=melt(X, id.vars = c('Sex','time','time_sqrd','Age','PTID','Sample','V2','cellTypeLabelList','N'))
	X$Group = ifelse(X$Sample %in% group1_samples,
									 'group1','group2')
	if(sum(y==0)>0){
		res = data.frame(
			Time_pval= fisher.test(X$time>0, X$value>0)$p.value,
			Cells_pval= fisher.test(X$N > 300, X$value>0)$p.value,
			Age_pval= fisher.test(X$Age > 40, X$value>0)$p.value,
			Sex_pval= fisher.test(X$Sex, X$value>0)$p.value,
			Time_ODDs= fisher.test(X$time>0, X$value>0)$estimate,
			Cells_ODDs= fisher.test(X$N > 300, X$value>0)$estimate,
			Age_ODDs= fisher.test(X$Age > 40, X$value>0)$estimate,
			Sex_ODDS= fisher.test(X$Sex, X$value>0)$estimate
			)
		
		return(res)
	} else{
		return(data.frame(
			Time_pval= NA,
			Cells_pval= NA,
			Age_pval= NA,
			Sex_pval= NA,
			Time_ODDs= NA,
			Cells_ODDs= NA,
			Age_ODDs= NA,
			Sex_ODDS= NA
					)	)
	}
	
	
}

res_list = mclapply(1:49679,
				 function(x)
				 test_features(x),
				 mc.cores=8
)
res =rbindlist(res_list)
res = melt(res)
res$value[is.infinite(res$value)] <- 50

ggplot(res[value <51],
			 aes(x=value))+geom_histogram(bins=30)+
	facet_wrap(~variable, scales='free', ncol=4)

summary_res = res[, quantile(value, 0.4, na.rm=T), 
		by=variable]

require(dplyr)
summary_res %>% arrange(V1)

################################################################################
################################################################################

### Show an example promoter 

zeroes = colSums(tsam==0)
idx = which(zeroes > 0)
set.seed(1001)
examples = sample(idx, 3)
example_promoter = data.frame(tsam[,examples])
example_promoter = melt(example_promoter)
pdf('~/Desktop/Figure 3/0s_justification.pdf', height=4, width=10)
ggplot(example_promoter,
			 aes(x=value)
)+geom_histogram() + theme_minimal()+
	xlab('Log2 Promoter Accessibility')+
	ggtitle('0s in Log2 Promoter Accessibility')+
	theme(title = element_text(size=14, hjust = 0.5))+
	facet_wrap(~variable, ncol=5, scales='free')
dev.off()

### Show total # of Zeroes 

plot_heatmap <- function(tsam){
	tsam_mlt = data.table(melt(tsam, id.var='Sample'))
	tsam_mlt[tsam_mlt==0] <- NA
	#samplePromoters = sample(unique(tsam_mlt$variable), 490)
	samplePromoters = unique(tsam_mlt$variable)
	ggplot(tsam_mlt[ variable %in% samplePromoters,],
				 aes(y=Sample,
				 		x= reorder(variable, value, median),
				 		fill=value))+geom_tile()+
		theme(title=element_text(size=16),
					axis.text.y = element_text(size=11),
					axis.text.x = element_blank())+
		scale_fill_gradient2(low='purple', mid = 'black', high='yellow', midpoint=12,
												 na.value = 'white')+
		xlab('Promoter Tile')+
		ylab('Sample')
	
}

pdf('example_heatmap.pdf', height=12, width=9)
plot_heatmap(tsam)
dev.off()
