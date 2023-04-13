##############################################################################
##############################################################################

## Figure 5: analyses 
## 

##############################################################################
##############################################################################
rm(list=ls())
#devtools::install('~/Desktop/MOCHA/')
setwd('~/Desktop/Figure 6/')
source('figure5_helper_functions.R')

### load libraries 
require(data.table)
require(ggplot2)
require(RaggedExperiment)
require(MOCHA)
library(lmerTest)
source('figure5_helper_functions.R')

### load data objects 
Obj <- readRDS('~/Downloads/STM_ChromVAR.rds')

### We Center the 
### time variable in order 
### to minimize correlation
### between time and time^2
### ---> Ch8 Applied Linear Reg (Kuntner et al)

Obj$time <- Obj$days_since_symptoms-mean(Obj$days_since_symptoms, na.rm=T)
Obj$time_sqrd <- (Obj$time)^2
numCores = 1

##############################################################################
##############################################################################
### Define linear formula
linFormula <- exp ~ Age + sex_at_birth + time + (1|PTID)
quadFormula <- exp ~ Age + sex_at_birth + time + time_sqrd + (1|PTID)

##############################################################################
##############################################################################

### Get all covid
### subjects in FH
### Covid Cohort
Group= c(32209, 31207, 32140, 32245, 32054, 32255, 
				32131, 31945, 32416, 42409, 32220, 32038, 
				31874, 32124, 31924, 32251, 32196, 32415)   

dev_G1 <- subsetDev(Obj,
                    subsetBy =  'PTID', groupList = Group)


linearModel <- modelDeviations(dev_G1 ,formula = linFormula,
														 type = 'z', numCores = 4)

quadraticModel <- modelDeviations(dev_G1 ,formula = quadFormula,
															 type = 'z', numCores = 4)

quadraticResults <- lapply(1:length(quadraticModel),
														 function(x){
														 	fit = quadraticModel[[x]]
														  summary_fit = summary(fit)	
														  res = c(summary_fit$coefficients[,1],
														  				summary_fit$coefficients[,5])
														 	res
														 })
## Manipulate the 
quadraticResults = do.call(rbind,quadraticResults)
colnames(quadraticResults)[6:10] = paste(colnames(quadraticResults)[6:10],'_pval',sep='')
quadraticResults = data.frame(quadraticResults)
quadraticResults$FDR = p.adjust(quadraticResults$time_pval, 'fdr')
quadraticResults$TimeSignificant = quadraticResults$FDR < 0.1
quadraticResults$TF = rownames(Obj)

pdf('Volcano_TFs.pdf', height=9, width=9)
ggplot(quadraticResults,
			 aes(x=time,
			 		 y=-log10(FDR),
			 		 label=TF)) + geom_point()+
	theme_minimal()+
	ggrepel::geom_text_repel(
		size=2.5
	)+
	theme(title=element_text(size=18))+
	xlim(NA,0.11)+
	xlab('Slope') + ylab('-log10(FDR)')
dev.off()
	
pvalues = sapply(1:length(linearModel),
								 function(x)
								 	as.numeric(try(anova(linearModel[[x]], quadraticModel[[x]])$Pr[2])
								 	)
)

hist(pvalues, main = 'Quadratic Term in ChromVar Deviations')

z_scores = data.frame(t(dev_G1@assays@data$z))
z_scores$Sample = row.names(z_scores)
meta = as.data.table(dev_G1@colData)
X = meta[,c('Sample','PTID','time','time_sqrd','Age','Sex')]
##############################################################################
##############################################################################

require(dplyr)
quadraticResults= quadraticResults %>% arrange (time_pval)

plot_examples<- function(tf_name){
	mat = dplyr::left_join(X, z_scores[,c(tf_name,'Sample')])
	mat$PTID = factor(mat$PTID)
	mat= melt(mat, measure.vars = tf_name)
	coefficients=quadraticResults[quadraticResults$TF==tf_name,]
	
	
	mat$PopulationPrediction = as.numeric(coefficients[1]) + 
		as.numeric(coefficients[2])*mean(mat$Age)+
		mat$time * as.numeric(coefficients[4])+
		mat$time_sqrd * as.numeric(coefficients[5])
	
	mat = mat %>% arrange(time)
	mat$time = mat$time+39
	#mat = melt(mat, measure.vars = c('male','female'))

	a = ggplot(mat,
				 aes(x=time,
				 		y=value,
				 		col=PTID)) +geom_line(linewidth=1, alpha=0.3)+
		theme_minimal()+
		geom_point(data=mat,
							 aes(x=time,
							 		y=value
							 		))+
		geom_line(data=mat,
							aes(x=time,
									y=PopulationPrediction
							), linewidth=2, col='black')+		
				ggtitle(tf_name)+
		theme(legend.position = 'none',
					title=element_text(size=24),
					axis.text=element_text(size=12))+
		xlim(0,100)+
		ylab('Z-score')+
		xlab('days')
	
	return(list(Plot = a,
							Df= mat))
	
	
}


### plot examples
tfs = c('JUN','ATF7','IRF4','HIVEP1','FOXO1','ATF5')
plotList=lapply(tfs,
			 function(x) 
			 	plot_examples(x)$Plot)
			 
pdf('lineplots.pdf',
		width=9, height=6)
ggpubr::ggarrange(plotlist = plotList,
									ncol=3,nrow=2)
dev.off()


### extract data
matList=lapply(tfs,
								function(x) 
									plot_examples(x)$Df)
exemplar_TFs = rbindlist(matList)

write.csv(exemplar_TFs,
					file='exemplar_TFs.csv')

write.csv(quadraticResults,
					file='FullCovid_TF_results.csv')
