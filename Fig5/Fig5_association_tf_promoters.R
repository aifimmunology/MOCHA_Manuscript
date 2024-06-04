##############################################################################
##############################################################################

## Figure 5: analyses 
## 
##   - integrate chromVar
##   - with Promoter Accessibility
##   - via correlations

##############################################################################
##############################################################################
#install.packages('glmmTMB')
rm(list=ls())
require(data.table)
require(parallel)
require(stringi)
require(stringr)
require(ggplot2)
require(glmmTMB)
setwd('/home/jupyter/MOCHA_Manuscript/Fig5/Fig5_integration')

## load data
pathway_hits = fread('covid_pathway_hits_annotated_v2.csv')
tsam = fread('TSAM_promoters.csv')
colnames(tsam)[1] = 'Sample'
tf_mat = read.csv('tf_matrix.csv')
colnames(tf_mat)[1] = 'Sample'
res_prom = fread('../gene_hits.csv')
res_tf = fread('TF_quadratic_model.csv')

### Extract Significant Genes 
unique_genes <- unique(res_prom$Gene[res_prom$FDR <0.1])

### some promoters span
### multiple genes 
### so we extract them 

group1_genes = sapply(unique_genes,
                      function(x)
                          unlist(strsplit(x, ', '))
                      )
unlist(group1_genes) -> genes

##### Predict Pathway-TF association 

predict_pathways <- function(pathway, plot_hits=F, tf='JUN'){
    
	pathway_genes = strsplit(innate_immune_pathways[description==pathway]$userId,
								';')[[1]]
	
    ## 
	promoters = res_prom[intersect(grep(paste(pathway_genes, collapse = '|'), res_prom$Gene),
                                   which(res_prom$FDR < 0.1)),]
	
	pathway_tsam = tsam[, promoters$Promoter, with=F]
    tf_vector = tf_mat[,tf]
	
	predict_tf = function(gene, plot=F, tf_vector=tf_vector){
		df =  data.frame(
			promoter= pathway_tsam[,gene,with=F][[1]], 
			TF =tf_vector
		)
		
		if(sum(df$promoter==0)>0){
			fit =glmmTMB::glmmTMB(promoter~TF,family = gaussian(),
														ziformula = ~.,
														data=df)
		} else { 
			fit =glmmTMB::glmmTMB(promoter~TF,
														family = gaussian(),
														ziformula = ~0,
														data=df)
		}
		
		sum_fit = summary(fit)
		model_coefs = sum_fit$coefficients$cond[,1]
		df$prediction= model_coefs[1]+
			model_coefs[2]*df$TF
		
		if(plot){
			p=ggplot(df[df$promoter>0,],
							 aes(x=TF,
							 		y=promoter))+geom_point(size=0.5)+
				geom_line(data=df,
									aes(x=TF,
											y=prediction),
									linewidth=2)+
				theme_minimal()+
				ylab(promoters$Gene[gene])+
				xlab(tf)
			return(p)
		}
		
		res = as.data.frame(t(sum_fit$coefficients$cond['TF',]))
		res$Gene = promoters$Gene[gene]
		res$Promoter=promoters$Promoter[gene]
		return(res)
	}
	
	res_p1 = mclapply(1:length(promoters$Gene),
									function(x) predict_tf(x,plot=F,
                                                           tf_vec=tf_vector),
                      mc.cores=50
	)
	
	res_p1=rbindlist(res_p1)
    res_p1$Pathway=pathway
    res_p1$Unique_genes = length(pathway_genes)
	res_p1$Pval = res_p1$`Pr(>|z|)`			
    res_p1$TF = tf
	
	if(plot_hits){
		p=ggpubr::ggarrange(plotlist=lapply(which(res_p1$`Pr(>|z|)`<0.05),
					function(x) predict_tf(x, plot=T
                                                                                                  )
                                            
                                            )
                                            )
        
        print(p)
	}

	return(res_p1)	
}

#####
require(dplyr)           
res_tf$FDR = p.adjust(res_tf$time_pval, 'fdr')
top_tfs = res_tf[FDR < 0.1,]
top_positive = top_tfs[time>0] %>% arrange(-time)
top_negative = top_tfs[time<0] %>% arrange(time)

res_tf$TimeSignificant = res_tf$FDR < 0.1
write.csv(res_tf,
          'tf_usage_model_results.csv')

## expand list to entire innate immune pathways 
innate_immune_pathways = pathway_hits[Level2=='Innate Immune System']
rm(pathway_hits)

## remove non innate immune
## pathway hits 

tfs_to_plot = c(top_positive$TF[1:5],
                top_negative$TF[1:5])


tf_pathway_interactions=lapply(tfs_to_plot,
       function(y)
           rbindlist(lapply(innate_immune_pathways$description,
			 function(x) 
			 	predict_pathways(x, plot_hits = F,tf=y
                                ))
                                )
                      
)


quantify_pathway <- function(tf_pathway_list, pathway){
    
    genes_tested = strsplit(innate_immune_pathways[description==pathway]$userId,
								';')[[1]]
    
    tmp = tf_pathway_list[tf_pathway_list$Pathway==pathway,]
    tmp = tmp[tmp$Pval < 0.05]
    gene_hits = unique(unlist(strsplit(unique(tmp$Gene),', ')))
    sum(gene_hits %in% genes_tested) / length(genes_tested)
}

tf_percentages = lapply(1:length(tf_pathway_interactions),
               function(y) 
              data.frame(
                  pathway = innate_immune_pathways$description,
                  TF = rep(tfs_to_plot[y],18),
    Perc_gene = sapply(1:length(innate_immune_pathways$description),
       function(x)
           quantify_pathway(tf_pathway_interactions[[y]], innate_immune_pathways$description[x]
                           )
                       )
    )
               )
                               
tf_percentages = rbindlist(tf_percentages)
                               
           
pdf("tf_pathway_results2.pdf", width=9, height=6)
ggplot(tf_percentages,
       aes(y=reorder(pathway, Perc_gene,mean),
           x=Perc_gene*100,
          fill=TF,
          col=TF))+geom_bar(stat='identity', position='dodge')+
    xlab('Percent of Genes predicted by a TF')+
    ylab('Pathway')+
    theme_minimal()+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
dev.off()
                     
    

write.csv(tf_percentages,
          'tf_promoter_prediction.csv')
                     
write.csv(rbindlist(tf_pathway_interactions),
          file='tf_gene_links.csv')
                     
                     