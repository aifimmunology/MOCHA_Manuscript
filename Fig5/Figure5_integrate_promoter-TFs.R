##############################################################################
##############################################################################

## Figure 5: analyses 
## 
##   - integrate chromVar
##   - with Promoter Accessibility
##   - via correlations

##############################################################################
##############################################################################
rm(list=ls())
require(data.table)
require(parallel)
require(stringi)
require(stringr)
require(ggplot2)
setwd('/home/jupyter/MOCHA_Manuscript/Fig5_integration')

## load data
pathway_hits = fread('covid_pathway_hits_annotated_v2.csv')
tsam = fread('TSAM_promoters.csv')
colnames(tsam)[1] = 'Sample'
tf_mat = read.csv('tf_matrix.csv')
colnames(tf_mat)[1] = 'Sample'
res_prom = fread('../Fig5/gene_hits.csv')
res_tf = fread('TF_quadratic_model.csv')

## 
predict_pathways <- function(pathway, plot_hits=F, tf='JUN'){
    
	p1 = strsplit(pathway_hits[description==pathway]$userId,
								';')[[1]]
	
	promoters = res_prom[intersect(grep(paste(p1, collapse = '|'), 
                                        res_prom$Gene),
                                   which(res_prom$cont_pval < 0.05)),]
	
	pathway_tsam = tsam[, promoters$Promoter, with=F]
    tf_vector = tf_mat[,tf]
	
	predict_tf = function(gene, plot=F, tf_vec=tf_vector){
		df =  data.frame(
			promoter= pathway_tsam[,gene,with=F][[1]], 
			TF =tf_vec
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
    res_p1$Unique_genes = length(p1)
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
pathway_list = pathway_hits[Level2=='Innate Immune System']

tfs_to_plot = c(top_positive$TF[1:5],
                top_negative$TF[1:5])
tf_pathway_interactions=lapply(tfs_to_plot,
       function(y)
           rbindlist(lapply(pathway_list$description,
			 function(x) 
			 	predict_pathways(x, plot_hits = F,tf=y
                                ))
                                )
                      
)


quantify_pathway <- function(tf_pathway_list, pathway){
    
    genes_tested = strsplit(pathway_hits[description==pathway]$userId,
								';')[[1]]
    
    tmp = tf_pathway_list[Pathway==pathway]
    tmp = tmp[tmp$Pval < 0.05]
    gene_hits = unique(unlist(strsplit(unique(tmp$Gene),', ')))
    sum(gene_hits %in% genes_tested) / length(genes_tested)
}

tf_percentages = lapply(1:length(tf_pathway_interactions),
               function(y) 
              data.frame(
                  pathway = pathway_list$description,
                  TF = rep(tfs_to_plot[y],18),
    Perc_gene = sapply(1:length(pathway_list$description),
       function(x)
           quantify_pathway(tf_pathway_interactions[[y]], pathway_list$description[x]
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
                     
                     
tf_percentages$Group = 'UpRegulated'
tf_percentages$Group[tf_percentages$TF %in% c('JUN',
                                              'FOSL1',
                                              'SMARCC1',
                                              'FOSL2',
                                              'JUNB')] = 'DownRegulated'

tf_percentages$pathway[tf_percentages$pathway=="Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways"] <- "Nucleotide-binding domain, leucine rich \nrepeat containing receptor (NLR) signaling pathways"
                     
pdf("tf_pathway_results3.pdf", width=9, height=8)
    ggplot(tf_percentages,
       aes(x=reorder(pathway, Perc_gene,mean),
           y=Perc_gene*100,
           col=Group,
           fill=Group))+geom_boxplot(scale='width')+
                                          theme_minimal()+
                     #ggpubr::stat_compare_means()+
                     theme(axis.text.x=element_text(size=12, angle=90),
                           text=element_text(size=14))+
                     ylab('% Significant Genes\nAssociated with TF')+
                     xlab('Pathway')+
                     geom_hline(yintercept=20)

                     
dev.off()
                     
### quantify TSAM zeroes 
zeroes = colMeans(tsam==0)         
summary(zeroes)

### unique genes                     
########################################################################               