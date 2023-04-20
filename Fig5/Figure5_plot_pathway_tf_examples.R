# #############################################################################
# #############################################################################

# # Figure 5: analyses 
# # 
# #   - integrate chromVar
# #   - with Promoter Accessibility
# #   - via correlations

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
res_prom = fread('../gene_hits.csv')
res_tf = fread('TF_quadratic_model.csv')

## 
plot_pathways <- function(pathway, plot_hits=F, tf='JUN'){
    
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
							 		y=promoter))+geom_point(size=0.5,
                                                           alpha=0.6)+
				geom_line(data=df,
									aes(x=TF,
											y=prediction),
									linewidth=1.5)+
				theme_minimal()+
				ylab(paste(promoters$Gene[gene],
                           promoters$Promoter[gene],
                           sep='\n'))+
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
	
	p=ggpubr::ggarrange(plotlist=lapply(which(res_p1$`Pr(>|z|)`<0.05),
					function(x) predict_tf(x, plot=T
                                                                                                  )
                                            
                                            )
                                            )
    p= ggpubr::annotate_figure(p,
                            top=ggpubr::text_grob(pathway,
                                          color='black'))
    return(p)    
    
}
#####
require(dplyr)                    
top_tfs = res_tf[time_pval < 0.05,]
top_positive = top_tfs[time>0] %>% arrange(-time)
top_negative = top_tfs[time<0] %>% arrange(time)

## expand list to entire innate immune pathways 
pathway_list = pathway_hits[Level2=='Innate Immune System']

tfs_to_plot = c(top_positive$TF[1:5],
                top_negative$TF[1:5])

pdf('positive_tfs.plots.pdf')
positive_tf_plots=lapply(top_positive$TF[1:2],
       function(y)
           lapply(pathway_list$description,
			 function(x) 
			 	plot_pathways(x, plot_hits = T,tf=y
                                )
                            )
                                )
                        
print(positive_tf_plots)
dev.off() 


pdf('JUNB_tfs.plots.pdf',
   height=10, width=10)

## Show Toll like Receptor 4 pathway
a = plot_pathways(pathway_list$description[2], plot_hits = T,tf=top_negative$TF[1])
print(a)
dev.off()
                                
           
                
                     

write.csv(rbindlist(tf_pathway_interactions),
          file='tf_gene_links.csv')

########################################################################               
