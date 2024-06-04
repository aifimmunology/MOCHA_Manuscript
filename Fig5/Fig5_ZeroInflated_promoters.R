##############################################################################
##############################################################################
## GLMMTMBs
##
##
##############################################################################

###############################################################################
###############################################################################
#install.packages(c('lme4','glmmTMB', 'ggbreak','lmerTest','WebGestaltR'))
require(glmmTMB)
require(data.table)
require(MOCHA)
require(ggplot2)
require(parallel)
require(lmerTest)
require(scattermore)
### laod data 
stm = readRDS('/home/jupyter/MOCHA_Manuscript2/Fig5/SampleTileObject-2.rds')

### extract CD16 Promoters 
tmp <- plyranges::filter(SummarizedExperiment::rowRanges(stm), tileType == 'Promoter')
cd16_promoters = tmp[tmp$`CD16 Mono`==T]
cd16_tsam = MOCHA::getCellPopMatrix(stm, cellPopulation = 'CD16 Mono')
cd16_promoter_tsam <- cd16_tsam[row.names(cd16_tsam) %in% MOCHA::GRangesToString(cd16_promoters),]

### Extract Metadata 
meta = as.data.table(stm@colData)
cellCounts = data.table(stm@metadata$CellCounts)

###############################################################################
###############################################################################

### 
meta$time <- meta$days_since_symptoms-mean(meta$days_since_symptoms, na.rm=T)
meta$time_sqrd <- (meta$time)^2
numCores = 1

### Establish group1 
### and group 2
Covid_subjects = c(32209, 31207, 32140, 32245, 32054, 
                   32255, 32131, 31945, 32416, 42409, 
                   32220, 32038, 31874, 32124, 31924, 
                   32251, 32196, 32415)   

### Extract samples
### from group 1 and 2
covid_samples <- meta$Sample[meta$PTID %in% Covid_subjects]

# ## create separate 
# ## TSAMs for each 

tsam = cd16_promoter_tsam[, colnames(cd16_promoter_tsam) %in% covid_samples]
tsam = as.data.frame(t(tsam)); tsam$Sample = row.names(tsam)


meta <- meta[Sample %in% tsam$Sample, c('Sex','time','time_sqrd','Age','PTID','Sample')]
cellCounts=cellCounts[cellTypeLabelList=='CD16 Mono']


metadf = meta
promot_df = tsam

cellCounts = cellCounts[cellTypeLabelList=='CD16 Mono' & V2 %in% covid_samples,]
cellCounts$Sample = cellCounts$V2

###############################################################################
###############################################################################

# ##############################################################################
# ##############################################################################

fit_glmmTMB <- function(y){

	X = metadf
    X$y = y
    X= dplyr::left_join(X, cellCounts)	
    X$CellCounts = X$N
    ### Zero-Inflated LMM
    ### Gaussian Identity 
    ### Provides Linear
    ### Regression as 
    ### 'identity link function'
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
                                ziformula = ~ 0+CellCounts,
                                data=X,
                                family = gaussian(),
                                REML=T))
        }
    

	fit

}

# ##############################################################################
# ##############################################################################

# ## train LMMs for all covid subjects

cl <- makeCluster(60)
clusterExport(cl, c("fit_glmmTMB", "promot_df","cellCounts",
                    "glmmTMB","metadf","glmmTMBControl"),
                      envir=environment())

start = Sys.time()
models=parLapply(1:49679,
			 function(x)
			 	fit_glmmTMB(promot_df[,x]),
                      cl=cl)
Sys.time()-start

stopCluster(cl)

###############################################################################
###############################################################################
setwd('/home/jupyter/MOCHA_Manuscript/Fig5/')
######## Analyze purely continuous shift 

# ## Quantify # of Promoter hits 
# ## Per group 

extract_coefs <-  function(x,group1_models){
    
    
    print(x)
    
    fit = group1_models[[x]]
    
    if(class(fit)=='glmmTMB'){
            summary_fit = summary(fit)
           
            res = data.frame(
                cont_pval = summary_fit$coefficients$cond['time',4],
                cont_coef = summary_fit$coefficients$cond['time',1])

         
        } else{ 
          res = data.frame(
                cont_pval = NA,
                cont_coef = NA)
        }
       return(res)
        
}

summarize_results <- function(models, fname){
   
    res_list_group1 = mclapply(1:length(models),
                                 function(x)
                                     try(extract_coefs(x,models)),
                               mc.cores=5
                             )
    
    failed_to_converge <- which((sapply(res_list_group1, class))=='try-error')
    res_list_group2 = res_list_group1
    
    tmp_res <-  data.frame(
                cont_pval = NA,
                cont_coef = NA)
    
    res_list_group2[failed_to_converge] =lapply(failed_to_converge,
           function(x)
               res_list_group2[[x]] <- tmp_res
           )
    
    res_group1 = rbindlist(res_list_group2)
    res_group1$Gene = cd16_promoters$Gene
    res_group1$Promoter = row.names(cd16_promoter_tsam)

    return(res_group1)
}

### plot results 
res_group1 = summarize_results(models)
res_group1 = res_group1 %>% arrange(cont_pval)
res_group1$FDR = p.adjust(res_group1$cont_pval, 'fdr')
### Extract Interesting Genes 
unique_genes <- unique(res_group1$Gene[res_group1$FDR <0.1])

### some promoters span
### multiple genes 
### so we extract them 
group1_genes = sapply(unique_genes,
                      function(x)
                          unlist(strsplit(x, ', '))
                      )
unlist(group1_genes) -> genes
###############################################################################
###############################################################################
#install.packages('WebGestaltR')
require(WebGestaltR)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
require(ggbreak)
require(dplyr)
### Confirm Pathway hits 
enrichDataBaseList <- WebGestaltR::listGeneSet()  

simplifiedORA <- function(database, foreground, background){
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                             fdrThr=0.1,
                            referenceGeneType= "genesymbol")
}

allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE))
refList <- mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID")

group1_reactome_pathways <- simplifiedORA(enrichDataBaseList$name[9],
                               genes, refList)

pdf('all_subjects-reactome.pdf', width=19, height=19)
    ggplot(group1_reactome_pathways,
       aes(y=reorder(description, FDR, mean),
           x=-log10(FDR)))+geom_point()+
        theme_minimal()+ylab('Reactome Pathway')+
        theme(text=element_text(size=22))
dev.off()

res_group1 = res_group1 %>% arrange(cont_pval)

## assign 0 pvalues
## to smallest value detected
## in model
#res_group1$cont_pval[res_group1$cont_pval==0] <- 6.684058e-267

pdf('all_subjects-Volcano.pdf', width=9, height=6)
p1 =ggplot(res_group1,
       aes(x=cont_coef,
           label=Gene,
           y=-log10(FDR)))+geom_scattermore()+
        theme_minimal()+ylab('-log(FDR)')+
        xlab('Slope')+
        theme(text=element_text(size=18))+
 ggrepel::geom_text_repel(
    data = res_group1[FDR < 0.001,],
    aes(label = Gene),
    size = 3,
    point.padding = unit(0.3, "lines")
  )
print(p1)
dev.off()


# ####### Analyze purely binary shift 

write.csv(group1_reactome_pathways,
          file='covid_pathway_hits.csv')

write.csv(res_group1,
          file='gene_hits.csv')



## Plot examples
plot_examples <- function(i){
	X = metadf
    X$y = promot_df[,i]
    X= dplyr::left_join(X, cellCounts)	
    X$CellCounts = X$N
    X$PTID = factor(X$PTID)
    
    X = X[PTID %in% c(32124, 32140, 31924, 32038, 32209)]
    X$time =X$time +39
    p1=ggplot(X,
           aes(x=time,
               y=y,
               col=PTID))+
    geom_point()+theme_minimal()+
    theme(legend.position='none')+
    facet_wrap(~PTID, ncol=1)+
    geom_smooth(method='lm', se=F)+
    ylim(0,NA)
    
    X2 = X
    X2 = X2[X2$y>0,] 
    p2=ggplot(X2,
           aes(x=time,
               y=y,
               col=PTID))+geom_smooth(method='lm', se=F)+
    geom_point()+theme_minimal()+facet_wrap(~PTID, ncol=1)+
    theme(legend.position='none')+
        ylim(0,NA)
    ggpubr::ggarrange(p1, p2, ncol=2)
    
}

pdf('linePlot_examples.pdf')
lapply(rpois(20, lambda=45000),
      function(x)
          try(plot_examples(x)))
dev.off()
       
       
       

       
### Extract model for 
### a subset of genes 
### of interest 
require(dplyr)
genes = c('HIVEP1','IKBK','NFKB1','NFKBIE','RELA','RELB')
interesting_genes = res_group1[Gene %in% genes]    
interesting_genes = interesting_genes %>% arrange(cont_pval)
promoter_hits = interesting_genes[, first(.SD),by=Gene, .SDcols=c(1:4)]      
promoter_hits=promoter_hits[,-1]       

coefs =lapply(1:5,
       function(x){
           
           prom_idx=  which(colnames(promot_df) %in% promoter_hits$Promoter[x] )
           prom_t=promoter_hits$Promoter[x]
           res = data.frame(t(summary(models[prom_idx][[1]])$coefficients$cond[,1]))
           res$Promoter =prom_t
            res
           }
       )
coefs = data.frame(do.call(rbind, coefs) )   
time_grid = metadf$time
       
predict_matrix <- 
       lapply(1:5,
              function(x) 
                  as.numeric(coefs[x,1])+as.numeric(coefs[x,3])*mean(X$Age)+
       as.numeric(coefs[x,5])*time_grid
              )
TF_family =data.table(do.call(cbind, predict_matrix))
colnames(TF_family) <- promoter_hits$Gene
TF_family$Time = time_grid +39      

png('tf_family.png')
ggplot(
    melt(TF_family,id.var='Time'),
    aes(x=Time,
        y=value,
        col=variable))+geom_point()+geom_line()+
       theme_minimal()+
       theme(text=element_text(size=16))
dev.off()
       
       
       
       
######### Extract Full Model Coefficients 
# Per group 

extract_full_coefs <-  function(x,models){
            print(x)
            fit = models[[x]]
            summary_fit = summary(fit)
           
            res = as.data.frame(t(summary_fit$coefficients$cond[,1]))
            res2=as.data.frame(t(summary_fit$coefficients$cond[,4]))
            colnames(res) = paste(colnames(res)[1:4],'coef', sep='_')
            colnames(res2) = paste(colnames(res2)[1:4],'pval', sep='_')
            mat = cbind(res, res2)
       
            if(is.null(summary_fit$coefficients$zi)){
                  
               mat$ZI_Estimate = NA
               mat$ZI_pval= NA
            } else {
               mat = cbind(mat, t(summary_fit$coefficients$zi[,c(1,4)]))
               colnames(mat)[9:10] <- c('ZI_Estimate','ZI_pval')
            }
            
            return(mat)

}

summarize_results2 <- function(models, fname){

    res_list_group1 = mclapply(1:length(models),
                                 function(x)
                                     try(extract_full_coefs(x,models)),
                               mc.cores=5
                             )
    idx = which(sapply(res_list_group1, class)!='try-error')
    res_group1 = rbindlist(res_list_group1[idx])
    res_group1$Gene = cd16_promoters$Gene[idx]
    res_group1$Promoter = row.names(cd16_promoter_tsam)[idx]

    return(res_group1)
}

### save results 
    
full_model_coefs = summarize_results2(models)  
full_model_coefs = full_model_coefs %>% arrange(time_pval)
full_model_coefs$FDR = p.adjust(full_model_coefs$time_pval, 'fdr')

write.csv(full_model_coefs,
          file='zero_inflated_coefficients.csv')