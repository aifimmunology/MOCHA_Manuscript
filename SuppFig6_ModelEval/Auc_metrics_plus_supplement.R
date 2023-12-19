# ###########################################################
# ###########################################################

# Author: Samir Rachid Zaim
# Script: Fig. 1 D Downsampling
#         evaluation matrics 

# ###########################################################
# ###########################################################

require(TxDb.Hsapiens.UCSC.hg38.refGene)
require(cutpointr)
require(data.table)
require(ggplot2)
require(MOCHA)
library(ArchR)
require(parallel)
library(GenomicRanges)
library(plyranges)
require(GGally)
source('/home/jupyter/theme.R')
source('/home/jupyter/MOCHA_Manuscript/Fig2/helper_granges.R')
#require(arules)
require(precrec)
require(PRROC)

# ###########################################################
# ###########################################################

## load ArchR Project 
ArchRProj <- ArchR::loadArchRProject('/home/jupyter/FullCovid/')

## Extract Metadata and 
## required params for 
## calling open tiles 
metadf <- getCellColData(ArchRProj) 
studySignal = median(metadf$nFrags)
celltypes <- c('B naive','cDC2',
               'NK_CD56bright', 'CD14 Mono','NK')
############################################################
############################################################

### load scMACS peaks
### for plotting and 
### analyses 
setwd('/home/jupyter/covid/scMACS_manuscript_analyses/data')

getROC<- function(cell){
  
    ############################################################
    # Get our fragments for this cellPop
    frags <- MOCHA:::getPopFrags(
          ArchRProj = ArchRProj,
          metaColumn = 'predictedGroup_Co2',
          cellSubsets = cell,
          region = NULL,
          numCores = 30,
          sampleSpecific = FALSE,
          NormMethod = "nfrags",
          blackList = NULL,
          verbose = FALSE,
          overlapList = 50
      )
    ############################################################    
    
    ############################################################     
    ### get macs2 groundtruth 
    numCells = length(unique(frags[[1]]$RG))
    
    if(cell=='NK'){
        macs2_truth=fread('/home/jupyter/covid/scMACS_manuscript_analyses/trainModel/NK_trainingset.csv')
        
    } else{
        cell2 = gsub(' |_','.', cell)
        homeDir = '/home/jupyter/covid/scMACS_manuscript_analyses/data/'
        fname = paste(homeDir,cell,'/bulk/downsample_macs2_peaks/', 
                    max(numCells),'_cells_sample_1_',
                      cell2,
                      '_peaks.broadPeak',sep='')

        macs2_groundtruth <- fread(fname)

        macs2_truth <- peak_to_tile_macs2_bulk(1, list(macs2_groundtruth))
    }
    ############################################################    
    
    subsampleCells <- function(n){
            unique_cells <- unique(frags[[1]]$RG)
            total_cells = length(unique_cells)
            sampleCells <- sample(unique_cells, min(n,total_cells))
            sampleCells
    }

    get_samplePredictions <- function(n){
        
        cells <- subsampleCells(n)
        cellsSample <- which(ArchRProj$cellNames %in% cells)
        subArchRProj <- ArchRProj[cellsSample, ]

        # Parameters for calling open tiles.
        cellPopLabel <- "CellSubsets" 
        cellPopulations <- c(cell)
        numCores <- 40

        # ###################################################
        # 2. Call open tiles (main peak calling step)
        #    and get sample-tile matrices
        #    for all specified cell populations
        # ###################################################


        tmpfrags <-frags[[1]]
        tmpfrags =  tmpfrags[tmpfrags$RG %in% cells]
        
        study_prefactor <- 3668 / studySignal
        normalization_factors <- length(tmpfrags)
        
        ## 
        library(TxDb.Hsapiens.UCSC.hg38.refGene)
        library(org.Hs.eg.db)
        TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
        Org <- org.Hs.eg.db

        tilesGRangesList <- MOCHA:::callTilesBySample(
              blackList = ArchR::getBlacklist(subArchRProj),
              returnAllTiles = TRUE,
              totalFrags = normalization_factors,
              fragsList = tmpfrags,
              verbose = TRUE,
              StudypreFactor = study_prefactor
            )
        
        preds <- as.data.frame(tilesGRangesList)
        preds$Macs2 = ifelse(preds$tileID %in% macs2_truth$tileID, 1,0)
        preds = as.data.table(preds)
        preds$Ncells = n
        
        TP <- sum(preds$peak & preds$Macs2==1)
        FP <- sum(preds$peak & preds$Macs2==0)
        FN <- sum(!preds$peak & preds$Macs2==1)
        TN <- sum(!preds$peak & preds$Macs2==0)
        
        sensitivity = TP / (TP+FN)
        Specificity = TN/(TN+FP)
        PPV = TP/(TP+FP)

        ### calculate the optimal 
        ### threshold based on maximized
        ### youden index
        youden_thresholds <- cutpointr(data=preds, class=Macs2, 
                                     Prediction, 
                                     direction = ">=", 
                           pos_class = 1,
                           neg_class = 0,
                           method = cutpointr::maximize_metric,
                           metric = cutpointr::youden
                                    )

        youden_thresholds = as.data.table(youden_thresholds)
        youden_thresholds$Ncells <- n
        youden_thresholds$PPV = PPV
        youden_thresholds$sensitivity =sensitivity
        youden_thresholds$Specificity = Specificity
        
        pryr::mem_used()
        
        return(list(preds,
                   youden_thresholds))
    }
    
    
    ### Select Downsamples
    cellQuants <- unique(c(5,20, 25, 50, 100,200, 500,750, 1000, 2000, 5000,
                          min(numCells, 10000),
                          min(numCells, 15000),                           
                          min(numCells, 25000),
                          min(numCells, 50000),
                          min(numCells, 100000),
                          min(numCells,135000))
                         )
                    
    ### Parallelize Downsampling 
    preds <- parallel::mclapply(cellQuants,
           function(x)
                    get_samplePredictions(x),
                      mc.cores=20,
                                mc.preschedule=F
           )
  
    ## Extract matrices 
    ## with metrics
    youden <- rbindlist(lapply(preds,
                     function(x) 
                     x[[2]]))
    
    
    mat_list <- lapply(preds,
                       function(x)
                           x[[1]]$Prediction
                       )
    labels <- lapply(preds,
                       function(x)
                           x[[1]]$Macs2
                       )    
    modnames = sapply(preds,
                      function(x) x[[1]]$Ncells[1]
                      )
                               
    modnames.txt <- as.character(modnames)
                       
    mdat <- mmdata(mat_list,
               labels, 
               modnames = modnames.txt,
               dsids = modnames)
                               
    setwd('/home/jupyter/MOCHA_Manuscript/SuppFig2_ModelEval/')

    #############################################################################
    #############################################################################   
    if(cell=='NK'){
        
       idx = which(cellQuants %in% c(500, 2000, 5000, 10000, 25000))
        
       subsetted_mat_list = mclapply(idx, 
               function(x){
                        df = data.table(
                            pred = mat_list[[x]],
                            Macs2_label= as.factor(labels[[x]]),
                            Cells = cellQuants[x]
                            )
                        df = df[sample(nrow(df), min(nrow(df), 100000)),]
               },
                                     mc.cores=5
               )
        subsetted_mat = rbindlist(subsetted_mat_list)
        
        pdf('test.pdf', width=)
        p= ggplot(subsetted_mat,
               aes(x=pred,
                   fill=Macs2_label))+geom_histogram(binwidth=0.05)+theme_minimal()+
                   facet_wrap(~Cells, scale='free_y', ncol=2)+
                    theme(text=element_text(size=12),
                         legend.position='none')+
                          xlab("MOCHA Prediction") + ylab('Distribution')
                   
        print(p)
        dev.off()
        
        write.csv(subsetted_mat,
                  file='NK_probabilities.csv')
                               
    #############################################################################
    #############################################################################   
                               
                               
    setwd('/home/jupyter/MOCHA_Manuscript/Fig1/')
                        
    ## Generate an mscurve object that contains ROC and Precision-Recall curves
    mscurves <- evalmod(mdat)
    ## Shows AUCs
    png(paste('auc_plots/',cell,'.png',sep=''),
        width=400,
        height=400
   )
    plot(mscurves, curvetype='ROC', show_legend=F)
    dev.off()
    
    aucs <- lapply(preds,
         function(x){
             x =x[[2]]
             x$AUC
             }
         )
    roc_df <- data.frame(AUC = unlist(aucs),
                     Cells= cellQuants)
       
    #### Generate Plots:
    
    png(paste('auc_plots/',cell,'_aucs.png',sep=''),
        width=400,
        height=400
   )     

    p <- ggplot(roc_df,
           aes(x=Cells,
               y=AUC))+geom_point()+geom_smooth(span=0.5, se=FALSE)+ylim(c(0,1))+
            theme(legend.title=element_text(size=16),
                legend.text=element_text(size=16),
                    axis.title.y = element_text(size=20),
                    axis.title.x = element_text(size=20),
                    strip.text.x = element_text(size=20),
                    axis.text.x = element_text(size=17, angle=90),
                    axis.text.y = element_text(size=17))+
                  xlab('Cells Sampled') + ylab('AUC')+
                               theme_minimal()
    print(p)
    dev.off()
    
    ### ALl Metrics                               
    png(paste('auc_plots/',cell,'_metrics.png',sep=''),
        width=400,
        height=400
    )     

    youden <- melt(youden, measure.vars=c('Specificity',
                                          'sensitivity',
                                          'AUC',
                                          'youden'))
    p <- ggplot(youden,
           aes(x=Ncells,
               y=value,
               col=variable))+geom_point()+
                               geom_line(linewidth=2)+
                               ylim(c(0,1))+
                               #facet_wrap(~variable)+
            theme( axis.title.y = element_text(size=20),
                    axis.title.x = element_text(size=20),
                    strip.text.x = element_text(size=20),
                    axis.text.x = element_text(size=17, angle=90),
                    axis.text.y = element_text(size=17))+
                  xlab('Cells Sampled') + ylab('Value')+
                               theme_minimal()+
                               theme(text=element_text(size=20),
                                     legend.position='none',
                                    axis.text.x=element_text(size=20, angle=90))
    print(p)
    dev.off()
              
    #### CD16 is main Plot 
    if(cell=='CD16 Mono'){
        
    pdf(paste('plots/Fig1D',cell,'_metrics.pdf',sep=''),
        width=4,
        height=6
    )     

    p <- ggplot(youden,
           aes(x=Ncells,
               y=value,
               col=variable))+geom_point()+
                               geom_line(linewidth=2)+
                               ylim(c(0,1))+
                               #facet_wrap(~variable)+
            theme( axis.title.y = element_text(size=20),
                    axis.title.x = element_text(size=20),
                    strip.text.x = element_text(size=20),
                    axis.text.x = element_text(size=17, angle=90),
                    axis.text.y = element_text(size=17))+
                  xlab('Cells Sampled') + ylab('Value')+
                               theme_minimal()+
                               theme(text=element_text(size=20),
                                     legend.position='none',
                                    axis.text.x=element_text(size=20, 
                                                             angle=90))
    print(p)
    dev.off()
        
    }
    #############################################################################       
    final_results = youden[, c('Ncells','variable','value')]
    write.csv(final_results,
               file= paste('supplementalFiles/',cell,'_metrics.csv',sep=''))

    rm(frags, aucs)
    pryr::mem_used()
    gc()

}


lapply(celltypes, 
       function(x) getROC(x)
       )

# ########################################################################
# ########################################################################

library(pROC)
require(PRROC)
require(ROCR)
                               
                               
probs_subset <- scMACS_peaks2[seq(1, 83, by=3)]
probs_subset <- rbindlist(probs_subset)
probs_subset <- probs_subset[,c('tileID','Prediction','numCells')]


macs2_groundtruth <- fread('/home/jupyter/covid/scMACS_manuscript_analyses/data/CD14 Mono/bulk/downsample_macs2_peaks/135949_cells_sample_1_CD14.Mono_peaks.broadPeak')

macs2_truth <- peak_to_tile_macs2_bulk(1, list(macs2_groundtruth))
probs_subset$Peak <- ifelse(probs_subset$tileID %in% macs2_truth$tileID,
                            1,0)


numCells = unique(probs_subset$numCells)
aucs <- mclapply(numCells,
         function(x){
         mat= probs_subset[numCells==x]
         auc =roc.curve(scores.class0 = mat$Prediction[mat$Peak==1],
               scores.class1 = mat$Prediction[mat$Peak==0])$auc
             },
         mc.cores=15
         )
roc_df <- data.frame(AUC = unlist(aucs),
                     Cells=numCells)

setwd('/home/jupyter/covid/all_methods')
png('aucs_cd14.png',
   width=400,
   height=400)

ggplot(roc_df,
       aes(x=Cells,
           y=AUC))+geom_point()+geom_smooth(se=FALSE)+ylim(c(0,1))+
        theme(legend.title=element_text(size=16),
            legend.text=element_text(size=16),
                axis.title.y = element_text(size=20),
                axis.title.x = element_text(size=20),
                strip.text.x = element_text(size=20),
                axis.text.x = element_text(size=17, angle=90),
                axis.text.y = element_text(size=17))+
              xlab('Cells Sampled') + ylab('AUC')

dev.off()



aucs <- mclapply(numCells,
         function(x){
         mat= probs_subset[numCells==x]
         auc =roc.curve(scores.class0 = mat$Prediction[mat$Peak==1],
               scores.class1 = mat$Prediction[mat$Peak==0],curve=TRUE)
             },
         mc.cores=15
         )



mat_list <-  mclapply(numCells,
         function(x){
         mat= probs_subset[numCells==x]
         preds = mat$Prediction
             },
         mc.cores=15
         )

labels <-  mclapply(numCells,
         function(x){
         mat= probs_subset[numCells==x]
         preds = mat$Peak
             },
         mc.cores=15
         )


require(precrec)

idx = seq(3, length(numCells), by=3)

mdat <- mmdata(mat_list[idx],
               labels[idx], 
               modnames = as.character(numCells)[idx],
               dsids = numCells[idx])

## Generate an mscurve object that contains ROC and Precision-Recall curves
mscurves <- evalmod(mdat)
## Shows AUCs
png('/home/jupyter/covid/all_methods/cd14_roc.png',
    width=400,
    height=400
   )
plot(mscurves, curvetype='ROC', show_legend=F)
dev.off()






nk_mat = fread('/home/jupyter/MOCHA_Manuscript2/SuppFig2_ModelEval/NK_metrics.csv')                           
pdf('NK_performance.pdf')
        p= ggplot(nk_mat,
               aes(x=Ncells,
                   y=value,
                   col=variable))+ geom_line(lwd=3)+theme_minimal()
                   
        print(p)
        dev.off()