###############################################################################
#install.packages(c('lme4','glmmTMB' 'ggbreak','lmerTest','WebGestaltR'))
require(glmmTMB)
require(data.table)
require(MOCHA)
require(ggplot2)
require(parallel)

### laod data 
stm = readRDS('/home/jupyter/MOCHA_Manuscript2/Fig5/SampleTileObject-2.rds')

### extract CD16 Promoters 
tmp <- plyranges::filter(SummarizedExperiment::rowRanges(stm), tileType == 'Promoter')


summarize_dropout <- function(cellType){
    
    tmpDF = as.data.table(tmp)
    colnames(tmpDF) = gsub('\\.', ' ', colnames(tmpDF))
    tmpDF$Feature = tmpDF[, cellType,with=F]
    
    promoters = tmp[which(tmpDF$Feature==T),]

    tsam = MOCHA::getCellPopMatrix(stm, cellPopulation = cellType)
    promoter_tsam <- tsam[row.names(tsam) %in% MOCHA::GRangesToString(promoters),]

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

    promoter_tsam = promoter_tsam[, colnames(tsam) %in% covid_samples]
    promoter_tsam = data.frame(t(promoter_tsam)) 
    promoter_tsam$Sample = row.names(promoter_tsam)


    meta <- meta[Sample %in% promoter_tsam$Sample, c('Sex','time','time_sqrd','Age','PTID','Sample')]
    cellCounts=cellCounts[cellTypeLabelList==cellType]


    metadf = meta
    promot_df = promoter_tsam

    cellCounts = cellCounts[cellTypeLabelList==cellType & V2 %in% covid_samples,]
    cellCounts$Sample = cellCounts$V2

    ###############################################################################
    ###############################################################################

    # ##############################################################################
    # ##############################################################################

    fitTechnicalNoise <- function(y,
                                 index){

        X = metadf
        X$y = y
        X= dplyr::left_join(X, cellCounts)	
        X$CellCounts = X$N
        X$Sex = ifelse(X$Sex=='Female',1,0)

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

            res = data.frame(
                Estimate = rep(1,4),
                sigma=rep(1,4),
                Z = rep(1,4),
                Pval=rep(1,4))





        } else{
            fit = try(glmmTMB(y ~  Age + Sex +
                              time + (1|PTID) ,                     
                                    ziformula = ~0+CellCounts+Age+Sex,
                                    data=X,
                                    family = gaussian(),
                                    REML=T))

            a = summary(fit)

            res = data.frame(a$coefficients$zi)
            res$Index=index


            colnames(res) = c('Estimate',
                              'Std. Error',
                               'z value',
                              'Pr(>|z|)',
                             'Index')    

        return(res)
            }

    }

    # ##############################################################################
    # ##############################################################################

    ## test zeroes on cell counts 
    setwd('../MOCHA_Revision')

    cl <- makeCluster(60)
    clusterExport(cl, c("fitTechnicalNoise", "promot_df","cellCounts",
                        "glmmTMB","metadf","glmmTMBControl"),
                          envir=environment())

    start = Sys.time()
    technical_zeroes=parLapply(1:(ncol(promot_df)-1),
                 function(x)
                    try(fitTechnicalNoise(promot_df[,x], index=x)),
                          cl=cl)
    Sys.time()-start

    stopCluster(cl)


    ##### remove peaks
    ##### that had convergence errors 
    technical_zeroes = technical_zeroes[sapply(technical_zeroes, class)!='try-error']

    ##### remove peaks that did not 
    ##### have zeroes 
    technical_zeroes = technical_zeroes[sapply(technical_zeroes, class)!='character']

    technical_zeroes = rbindlist(technical_zeroes, fill=TRUE)
    technical_zeroes = technical_zeroes[ is.na(sigma)]

    technical_zeroes = technical_zeroes[, c("Estimate",
                                        "Std. Error",
                                        "z value",
                                        "Pr(>|z|)",
                                        "Index")]
    
    colnames(technical_zeroes)[4] = 'Pvalue'

    technical_zeroes$Variable = rep(c('Cellcounts','Age','Sex'),
                                    nrow(technical_zeroes)/3)
    

    pdf(paste(cellType,'technical0s.pdf', sep='_'))
    p1= ggplot(technical_zeroes,
           aes(x=Pvalue))+geom_histogram(bins=100)+facet_wrap(~Variable, 
                                                      scales='free',
                                                      nrow=3)+
        theme_minimal()+theme(text=element_text(size=12))+
        ggtitle(paste('Testing for 0-inflation in',cellType))
    print(p1)
    dev.off()


    rowMeans(promot_df==0)

    zi_summary = data.frame(
        Sample = row.names(promot_df),
        ZI =rowMeans(promot_df==0))

    zi_summary = dplyr::left_join(zi_summary, cellCounts[cellTypeLabelList==cellType,
                                                  c('N','Sample')],
                           by='Sample')

    zi_summary = dplyr::left_join(zi_summary,
                                  metadf)
    pdf(paste(cellType, 'zeros_vs_cellCounts.pdf',sep='_'))

    p = ggplot(zi_summary,
           aes(x=N,
               y=ZI))+geom_point()+geom_smooth() +theme_minimal()+
    theme(text=element_text(size=16)) + 
     scale_x_log10()+

    xlab('Log 10 Cell Counts')+ylab('Proportion of Tiles with 0s')+
    ggtitle('Proportion of Promoter Tiles with 0s by Cell Counts')

    print(p)
    dev.off()

    plot(x=zi_summary$Age,
         y=jitter(promot_df[,3]))
    dev.off()
}


##################################################################
##################################################################
summarize_dropout("CD16 Mono")
summarize_dropout("B Naive")
summarize_dropout("CD4 Naive")
