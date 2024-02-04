
library(MOCHA)
library(SummarizedExperiment)
library(tidyverse)
library(plyranges)
library(lme4)
setwd('COVID_scATAC_Manuscript')

STM <- readRDS('CD16_Full_SampleTileMatrix.rds')

### Function for parallelized decomposition   
runDecomposition <- function(metaData, countMat, variableList, numCores = 1){

    form <- paste(paste("(1|", variableList, ")", sep = ""), collapse = " + ")
    form1 <- as.formula(paste("exp ~ ", form, sep = ""))

    rowN <- dim(countMat)[1]
    op <- pbapply::pboptions(type = "timer")

    suppressMessages(lmem_res <- pbapply::pblapply(c(1:rowN),
        function(x) {
            df <- data.frame(exp = as.numeric(countMat[x, ]),
                metaData, stringsAsFactors = FALSE)
            lmem <- lme4::lmer(formula = form1, data = df)

            lmem_re <- as.data.frame(VarCorr(lmem))
            row.names(lmem_re) <- lmem_re$grp

            lmem_re <- lmem_re[row.names(lmem_re), ]
            fix_effect <- fixef(lmem)  #get fixed effect
            lmem_re$CV <- lmem_re$sdcor/fix_effect

            normVar <- (lmem_re$vcov)/sum(lmem_re$vcov)
            names(normVar) <- row.names(lmem_re)

            output_df <- data.frame(Tiles = rownames(countMat)[x], 
                                    Mean = mean(df$exp, na.rm = TRUE),
                                    SD = sd(df$exp, na.rm = TRUE))
            cbind(output_df, t(as.data.frame(normVar)))

        }, cl = numCores), classes = "message")
    pbapply::pboptions(op)
    lmem_res <- do.call(rbind, lmem_res)
    return(lmem_res)
}
           
i = 'CD16 Mono'
accMat <- getCellPopMatrix(STM[,!is.na(STM$days_since_symptoms)], i, NAtoZero = TRUE)
accMat2 <- log2(accMat+1)

metaData <- STM[,!is.na(STM$days_since_symptoms)]@colData

## Look for regions that are 80% non-zeros within at least one infection stage
listReg <- group_by(as.data.frame(metaData), InfectionStages) %>% summarize(SampleList = list(Sample)) 
strongPeaks <- unique(unlist(lapply(listReg$SampleList, function(x){
    
    which(apply(accMat2[,x], 1, function(x){ sum(x !=0)/length(x) >= 0.8 }))
    
})))

deComp <- runDecomposition(metaData, accMat2[strongPeaks,], 
                    variableList = c('Age', 'Sex', 'days_since_symptoms','PTID'),
                               numCores = 60)
write.csv(deComp, 'CD16_Mono_VarianceDecomposition.csv')

topVar <- colnames(deComp)[c(4:8)][unlist(apply(deComp, 1, function(XX) which.max(XX[c(4:8)])))]
                                                
varTable <-  as.data.frame(table(topVar))
                                                   
varDF <- t(data.frame(TileNumber = varTable$Freq,
           PercentOfTiles = varTable$Freq/sum(varTable$Freq)))
colnames(varDF) <- as.data.frame(table(topVar))$topVar
                                                   
write.csv(varDF, 'SupplementalTable3_VarianceDecomposition.csv')
                                                          
                                                   