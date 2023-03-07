### Load in enhancer-promoter links and convert to Hg38 (from Hg19)

library(plyranges)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)


allEnhancers <- read.table('ActivePromoterEnhancerLinks.tsv')

colnames(allEnhancers) <- gsub("\\(|\\)","", allEnhancers[1,])
allEnhancers <- as.data.frame(allEnhancers[-1,])

bait1 <- makeGRangesFromDataFrame(allEnhancers[,c(1:3)], seqnames.field = 'baitChr',
                                 start.field = 'baitSt', end.field = 'baitEnd')
cap1 <- makeGRangesFromDataFrame(allEnhancers[,c(5:7)], seqnames.field = 'oeChr',
                                 start.field = 'oeSt', end.field = 'oeEnd')
## Convert to Hg38

library(liftOver)
path = system.file(package="liftOver", "hg19ToHg38.over.chain")
ch = import.chain("hg19ToHg38.over.chain")

bait1_Hg38 = liftOver(bait1, ch)  %>% parallel::mclapply(., mergeAllRanges, mc.cores = 35)
length(bait1_Hg38)

cap1_Hg38 = liftOver(cap1, ch) %>% parallel::mclapply(., mergeAllRanges, mc.cores = 60)
length(cap1_Hg38)

all(names(cap1_Hg38) == names(bait1_Hg38))

Tile1 = MOCHA::GRangesToString(unlist(as(bait1_Hg38, 'GRangesList')))
Tile2 = MOCHA::GRangesToString(unlist(as(cap1_Hg38, 'GRangesList')))
databaseLinks <- data.frame(Tile1 = Tile1[-which(lengths(cap1_Hg38) == 0)],
                            Tile2 = Tile2,
                            cellTypes = allEnhancers$cellTypes[as.numeric(names(bait1_Hg38))][
                                -which(lengths(cap1_Hg38) == 0)])
write.csv(databaseLinks, 'ActivePromoterEnhancerLinks_Hg38.csv')
databaseLinks <- read.csv('ActivePromoterEnhancerLinks_Hg38.csv')


#######################################################################

#### function

#######################################################################

mergeAllRanges <- function(range){

    if(length(range) > 1){

        makeGRangesFromDataFrame(data.frame(seqnames = seqnames(range)[1],
                                        start = min(start(range)),
                                        end = max(end(range)),
                                        strand = strand(range)[1]))
    }else {range}

}

#' @title \code{varDecomp}
#'
#' @description \code{varDecomp} Runs PALMO to model variance decomposition along the variables
#'
#' @param Obj A RangedSummarizedExperment generated from getSampleTileMatrix
#' @param Donor_col The metadata column with donor information.
#' @param Time_col The metadata column with longitudinal information
#' @param Group_col The metadata column with group information. Default is NULL, in which case i
#' @param returnObj Boolean flag for whether or not you want the full PALMO object. Default is F
#'
#' @return variance decomposition matrix
#'
#'
#' @examples
#' \dontrun{
#'      vd_mat <- varDecomp(SampleTileObj, Donor_col = 'PTID', Time_col = 'Days')
#' )
#' }
#'
#' @export
#'
#'

varDecomp <- function(Obj, variableList, groupCol = NULL, percNonZero = 0, numCores = 1){

    if(length(SummarizedExperiment::assays(Obj)) > 1 ){

        variableList <- append(variableList, 'CellType')

        allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
        allMatrices[is.na(allMatrices)] = 0

        colnames(allMatrices) <- apply(expand.grid(names(SummarizedExperiment::assays(Obj)), uni
                                paste, collapse="__") %>% gsub(" ", "_", .)
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]
        meta <- parallel::mclapply(1:length(SummarizedExperiment::assays(Obj)), function(x){

                    SummarizedExperiment::colData(Obj) %>% as.data.frame() %>%
                        dplyr::mutate(Sample = paste(names(SummarizedExperiment::assays(Obj))[x]
                                    CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
                        dplyr::mutate(Sample = gsub(" ","_", Sample),
                                    CellType = gsub(" ", "_", CellType))

        }, mc.cores = numCores) %>% do.call('rbind', .)

    }else if(length(SummarizedExperiment::assays(Obj)) == 1 ){
            
        allMatrices <- SummarizedExperiment::assays(Obj)[[1]]
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]

        meta <- SummarizedExperiment::colData(Obj)

    }

    meta_f <- meta[meta$Sample %in% colnames(subMatrix),]

    #Filter out rows that are too sparse for modeling, and/or may cause issues. 
    if(is.null(groupCol)){

        subMatrix <- subMatrix[rowSums(subMatrix != 0)/ncol(subMatrix) > percNonZero,]

    }else if(groupCol %in% colnames(meta_f)){

        groupList <- unique(meta_f[,groupCol])
        keepRows <- mclapply(groupList, function(x){

            group_tmp <- meta_f$Sample[meta_f[,groupCol] == x]
            rowSums(subMatrix[,colnames(subMatrix) %in% group_tmp] != 0)/length(group_tmp) > per

        }, mc.cores = numCores) %>% do.call(cbind,.) %>% rowSums(.)

        keptTiles <- keepRows > 1

        subMatrix <- subMatrix[keptTiles,]

    }


    form <- paste(paste("(1|", variableList, ")", sep = ""), collapse = " + ")
    form1 <- as.formula(paste("exp ~ ", form, sep = ""))

    rowN <- dim(subMatrix)[1]
    op <- pbapply::pboptions(type = "timer")
                                   
    suppressMessages(lmem_res <- pbapply::pblapply(c(1:rowN),
        function(x) {
            df <- data.frame(exp = as.numeric(subMatrix[x, ]),
                meta_f, stringsAsFactors = FALSE)
            lmem <- lme4::lmer(formula = form1, data = df)

            lmem_re <- as.data.frame(VarCorr(lmem))
            row.names(lmem_re) <- lmem_re$grp

            lmem_re <- lmem_re[row.names(lmem_re), ]
            fix_effect <- fixef(lmem)  #get fixed effect
            lmem_re$CV <- lmem_re$sdcor/fix_effect

            normVar <- (lmem_re$vcov)/sum(lmem_re$vcov)*100
            names(normVar) <- row.names(lmem_re)

            output_df <- data.frame(Tiles = rownames(subMatrix)[x], Mean = mean(df$exp, na.rm =
                     Median = median(df$exp, na.rm = TRUE), SD = sd(df$exp, na.rm = TRUE),
                     Max = max(df$exp, na.rm = TRUE))
            cbind(output_df, t(as.data.frame(normVar)))

        }, cl = numCores), classes = "message")
    pbapply::pboptions(op)
    lmem_res <- do.call(rbind, lmem_res)
    return(lmem_res)

}
