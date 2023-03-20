#' @title \code{varDecomp}
#'
#' @description \code{varDecomp} Runs PALMO to model variance decomposition along the variables needed. 
#'
#' @param Obj A RangedSummarizedExperment generated from getSampleTileMatrix
#' @param Donor_col The metadata column with donor information.
#' @param Time_col The metadata column with longitudinal information
#' @param Group_col The metadata column with group information. Default is NULL, in which case it will attempt to use celltype as a proxy. Do not run group and cell type at the same time. 
#' @param returnObj Boolean flag for whether or not you want the full PALMO object. Default is FALSE, in which case it returns the variance decomposition data.frame. 
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

linearModeling <- function(Obj, formula, CellType, threshold = 0, NAtoZero = FALSE, numCores = 1){
  
  meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))
  
  if(length(CellType) == 1){
    
    mat1 <- MOCHA::getCellPopMatrix(Obj,CellType,NAtoZero = NAtoZero)
    
    meta <- meta1[meta1$Sample %in% colnames(mat1),]
    rowsToKeep <- rowSums(!is.na(mat1))/dim(mat1)[2] > threshold
    
    mat1 <- mat1[rowsToKeep,]
    
  }else{
    
    allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
    if(NAtoZero){
      allMatrices[is.na(allMatrices)] = 0
    }
    
    
    colnames(allMatrices) <- apply(expand.grid(names(SummarizedExperiment::assays(Obj)), unique(colnames(allMatrices))), 1, 
                                   paste, collapse="__") %>% gsub(" ", "_", .)
    mat1 <- allMatrices[,colSums(allMatrices) != 0]
    
    meta <- parallel::mclapply(1:length(SummarizedExperiment::assays(Obj)), function(x){
      
      meta1 %>% as.data.frame() %>%
        dplyr::mutate(Sample2 = paste(names(SummarizedExperiment::assays(Obj))[x], Sample, sep = "__"),
                      CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
        dplyr::mutate(Sample2 = gsub(" ","_", Sample2), 
                      CellType = gsub(" ", "_", CellType))
      
    }, mc.cores = numCores) %>% do.call('rbind', .)
    
    meta <- meta[meta$Sample2 %in% colnames(mat1),]
    
    rowsToKeep <- rowSums(is.na(mat1))/dim(mat1)[2] > threshold
    
    mat1 <- mat1[rowsToKeep,]
    
  }
  
  #    lmem_res <- lapply(c(1:nrow(mat1)),
  ##        function(x) {
  #           print(paste(x, any(is.na(mat1[x,]))))
  #
  #            df <- data.frame(exp = as.numeric(mat1[x, ]), 
  #                meta, stringsAsFactors = FALSE)
  #            lmerTest::lmer(formula = formula, data = df)
  #        })
  
  #    return(lmem_res)
  
  suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
                                                 function(x) {
                                                   df <- data.frame(exp = as.numeric(mat1[x, ]), 
                                                                    meta, stringsAsFactors = FALSE)
                                                   lmerTest::lmer(formula = formula, data = df)
                                                 }, cl = numCores), classes = "message")
  
  names(lmem_res) = rownames(mat1)
  return(lmem_res)
  
}


modelDeviations <- function(Obj, formula, type = 'z', numCores = 1){
  meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))
  mat1 <- SummarizedExperiment::assays(Obj)[[type]] 
  meta <- meta1[meta1$Sample %in% colnames(mat1),]
  suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
                                                 function(x) {
                                                   df <- data.frame(exp = as.numeric(mat1[x, ]), 
                                                                    meta, stringsAsFactors = FALSE)
                                                   lmerTest::lmer(formula = formula, data = df,
                                                   							 REML = F)
                                                 }, cl = numCores), classes = "message")
  return(lmem_res)
}


subsetDev <- function(Object,
                      subsetBy,
                      groupList){
  sampleData <- SummarizedExperiment::colData(Object)
  if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
    }
  }
  keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]]))]
  return(Object[, keep])
}
