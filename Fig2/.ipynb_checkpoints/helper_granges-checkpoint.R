############################################################
############################################################
## peak_to_tile 
## 
## 

require(dplyr)

peak_to_tile_macs2_sample <- function(i, tile_list){

        mat <- tile_list[[i]]

        colnames(mat)[1:9] <- c('seqnames','start',
                                        'end','name', 
                                        'score', 'strand',
                                       'signalValue','pValue',
                                       'qValue')


        ## add sample 
        ## specific information
        ## for each peak-call
        mat <- as.data.table(mat)        
        
        diff_start <- mat$start %% 500
        diff_end <- mat$end %% 500

        mat$new_start <- mat$start
        idx_start <- which(diff_start <= 75)
        mat$new_start[idx_start] <- mat$start[idx_start] + 75
        
        mat$new_end <- mat$end 
        idx_end <- which(diff_end <= 75)
        mat$new_end[idx_end] <- mat$new_end[idx_end] - 75
        
        colnames(  mat )[2:3] <-c('old_start','old_end')

        mat <- mat[,c('seqnames', 
                      'old_start','new_start',
                      'old_end','new_end'),
                   with=F]


        tmp <- makeGRangesFromDataFrame(mat, start.field='new_start',
                                 end.field='new_end', keep.extra.columns=TRUE)

        ### re-align to 
        ### nearest 500bp by
        ### rounding up and down 
        ### using MP's alignment
        ### from dynamic bins

        blackList <- getBlacklist(ArchRProj)

        #Let's subtract out areas of fragments that overlap with blacklist regions.
        #Let's not remove the region entirely, because the regions may be long or only
        binSize=500
        TotalRangesFilt<- plyranges::setdiff_ranges(tmp, blackList)

        RangeBins <- plyranges::stretch(plyranges::anchor_end(TotalRangesFilt),
                                      extend = GenomicRanges::start(TotalRangesFilt)%% binSize)

        FinalBins <- plyranges::stretch(plyranges::anchor_start(RangeBins),
                                      extend = (binSize - end(RangeBins)%% binSize)) %>%
        plyranges::reduce_ranges() %>% plyranges::slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

        tmp2 <- as.data.table(FinalBins)

        tmp2$tileID <- paste(tmp2$seqnames,':',
                               tmp2$start,'-',
                               tmp2$end,
                               sep='')

        return(tmp2)

    }


peak_to_tile_macs2_bulk<- function(i, tile_list){

        mat <- tile_list[[i]]

        colnames(mat)[1:9] <- c('seqnames','start',
                                        'end','name', 
                                        'score', 'strand',
                                       'signalValue','pValue',
                                       'qValue')


        ## add sample 
        ## specific information
        ## for each peak-call
        mat <- as.data.table(mat)
        mat$name <-  gsub('../downsample_macs2_peaks/',
                          '',
                          mat$name)
        
        
        diff_start <- mat$start %% 500
        diff_end <- mat$end %% 500

        mat$new_start <- mat$start
        idx_start <- which(diff_start <= 75)
        mat$new_start[idx_start] <- mat$start[idx_start] + 75
        
        mat$new_end <- mat$end 
        idx_end <- which(diff_end <= 75)
        mat$new_end[idx_end] <- mat$new_end[idx_end] - 75
        
        
        
        colnames(  mat )[2:3] <-c('old_start','old_end')

        mat <- mat[,c('seqnames', 
                      'old_start','new_start',
                      'old_end','new_end'),
                   with=F]

        tmp <- makeGRangesFromDataFrame(mat, start.field='new_start',
                                 end.field='new_end', keep.extra.columns=TRUE)

        ### re-align to 
        ### nearest 500bp by
        ### rounding up and down 
        ### using MP's alignment
        ### from dynamic bins

        blackList <- getBlacklist(ArchRProj)

        #Let's subtract out areas of fragments that overlap with blacklist regions.
        #Let's not remove the region entirely, because the regions may be long or only
        binSize=500
        TotalRangesFilt<- plyranges::setdiff_ranges(tmp, blackList)

        RangeBins <- plyranges::stretch(plyranges::anchor_end(TotalRangesFilt),
                                      extend = GenomicRanges::start(TotalRangesFilt)%% binSize)

        FinalBins <- plyranges::stretch(plyranges::anchor_start(RangeBins),
                                      extend = (binSize - end(RangeBins)%% binSize)) %>%
        plyranges::reduce_ranges() %>% plyranges::slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

        tmp2 <- as.data.table(FinalBins)
        tmp2$tileID <- paste(tmp2$seqnames,':',
                               tmp2$start,'-',
                               tmp2$end,
                               sep='')

        return(tmp2)

    }



############################################################
############################################################
## peak_to_tile 
## 
## 

peak_to_tile_homer_sample <- function(i, tile_list){

        mat <- tile_list[[i]]

        colnames(mat)[1:6] <- c('peakID', 'seqnames','start',
                                'end','strand','score')
    
        #mat$score <- mat$score / 10000


        ## add sample 
        ## specific information
        ## for each peak-call
        mat$sampleName <- samples[i]
        mat <- as.data.table(mat)
           
        diff_start <- mat$start %% 500
        diff_end <- mat$end %% 500

        mat$new_start <- mat$start
        idx_start <- which(diff_start <= 75)
        mat$new_start[idx_start] <- mat$start[idx_start] + 75
        
        mat$new_end <- mat$end 
        idx_end <- which(diff_end <= 75)
        mat$new_end[idx_end] <- mat$new_end[idx_end] - 75
        
        colnames(  mat )[3:4] <-c('old_start','old_end')

        mat <- mat[,c('seqnames', 'sampleName', 
                      'old_start','new_start',
                      'old_end','new_end'),
                   with=F]


        tmp <- makeGRangesFromDataFrame(mat, start.field='new_start',
                                 end.field='new_end', keep.extra.columns=TRUE)

        ### re-align to 
        ### nearest 500bp by
        ### rounding up and down 
        ### using MP's alignment
        ### from dynamic bins

        blackList <- getBlacklist(ArchRProj)

        #Let's subtract out areas of fragments that overlap with blacklist regions.
        #Let's not remove the region entirely, because the regions may be long or only
        binSize=500
        TotalRangesFilt<- plyranges::setdiff_ranges(tmp, blackList)

        RangeBins <- plyranges::stretch(plyranges::anchor_end(TotalRangesFilt),
                                      extend = GenomicRanges::start(TotalRangesFilt)%% binSize)

        FinalBins <- plyranges::stretch(plyranges::anchor_start(RangeBins),
                                      extend = (binSize - end(RangeBins)%% binSize)) %>%
        plyranges::reduce_ranges() %>% plyranges::slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

        tmp2 <- as.data.table(FinalBins)
        tmp2$sampleName <- samples[i]  

        tmp2$tileID <- paste(tmp2$seqnames,':',
                               tmp2$start,'-',
                               tmp2$end,
                               sep='')

        return(tmp2)

    }


peak_to_tile_homer_bulk <- function(i, tile_list){

        mat <- tile_list[[i]]

        colnames(mat)[1:6] <- c('peakID', 'seqnames','start',
                                'end','strand','score')
    
        ## add sample 
        ## specific information
        ## for each peak-call
        mat <- as.data.table(mat)
           
        diff_start <- mat$start %% 500
        diff_end <- mat$end %% 500

        mat$new_start <- mat$start
        idx_start <- which(diff_start <= 75)
        mat$new_start[idx_start] <- mat$start[idx_start] + 75
        
        mat$new_end <- mat$end 
        idx_end <- which(diff_end <= 75)
        mat$new_end[idx_end] <- mat$new_end[idx_end] - 75
        
        colnames(  mat )[3:4] <-c('old_start','old_end')

        mat <- mat[,c('seqnames', 
                      'old_start','new_start',
                      'old_end','new_end'),
                   with=F]


        tmp <- makeGRangesFromDataFrame(mat, start.field='new_start',
                                 end.field='new_end', keep.extra.columns=TRUE)

        ### re-align to 
        ### nearest 500bp by
        ### rounding up and down 
        ### using MP's alignment
        ### from dynamic bins

        blackList <- getBlacklist(ArchRProj)

        #Let's subtract out areas of fragments that overlap with blacklist regions.
        #Let's not remove the region entirely, because the regions may be long or only
        binSize=500
        TotalRangesFilt<- plyranges::setdiff_ranges(tmp, blackList)

        RangeBins <- plyranges::stretch(plyranges::anchor_end(TotalRangesFilt),
                                      extend = GenomicRanges::start(TotalRangesFilt)%% binSize)

        FinalBins <- plyranges::stretch(plyranges::anchor_start(RangeBins),
                                      extend = (binSize - end(RangeBins)%% binSize)) %>%
        plyranges::reduce_ranges() %>% plyranges::slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

        tmp2 <- as.data.table(FinalBins)
        tmp2$tileID <- paste(tmp2$seqnames,':',
                               tmp2$start,'-',
                               tmp2$end,
                               sep='')

        return(tmp2)

    }




############################################################
############################################################
## identify adjacent regions 

# adjacency <- function(chr, cd14_scMACS){
#     tmp_mat <- cd14_scMACS[seqnames==chr]
#     tmp_mat <- tmp_mat[order(tmp_mat$start, decreasing=F),]
    
#     adjs <- which( (tmp_mat$end +1) %in% tmp_mat$start) 
    
#     tmp_mat$open_adjacent <- 0
#     tmp_mat$open_adjacent[adjs] <- 1
    
#     consecutive_peaks <- rle(tmp_mat$open_adjacent)
#     tmp_mat$diffs <- c(0,diff(tmp_mat$open_adjacent ))
    
#     i=1
#     k=0
    
#     j=numeric(nrow(tmp_mat))
    
#     while(i < nrow(tmp_mat)){
        
#                     if(tmp_mat$open_adjacent[i] ==0){
#                         i =i+1
#                     } else if(tmp_mat$open_adjacent[i]==1 & tmp_mat$diffs[i]==1 ){
#                         k=k+1
#                         j[i] <- k
#                         i=i+1
#                     } else if(tmp_mat$open_adjacent[i]==1 & tmp_mat$diffs[i]==0 ){
#                         j[i] <- k
#                         i=i+1


#                     }
        
#             }
        
#     tmp_mat$block_id <- paste(chr,'block',j,sep='-')
#     return(tmp_mat)    
# }

# chrNames <- unique(cd14_scMACS$seqnames)   
# cd14_scMACS_list <-  mclapply(chrNames,
#          function(x)
#              adjacency(x, cd14_scMACS),
#          mc.cores=20
#          )
# cd14_scMACS <- rbindlist(cd14_scMACS_list)
    
    
############################################################
############################################################



subdivideIRanges <- function(x,subsize=100) {
  if (length(x) == 0) {
    return(x)
  }
  start.pos <- start(x)
  end.pos <- end(x)
  widths <- width(x)
  nsubranges <- pmax(1,floor(widths/subsize))
  out.start.pos <- numeric(sum(nsubranges))
  out.end.pos <- numeric(sum(nsubranges))
  out.idx <- 1
  for (i in 1:length(x)) {
    if (widths[i] < 2*subsize) {
      out.start.pos[out.idx] <- start.pos[i]
      out.end.pos[out.idx] <- end.pos[i]
      out.idx <- out.idx + 1
    } else {
      sr <- subdivide.range(start.pos[i],end.pos[i],subsize)
      out.start.pos[out.idx:(out.idx+sr$length-1)] <- sr$start.pos
      out.end.pos[out.idx:(out.idx+sr$length-1)] <- sr$end.pos
      out.idx <- out.idx + sr$length
    }
  }
  IRanges(start=out.start.pos,end=out.end.pos)
}

subdivideGRanges <- function (x, subsize=100) {
  if (length(x) == 0) {
    return(x)
  }
  if (length(subsize) > 1) {
    stop("The subsize argument should be a single number: the desired width of the subdivided ranges")
  }
  x <- sort(reduce(x))
  gr_list <- lapply(levels(seqnames(x)), function(seqlvl) {
    if (!any(seqnames(x) == seqlvl)) {
      return(GRanges())
    }
    rg <- ranges(x[seqnames(x) == seqlvl])
    GRanges(seqnames = seqlvl, ranges = subdivideIRanges(rg,subsize), seqlengths = seqlengths(x))
  })
  do.call(c, gr_list)
}

subdivide.range <- function(start.pos,end.pos,subsize=100) {
  width <- end.pos - start.pos + 1
  if (width < 2*subsize) {
    stop("Width is less than 2 times subsize")
  }
  nchunks <- floor(width/subsize)
  relative.starts <- round(0:(nchunks-1)*width/nchunks)
  relative.ends <- round(1:nchunks*width/nchunks)-1
  return(list(start.pos=start.pos+relative.starts,end.pos=start.pos+relative.ends,length=nchunks))
}










############################################################
############################################################

### annotate tiles 
### by category

annotateTiles <- function(tileGRanges,
                          TxDb = NULL,
                          Org = NULL,
                          promoterRegion = c(2000, 100)) {
  
  txList <- suppressWarnings(GenomicFeatures::transcriptsBy(TxDb, by = ("gene")))
  names(txList) <- suppressWarnings(AnnotationDbi::mapIds(Org, names(txList), "SYMBOL", "ENTREZID"))

  txs <- IRanges::stack(txList) %>%
    GenomicRanges::trim() %>%
    S4Vectors::unique(.)

  promoterSet <- IRanges::stack(txList) %>%
    GenomicRanges::trim(.) %>%
    S4Vectors::unique(.) %>%
   GenomicRanges::promoters(., upstream = promoterRegion[1], downstream = promoterRegion[2])


  getOverlapNameList <- function(rowTiles, annotGR) {
    overlapGroup <- IRanges::findOverlaps(rowTiles, annotGR) %>% as.data.frame()
    overlapGroup$Genes <- as.character(annotGR$name[overlapGroup$subjectHits])
    last <- overlapGroup %>%
      dplyr::group_by(.data$queryHits) %>%
      dplyr::summarize(Genes = paste(unique(.data$Genes), collapse = ", "))
    return(last)
  }

  txs_overlaps <- getOverlapNameList(tileGRanges, txs)
  promo_overlaps <- getOverlapNameList(tileGRanges, promoterSet)

  tileType <- as.data.frame(GenomicRanges::mcols(tileGRanges)) %>%
    dplyr::mutate(Index = 1:nrow(.)) %>%
    dplyr::mutate(Type = dplyr::case_when(
      Index %in% promo_overlaps$queryHits ~ "Promoter",
      Index %in% txs_overlaps$queryHits ~ "Intragenic",
      TRUE ~ "Distal"
    )) %>%
    dplyr::left_join(promo_overlaps, by = c("Index" = "queryHits")) %>%
    dplyr::rename("Promo" = Genes) %>%
    dplyr::left_join(txs_overlaps, by = c("Index" = "queryHits")) %>%
    dplyr::rename("Txs" = Genes) %>%
    dplyr::mutate(Genes = ifelse(Type == "Promoter", Promo, NA)) %>%
    dplyr::mutate(Genes = ifelse(Type == "Intragenic", Txs, Genes))

  tileGRanges$tileType <- tileType$Type
  tileGRanges$Gene <- tileType$Genes

    return(tileGRanges)

}
