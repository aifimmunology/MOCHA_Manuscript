library(Seurat)
library(Signac)
library(tidyverse)

########## Run time analysis for Seurat
future::plan("multicore", workers = 60)
options(future.globals.maxSize = 400 * 1024 ^ 3)

## Load CD16 Signac object from Stability analysis script (Fig 3)
CD16s <- readRDS('CD16_Mono_SignacObject.RDS')
Idents(CD16s) = CD16s@meta.data$EarlyComparison

SeurDaps_RunTime <- lapply(seq_along(c(1:10)), function(x){

    print(x)
    startTime = Sys.time()

    seurM <- FindMarkers(
          object = CD16s,
          ident.1 = "CD16.Mono.Positive",
          ident.2 = "CD16.Mono.Negative",
          min.pct = 0.001,
          logfc.threshold = 10^-(x/3),
          test.use = 'LR',
      latent.vars = 'nFrags'
    )
    runTime <- Sys.time() - startTime
    seurM %>% mutate(RunTime = runTime, logfcThreshold = 10^-(x/3))

})

SeurRunTime <- do.call('rbind',SeurDaps_RunTime)

SeurRunSum <- SeurRunTime %>% group_by(logfcThreshold) %>%
                        summarize(TotalPeaks = dplyr::n(),
                                  TotalDAPs= sum(p_val_adj < 0.05),
                                 RunTime = unique(RunTime))

write.csv(SeurDaps_RunTime[[10]], 'Seurat_204K_PeaksTest.csv')
write.csv(SeurRunSum, 'Seurat_RunTimeAnalysis.csv')