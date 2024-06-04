

library(MOCHA)
library(SummarizedExperiment)
library(ArchR)


### Single Cell UMAPs
fullCOVID <- loadArchRProject('FullCovidf')

pdf('SingleCell_Batch.pdf')
plotEmbedding(fullCOVID, name = 'CellTypes', alpha = 0.3, size =0.05)

plotEmbedding(fullCOVID, name = 'AIFI.Batch', alpha = 0.3, size =0.05)
dev.off()

umap1 <- getEmbedding(fullCOVID)
metadata <- getCellColData(fullCOVID)
all(rownames(metadata) == rownames(umap1))
output <- cbind(umap1, metadata[, c('CellTypes', 'AIFI.Batch')])
write.csv(output, 'UMAP_DF.csv')

