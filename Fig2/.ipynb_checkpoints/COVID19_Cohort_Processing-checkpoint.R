#This code takes all the scATAC cells, runs them through dimensionality reduction and labels them according to major cell type IDs.
#After that, it subset out T cells, NK cells, B cells, and monocytes/other cells into seperate projects for further analysis.
# By subseting each major type(s) into a seperate project, we hope to identify 
# variable genes between subsets of those types and be able to pursue more accurate 
# cell type labeling. Additionally, it cut down on memory requirements.

library(dplyr)
library(ArchR)
library(Seurat)
library(readxl)
library(tidyverse)

source("../../scATAC_functions.R")

addArchRGenome('hg38')
addArchRThreads(45)

ource("../../scATAC_functions.R")

addArchRGenome('hg38')
addArchRThreads(45)

#Load metadata
moreMeta <- read.csv("20210318_final_meta_merged_df.csv", stringsAsFactors = F) %>% as.data.frame()

#Pull in pre-processed arrows, generate doublet scores, and save into an ArchRProject
ArrowFiles <- list.files(path = "CovidAllArrows", pattern = ".arrow")
doubScores <- addDoubletScores(paste("CovidAllArrows/", ArrowFiles,sep = ""), k = 10,knnMethod = "UMAP",LSIMethod = 1)
FullCovid <- ArchRProject(paste("CovidAllArrows/", ArrowFiles,sep = ""), outputDirectory ="FullCovid/", copyArrows=TRUE)

#Load in peakset fro Buenrostro et al for initial dimensionality reduction. See included GRanges objects ("hg38_peaks_gr.rds").
buenpeaks <- readRDS(file="hg38_peaks_gr.rds")
FullCovid <- addFeatureMatrix(FullCovid,features=buenpeaks)

#Filter out doublets
FullCovid <- filterDoublets(FullCovid, filterRatio = 8)


options(future.globals.maxSize= 1000485760000)
FullCovid <- addIterativeLSI(
    ArchRProj = FullCovid,
    useMatrix = "FeatureMatrix", 
    varFeatures = 10000,
    name = "IterativeLSI",
    iterations = 2,
    force = TRUE
)

FullCovid <- addClusters(
    input = FullCovid,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 3,
    force = TRUE
)

FullCovid <- addUMAP(ArchRProj = FullCovid, reducedDims = "IterativeLSI", force = TRUE)

saveArchRProject(FullCovid, "FullCovid/")

#Loads a reference Seurat object, from their standard reference dataset.
reference <- readRDS("reference.RDS")

#Conduct unconstrained integration to get approximate L1 labels
FullCovid <- addGeneIntegrationMatrix(
    ArchRProj = FullCovid,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    transferParams = list(dims = 1:10, k.weight = 20),
    seRNA = reference,
    addToArrow = FALSE,
    groupRNA = "celltype.l1",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un", force=TRUE
)

#Let's do a quick dive for QC of clusters and the dataset
totalList <- names(table(FullCovid$predictedGroup_Un))

NK <- totalList[grepl("NK", totalList)]
TCell <- totalList[grepl("CD4|CD8|T", totalList)]
MonoDC <- totalList[grepl("Mono|DC|other$",  totalList)]
BCell <- totalList[grepl("B|Plasma",  totalList)]
length(totalList)
length(NK) + length(TCell) +length(MonoDC) + length(BCell)

CellTypes <- list(TCell, NK, MonoDC, BCell)
CellTypeNames <- c("TCell","NK", "MonoDC", "BCell")
typePredic = c("predictedGroup_Un","predictedScore_Un")

QCThat(FullCovid, "FullCovid_Labeling", "Clusters",CellTypes, CellTypeNames, typePredic, "ridges")
dev.off()

#Now, let's calculate the ideal kmeans for cluster subsetting
registerDoParallel(45)
UMAP <- getEmbedding(ArchRProj = FullCovid, embedding = "UMAP", returnDF = TRUE)
best <- unlist(OptimalKMeans(FullCovid, UMAP, "KMeans", CellTypes, CellTypeNames, typePredic))
#Double check optimal mixing value
best[which.max(best)]
kClust <- kmeans(UMAP, which.max(best)+2, iter.max = 200, nstart = 3)
FullCovid <- addCellColData(FullCovid, data = as.character(kClust$cluster), name = "KMeans",cells = names(kClust$cluster), force = TRUE)

#Look at QC from best K-means clustering
QCThat(FullCovid, "FullCovid_KmeansLabeling", "KMeans", CellTypes, CellTypeNames, typePredic, "ridges")

saveArchRProject(FullCovid, "FullCovid/")
FullCovid <- loadArchRProject("FullCovid/")

#Now let's isolate the major label for each cluster.
#if the major label for a cluster constitutes less than 80% of the cells of that cluster,
#we will label that cluster as "Mixed"
cMM = as.matrix(confusionMatrix(FullCovid$KMeans, FullCovid$predictedGroup_Un))
preClust <- colnames(cMM)[apply(cMM, 1, which.max)]
names(preClust) = names(apply(cMM, 1, which.max))

#Now let's break out each major grouping we want for subsetting. 
#Cell label groups were defined above during QC.
#Then we identify which louvain clusters match to which labels.
#Then we create a list of cells in each group of clusters, with groups defined for subsetting.
id_clusters <- list("NK" = NK,
		     "MonoDC" = MonoDC,
                     "BCell" = BCell,
		     "TCell" = TCell)

atac_clusters <- lapply(id_clusters,
                        function(x) {
                          names(preClust)[preClust %in% x]
                        })

atac_cells <- lapply(atac_clusters,
                     function(x) {
                       FullCovid$cellNames[FullCovid$KMeans %in% x]
                     })

#use the following outputs to verify that atac_cells and groups were pulled in correctly before 
#the next several steps.
test <- lapply(atac_cells,
                     function(x) {
                       FullCovid$KMeans[FullCovid$cellNames %in% x]
                     })

test1 <- lapply(atac_cells,
                     function(x) {
                       FullCovid$predictedGroup_Un[FullCovid$cellNames %in% x]
                     })


#Let's subset to each major group, creating a new ArchR Project for each one, 
#according to the names under id_clusters.
ArchRList <- lapply(names(id_clusters),
                    function(x) {
                      subsetArchRProject(FullCovid, cells = atac_cells[x][[1]], 
		       outputDirectory = paste(x,"/", sep = ""), force = TRUE)
                    })

lapply(ArchRList, function(x) {table(x$predictedGroup_Un)})

#Let's save FullCovid, and then remove it to save on memory.
saveArchRProject(FullCovid, "FullCovid/")
rm(FullCovid)

#Let's run iteratureLSI, Louvain clustering, and UMAP on each subset.
ArchRList <- lapply(ArchRList,
                    function(x) {
                      addIterativeLSI(x,
		    	useMatrix = "TileMatrix", 
                    	name = "IterativeLSI", 
		   	iterations = 6,
                   	force = TRUE)
                    })


gc()

ArchRList <- lapply(ArchRList,
                    function(x) {
                      addClusters(x,
		     	reducedDims = "IterativeLSI", 
                     	name = "Clusters2", 
			resolution = 3,
                     	force = TRUE)
                    })


gc()

ArchRList <- lapply(ArchRList,
                    function(x) {
                      addUMAP(x,
		     	reducedDims = "IterativeLSI", 
                     	force = TRUE)
                    })

#Let's save our progress.
names(ArchRList) = names(id_clusters)

CellTypeNames <- c("NK", "MonoDC", "BCell","TCell")
ArchRList <- lapply(1:length(CellTypeNames), 
	function(x) { 
		loadArchRProject(paste(CellTypeNames[x],"/", sep = "")) 
	})

names(ArchRList) = CellTypeNames
gc()

#Let's subset our seurat reference by cell type for integration.
Idents(reference) <- reference@meta.data$celltype.l1
totalList <- names(table(Idents(reference)))
NK <- totalList[grepl("NK", totalList)]
TCell <- totalList[grepl("CD4|CD8|T", totalList)]
MonoDC <- totalList[grepl("Mono|DC|other$",  totalList)]
BCell <- totalList[grepl("B|Plasma",  totalList)]
id_clusters <- list("NK" = NK,
		     "MonoDC" = MonoDC,
                     "BCell" = BCell,
		     "TCell" = TCell)

SeuratList <- lapply(1:length(id_clusters),
		function(x) {
			subset(reference, idents = id_clusters[[x]])
		})
rm(reference)
saveRDS(SeuratList,file = "SeuratReferenceList.RDS")

#Update level 2 anotations to level 2.5 and remove cell types that don't appear in scATAC.
SeuratList1 <- lapply(SeuratList, function(x) {
			l2.5 <- as.character(x@meta.data$celltype.l2.5)
			l3 <- as.character(x@meta.data$celltype.l3)
			l2.5[l3 %in% c("CD8 TEM_4","CD8 TEM_5")] = "CD8 TEMRAS"
			x@meta.data$celltype.l2.5 = l2.5
			Idents(x) <- factor(l2.5)
			x[,which(!(Idents(SeuratList[[3]]) %in% c("Eryth","Platelet","Doublet")))]
		})

saveRDS(SeuratList2,file = "SeuratReferenceList.RDS")

#Let's integration (and thus label) according to level l2.5 identities for each subset.
ArchRList <- lapply(c(1:length(ArchRList)),
                    function(x) {
			 addGeneIntegrationMatrix(ArchRList[[x]],
				useMatrix = "GeneScoreMatrix",
				matrixName = "GeneIntegrationMatrix",
				reducedDims = "IterativeLSI",
				seRNA = SeuratList[[x]],
				addToArrow = FALSE,
    				transferParams = list(dims = 1:10, k.weight = 20),
    				nGenes = 4000,
				groupRNA = "celltype.l2.5",
				nameCell = "predictedCell_Co2",
				nameGroup = "predictedGroup_Co2",
				nameScore = "predictedScore_Co2",
        			force = TRUE)
			})

#Save each ArchRProject within ArchRList.
lapply(1:length(ArchRList), 
	function(x) { 
		saveArchRProject(ArchRList[[x]], 
				 outputDirectory = paste(names(ArchRList)[x],"/", sep = "")) 
	})
  
#Generate QC plots for the integration on each subset

GroupType = c("predictedGroup_Co2","predictedScore_Co2")
lapply(1:length(ArchRList), 
	function(x) { 
		QCThat(ArchRList[[x]], 
			paste(names(ArchRList)[x],"_InitialLabeling", sep = ""), "Clusters2",
			names(table(getCellColData(ArchRList[[x]], "predictedGroup_Co2")$predictedGroup_Co2)),
			names(table(getCellColData(ArchRList[[x]], "predictedGroup_Co2")$predictedGroup_Co2)),
			GroupType,
			"ridges")
		Sys.sleep(1)
	})
