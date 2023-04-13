##DAP pathway analysis: Reactome pathway
enrichDataBaseList <- WebGestaltR::listGeneSet() 
allGenes <- names(genes(TxDb.Hsapiens.UCSC.hg38.refGene, single.strand.genes.only = FALSE)) 
refList <- AnnotationDbi::mapIds(org.Hs.eg.db, allGenes, "SYMBOL", "ENTREZID") 

DAPs_geneSet <- plyranges::filter(daps, FDR < 0.2) %>% 
	MOCHA::annotateTiles(., TxDb.Hsapiens.UCSC.hg38.refGene, org.Hs.eg.db) %>% 
	plyranges::filter(tileType == 'Promoter') %>% mcols(.) %>% as.data.frame() %>% 
	dplyr::select(Gene) %>% unlist() %>% paste0(., collapse =", ") %>% 
	stringr::str_split(., pattern = ', ') %>% unlist() %>% unique() 

DAP_Reactome <- simplifiedORA(enrichDataBaseList$name[9], as.character(DAPs_geneSet), refList) 
Hierarch <- read.csv('ReactomePathwaysRelation.csv', col.names = c('V1', 'V2')) 
ReactID <- read.csv('ReactomePathways.csv', col.names = c('V1', 'V2','V3')) 
ReactID$V2 <- sub(" $","",ReactID$V2) 
PathAnnot = lapply(DAP_Reactome$description, function(y) 
		unlist(lapply(y, function(x) findTrunk(x, Hierarch, ReactID)))) %>% 
		unlist() 
DAP_Reactome$TopLevel = PathAnnotDAP_Reactome$SecondLevel = c(rep('Innate Immune System', times = 10), 
		rep('Cytokine Signaling in Immune System', n = 1), 
		rep('Innate Immune System', times = 6), 
		rep('Adaptive Immune System', times = 2), 
		'Metabolism', 
		rep('Signal Transduction', times = 1), 
		rep('Disease', times = 2), 
		'Innate Immune System', 
		'Gene Expression (Transcription)', 
		'Adaptive Immune System', 'Circadian Clock'))

DAP_Reactome <- arrange(DAP_Reactome, SecondLevel, enrichmentRatio) 
DAP_Reactome$description = factor(DAP_Reactome$description, levels = unique(DAP_Reactome$description)) 

write.csv(DAP_Reactome, 'SuppDataFile1_MOCHA_ReactomePathways.csv') 
PathAnnot = lapply(tmp1$description, function(y) unlist(lapply(y, function(x) findTrunk(x, Hierarch, ReactID)))) %>% unlist() 
PathAnnot[grepl('external stimuli',PathAnnot)] = 'Cellular responses to stimuli'
PathAnnot[grepl('Mitotic',PathAnnot)] = 'Cell Cycle' 
PathAnnot[grepl('Diseases of signal transduction',PathAnnot)] = 'Disease'

PathAnnot2 = lapply(tmp2$description, function(y) unlist(lapply(y, function(x) findTrunk(x, Hierarch, ReactID)))) %>% unlist() 
PathAnnot2[grepl('beta-catenin|CTNNB1',PathAnnot2)] = 'Disease'
PathAnnot2[grepl('Rho',PathAnnot2)] = 'Metabolism'


pdf('Fig3_b.pdf') 
ggplot(DAP_Reactome, aes(x = description, y = enrichmentRatio, fill = SecondLevel, Group = SecondLevel)) + 
	geom_col() + coord_flip() + theme_bw() + xlab('Term Description') + ylab('Term Enrichment') +
	 theme(legend.position = c(0.8, 0.8))+ ggtitle('Reactome Pathway Enrichment in All Genes with Differential Promoter') 
dev.off()