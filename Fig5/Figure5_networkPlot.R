### circlize network
require(data.table)
rm(list=ls())

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 1,
					 track.margin = c(-0.1, 0.1), 
					 points.overflow.warning = FALSE)
#par(mar = rep(0, 6))




mat = fread('~/Downloads/tf_promoter_prediction (9).csv')
mat=mat[,c('TF','pathway','Perc_gene')]

topPathways = mat[,median(Perc_gene),by=pathway]
topPathways = topPathways %>% arrange(V1)
topTFs =  mat[,mean(Perc_gene),by=TF]
topTFs = topTFs %>% arrange(V1)
### change text 

# color palette

mat2=mat[Perc_gene >= 0.33,]
mat2[,sum(Perc_gene >= 0.33), by=TF]
pathway_colors <-  viridis_pal(option = "G")(length(unique(mat2$pathway)))
tf_colors <- rainbow(length(unique(mat2$TF)))
mycolor = c(tf_colors, pathway_colors)


pdf('~/Desktop/Figure 6/tf_network_updated.pdf', height=8,
 		width=15)
 
chordDiagram(
	x = mat2,
	order = c(topTFs$TF, topPathways$pathway),
	grid.col = mycolor,
	transparency = 0.1,
	link.visible = mat2$Perc_gene > 0.2,
	directional = 1,
	#direction.type = c("arrows"), 
	diffHeight  = -0.01,
	annotationTrack = "grid", 
	annotationTrackHeight = c(0.05, 0.1),
	link.arr.type = "big.arrow", 
	link.sort = T, 
	link.arr.width = 0.002,
	link.largest.ontop = T)

# Add text and axis
circos.trackPlotRegion(
	track.index = 1, 
	bg.border = NA, 
	panel.fun = function(x, y) {
		
		xlim = get.cell.meta.data("xlim")
		sector.index = get.cell.meta.data("sector.index")
		
		# Add names to the sector. 
		circos.text(
			x = mean(xlim), 
			y = 3, 
			labels = sector.index, 
			facing = "downward", 
			niceFacing = T,
			cex = 0.75
		)
		
	}
)
dev.off()

