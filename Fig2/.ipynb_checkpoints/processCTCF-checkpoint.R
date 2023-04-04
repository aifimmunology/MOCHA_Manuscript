library(plyranges)

## Download CTCF ChIP-seq peak calls from all Blood cell types. 
#https://chip-atlas.org/peak_browser

#Unfortunately, they download with the same name, so I labeled the Hg19 one as Hg19. 


Hg38 <- plyranges::read_bed('Oth.Bld.05.CTCF.AllCell.bed')
plyranges::reduce_ranges(Hg38) %>% plyranges::write_bed('All_Blood_CTCF.bed')

Hg19 <- plyranges::read_bed('Oth.Bld.05.CTCF.AllCell_Hg19.bed')
plyranges::reduce_ranges(Hg19) %>% plyranges::write_bed('All_Blood_CTCF_hg19.bed')
