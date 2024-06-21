

## I'm going to explore the ChIP-seq data:
    #Write a function that will annotate peakset with ChIp-Seq validation
    
    
    
library(plyranges)
library(ArchR)

if(FALSE){
PedSen_d0 <- loadArchRProject(PedSen_d0)

peaks <- getPeakSet(PedSen_d0)

bloodMarks <- read_bed(" ", genome_info = "hg38", overlap_ranges)


CD14Mono <- read_bed("Oth.Bld.05.AllAg.Monocytes-CD14PULUS.bed", 
                        genome_info = "hg38")
}                      
############### Here are functions for processing useful in from the bed metadata file.
## Take in GRanges object
## Take out metadata column by name (column),
## Parse out metadata within that name by a splitting character
## Remove misc characters that might have been added. 
## Removes characters in order. 
processChipMeta <- function(ChIPGranges, column, splitCharacter, removeCharacters){

    tmp <- unique(mcols(ChIPGranges)[,column])
    tmp1 <- strsplit(tmp, split = splitCharacter)

    namesList <- lapply(c(1:length(tmp1)), function(x){

                tmp2 <- tmp1[[x]]
                names(tmp2) <- iterateRemove (tmp2, removeCharacters) %>%  gsub("=.*","", .)
                tmp2 <- iterateRemove(tmp2, removeCharacters) %>% gsub(".*=","", .)
                as.list(tmp2)[c(1:8)]

    })

    df <- rbindlist(namesList, fill = TRUE)  %>% as.data.frame(.)
    df <- df[,!is.na(colnames(df))]
    
    edit <- ChIPGranges
    meta <- mcols(edit) %>% as.data.frame(.)
    mcols(edit) <- dplyr::left_join(meta, cbind(OriginalName = tmp, df), 
                                    by= c('name' = 'OriginalName'))
    return(edit)
}

## This function iterates over a list of strings, removing them from a longer string in order. 
iterateRemove <- function(string, listReplace){

    tmp <- gsub(listReplace[1],"", string)
    if(length(listReplace) > 1){
        return(iterateRemove(tmp, listReplace[-1]))
    }else{
        
        return(tmp)
    
    }

}



if(FALSE){
## Now let's process 
toRemove = c("%20", "%","<br>", "2B|B3","B3")
CD14Mono1 <- processChipMeta(CD14Mono, "name", ";", toRemove)
head(CD14Mono1)
write_bed(CD14Mono1, file = "ChIPseq_TF_CD14Monocytes.bed")


Treg <- read_bed("Oth.Bld.05.AllAg.Treg.bed", 
                        genome_info = "hg38")
head(Treg$name)
toRemove = c("%20", "%","<br>", "2B|B3","B3")
Treg1 <- processChipMeta(Treg, "name", ";", toRemove)
head(Treg1)
write_bed(Treg1, file = "ChIPseq_TF_Treg.bed")

DCs <- read_bed("Oth.Bld.05.AllAg.Dendritic_Cells.bed", 
                        genome_info = "hg38")
DCs1 <- processChipMeta(DCs, "name", ";", toRemove)
head(DCs1)
write_bed(DCs1, file = "ChIPseq_TF_Dendritic_Cells.bed")

MemT <- read_bed("Oth.Bld.05.AllAg.Memory_T_cells.bed", 
                        genome_info = "hg38")
MemT1 <- processChipMeta(MemT, "name", ";", toRemove)
head(MemT1)
write_bed(MemT1, file = "ChIPseq_TF_MemoryTCells.bed")

BCells <- read_bed("Oth.Bld.05.AllAg.B_cells.bed", 
                        genome_info = "hg38")
BCells1 <- processChipMeta(BCells, "name", ";", toRemove)
head(BCells1)
write_bed(BCells1, file = "ChIPseq_TF_BCells.bed")

Mono <- read_bed("Oth.Bld.05.AllAg.Monocytes.bed", 
                        genome_info = "hg38")
Mono1 <- processChipMeta(Mono, "name", ";", toRemove)
head(Mono1)
write_bed(Mono1, file = "ChIPseq_TF_Mono.bed")
    
    }