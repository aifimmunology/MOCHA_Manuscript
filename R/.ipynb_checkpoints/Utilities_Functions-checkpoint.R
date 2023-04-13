#' @title \code{simplifiedORA}
#'
#' @description \code{simplifiedORA} Quick wrapper for using WebGestaltR's ORA methods. 
#'
#' @param database String, referencing which database you want to query. 
#' @param foreground foreground gene set. Should be gene symbols
#' @param background background gene set. Should be gene symbols. 
#'
#' @return
#'
#' @noRd
#'

simplifiedORA <- function(database, foreground, background){
    
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                            referenceGeneType= "genesymbol")
    
    }

#' @title \code{findTrunk}
#'
#' @description \code{findTrunk} Wrapper for leveragign Reactome's pathway hierarchy to annotate enriched pathways. 
#'
#' @param pathway pathway name
#' @param TreeStructure File from reactome containing the pathway hiearchy
#' @param IDMatrix The matrix that contains reactome pathway names and matching IDs. 
#' @param exportTree Boolean for whether to return just the last node, or the whole tree from each pathway. 
#' @return TrunKDescription
#'
#' @noRd
#'

findTrunk <- function(pathway, TreeStructure, IDMatrix, exportTree = FALSE){

    specID <- IDMatrix[match(pathway,IDMatrix[,2]),1]

    if(!(any(TreeStructure[,2] %in% specID))){

        return(pathway)

    }

    treeID <- TreeStructure[which(TreeStructure[,2] %in% specID)[1],]
    descriptionTree = pathway

    while(any(TreeStructure$V2 %in% treeID$V1)){

        treeID <- TreeStructure[which(TreeStructure$V2 %in% treeID$V1)[1],]
        descriptionTree <- paste(descriptionTree,  
                                 IDMatrix[match(treeID[,1],IDMatrix[,1])[1],2], sep = ", ")

    }

    TrunkDescription <- IDMatrix[match(treeID[,1],IDMatrix[,1]),2]


    if(exportTree){
    
        return(descriptionTree)
        
    }
    return(TrunkDescription)

}