 
simplifiedORA <- function(database, foreground, background){
    
    WebGestaltR::WebGestaltR(enrichMethods = 'ORA', organism ='hsapiens', 
                         enrichDatabase = database,
                         interestGene = foreground, 
                        interestGeneType ="genesymbol",
                        referenceGene =  background, 
                            referenceGeneType= "genesymbol")
    
    }
                                  
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