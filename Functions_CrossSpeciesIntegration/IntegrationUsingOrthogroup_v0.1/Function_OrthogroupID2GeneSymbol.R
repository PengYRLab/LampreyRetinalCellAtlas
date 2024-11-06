#' @description  This function converts the orthogroup IDs to gene symbols
#' @param mat Expression matrix
#' @param orthgroup_table The orthogroup table used
#' @param reference_species The species in the orthogroup table 
#' @author Junqiang Wang
#' @export
OrthogroupID2GeneSymbol<-function(mat,
                                  orthogroup_table,
                                  reference_species){
  
  require(stringr)
  
  common<-intersect(rownames(mat), rownames(orthogroup_table))
  
  mat<-mat[match(common, rownames(mat)), ]
  
  orthogroup_table<-orthogroup_table[match(common, rownames(orthogroup_table)), ]
  
  
  tmp<-orthogroup_table[, reference_species]
  
  if(length(reference_species)>1){
    
 tmp<-apply(tmp, 1, FUN = function(x){
    x<-as.character(x)
    x<-paste0(x, collapse=":")
    })
  
  }
  
  tmp<-as.vector(tmp)
  
  rownames(mat)<-tmp
  
  return(mat)
  
}

