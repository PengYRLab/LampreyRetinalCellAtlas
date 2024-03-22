

#' @description  This function converts the orthogroupIDs to gene symbols
#' @param mat Expression matrix
#' @param orthgroup_table The orthogroup table used
#' @param reference_species The species in the orthogroup table for which the symbols will used. Can be multiple species.
#' @author Junqiang Wang
#' @export
OrthogroupID2GeneSymbol<-function(mat,
                                  orthogroup_table,
                                  reference_species
                                  ){
  
  # orthogroup_table<-ortho_group_table_reference
  # orthogroup_ids<-rownames(mat)
  # reference_species<-c("Lamprey", "Zebrafish")
  
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

