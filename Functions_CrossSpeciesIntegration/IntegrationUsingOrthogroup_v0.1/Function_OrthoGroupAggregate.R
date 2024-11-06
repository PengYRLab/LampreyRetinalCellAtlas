#'@param exp_mat expression matrix
#'@param ortho_group_table the orthoFinder output ortho_group_table
#'@param columns_used the columns used for generating the aggregated expression matrix
#'@param aggregate_method the method used. anyone in c(mean, sum)
#'@export
AggregateOrthoGroup<-function(exp_mat,
                              aggregate_method="sum",
                              ortho_group_table,
                              column_species_index,
                              pattern=","
                              ){
  
  
  require(pbapply)
  require(hash)
  
  
  #-------------------- create the hash table 
  
  ortho_group_table <- as.data.frame(ortho_group_table) 
  
  hash_orthogroup<-hash::hash(keys=rownames(ortho_group_table), values=ortho_group_table[,column_species_index])
  
  hash_gene_in_expression <- hash::hash(keys=rownames(exp_mat), values=1:nrow(exp_mat))
  
  #-------------------- generate the aggregated expression matrix

  message("start generating the aggregated expression matrix...")
  
  
  ortho_id<-sort(keys(hash_orthogroup))
  
  exp_mat_orthogroup_list<-pblapply(ortho_id, 
                                  FUN=CreateExpList,
                                  aggregate_method=aggregate_method,
                                  hash_orthogroup=hash_orthogroup,
                                  hash_gene_in_expression = hash_gene_in_expression,
                                  exp_mat=exp_mat,
                                  pattern=pattern
                                  )
  
  names(exp_mat_orthogroup_list) <- ortho_id
  
  # remove NA values in the list and zero length
  exp_mat_orthogroup_list <- exp_mat_orthogroup_list[!is.na(exp_mat_orthogroup_list)]
  
  tmp<-lapply(exp_mat_orthogroup_list, FUN=length) %>% unlist() %>% unique()
  
  if(tmp!=ncol(exp_mat)){stop("Error: error occurs when aggregate the expressions of the orthogroup!")}
  
  exp_mat_orthogroup <- matrix(unlist(exp_mat_orthogroup_list), ncol = ncol(exp_mat), byrow = TRUE)

  rownames(exp_mat_orthogroup)<-names(exp_mat_orthogroup_list)
  
  colnames(exp_mat_orthogroup)<-colnames(exp_mat)
  
  return(exp_mat_orthogroup)
   
}



CreateExpList<-function(
    x,
    aggregate_method,
    hash_orthogroup,
    hash_gene_in_expression,
    exp_mat,
    pattern
){
  
  # x, key in hash_orthogroup
  
  value_i<-unname(hash::values(hash_orthogroup, keys=x))
  
  elements_i<-str_replace_all(value_i, " ", "")
  
  elements_i <- as.vector(str_split(elements_i, pattern=pattern, simplify = TRUE))
  
  # check whether the elements_i exists in the hash table
  
  tmp<-has.key(key=elements_i, hash=hash_gene_in_expression)
  
  tmp<-tmp[tmp==TRUE]
  
  if(length(tmp)>0){
    
    elements_i<-names(tmp)
    
    tmp<-hash::values(x=hash_gene_in_expression, keys=elements_i)
    
    exp_i<-exp_mat[tmp,]
    
    
    if(length(tmp)>1){
      
      nrow_i<-nrow(exp_i)
      
      exp_i<-colSums(exp_i)
      
      if(aggregate_method=="mean"){
        
        exp_i<-exp_i/nrow_i
        
      }
      
    }
  
  }else{
  
    exp_i<-NA
  }
  
  return(exp_i)
  
}




