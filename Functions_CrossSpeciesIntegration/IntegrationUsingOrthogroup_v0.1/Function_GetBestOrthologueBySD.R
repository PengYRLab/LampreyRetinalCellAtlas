#'@param exp_mat expression matrix
#'@param ortho_group_table the orthoFinder output ortho_group_table
#'@param columns_used the columns used for generating the aggregrated expression matrix
#'@param aggregate_method the method used. anyone in c(mean, sum)
#'@export
GetBestOrthologuesBySD<-function(
                              ortho_group_table=NULL,
                               exp_mat,
                              species=NULL,
                              pattern=","
                              ){

  
  require(pbapply)
  require(hash)
  
  #-------------------- create the hash table 
  
  ortho_group_table <- as.data.frame(ortho_group_table) 
  
  col_id<-match(species, colnames(ortho_group_table))
  
  orthologues<-ortho_group_table[,col_id ] %>% as.vector()
  
  #-------------------- generate the aggregated expression matrix

  message("starting get the best matched orthologues by SD ...")
  
  orthologues.best.match<-pbsapply(orthologues, simplify = TRUE, FUN=CompareExpSD, exp_mat=exp_mat, pattern=pattern) %>% unname()
  
  ortho_group_table[,col_id]<- orthologues.best.match
  
  return(ortho_group_table)
}



CompareExpSD<-function(
    x,
    exp_mat,
    pattern=","
){
  
  
  x<-str_replace_all(x, " ", "")
  x<-str_split(x,pattern=pattern, simplify=TRUE) %>% as.vector()

  
  if(length(x)>1){
    
  sd<-sapply(x, simplify=TRUE, FUN=function(x, exp_mat){
    exp<-exp_mat[match(x[1], rownames(exp_mat)), ]
    sd<-sd(exp)
    
  }, exp_mat=exp_mat)
  
  x<-which.max(sd) %>% names()
        
  }
  
  return(x)
}



