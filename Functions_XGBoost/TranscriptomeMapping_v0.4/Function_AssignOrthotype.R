#' @description This function assigns the orthotypes from the XGBoost prediction matrix
#' @param X Ordered normalized prediction matrix with columns are the training labels in the classifier
#' @param ann_data The annotation data
#' @param col_name_ref The cell type reference
#' @param percent_cutoff The cutoff used to assign the orthotypes.
#' @author Junqiang Wang
#' @export
AssignOrthotype<-function(X,
                          ann_data,
                          col_name_ref=NA,
                          col_name_add="orthotype",
                          percent_cutoff=50,
                          simplify_orthotype=TRUE,
                          simplify_orthotype_pattern="_"
                          ){
  
  require(reshape2)
  require(hash)
  
  
  if(is.na(col_name_ref)){stop("Error: please assign the col_name_ref parameter to select the column in the ann_data as the reference cell type!")}
  
  X0<-X
  
  # remove the NA column
  tmp<-match("NA", colnames(X))
  tmp<-setdiff(tmp, NA)
  
  if(length(tmp)>0){
    
    X<-X[, -tmp]
    
  }
  
  levels_label<-colnames(X)
  
  X<-melt(X)
  
  X$test_predlabels<-factor(X$test_predlabels, levels=levels_label)


  
  X<-X %>% 
    dplyr::group_by(test_predlabels) %>%
    dplyr::filter(value == max(value)) %>%
    dplyr::filter(value>percent_cutoff) %>%
    as.data.frame()
  

  # https://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters
  X[] <- lapply(X, as.character)
  
  # create the hash
  keys<-c(X$test_id, X$test_predlabels) 
  
  if(isTRUE(simplify_orthotype)){
    
    new.lables.1<-str_split(X$test_predlabels, pattern=simplify_orthotype_pattern, simplify = TRUE)[,1] 
    new.lables.2<-str_split(X$test_id, pattern=simplify_orthotype_pattern, simplify = TRUE)[,1] 
    values<-c(paste0(new.lables.1, ":", new.lables.2), paste0(new.lables.1, ":", new.lables.2))
    
  }else{
  values<-c(paste0(X$test_predlabels, ":", X$test_id), paste0(X$test_predlabels, ":", X$test_id))
  }
  
  hash_orthotype<-hash::hash(keys=keys, values=values)
  
  match_orders<-unique(values)
  
  # edit the ann_data
  orthotype<-ann_data[, match(col_name_ref, colnames(ann_data))]
  
  if(is.factor(orthotype)){
    
    orthotype<-as.character.factor(orthotype)
      
      }
  
  #
  tmp<-which(orthotype %in% keys(hash_orthotype))
  
  if(length(tmp)<1){stop("Error: no orthotype is assigned!")}
  
  values<-hash::values(x=hash_orthotype, keys=orthotype[tmp])
  
  values
  
  orthotype[tmp]<-values
  
  ann_data$orthotype<-orthotype
  
  if(!is.na(col_name_add)){
    
    tmp<-match("orthotype", colnames(ann_data))
    
    if(length(tmp)>0){
      
      colnames(ann_data)[tmp]<-col_name_add
      
      }
    
  }
  
  assigned_orthotype<-list("ann_data"=ann_data, "hash_orthotype"=hash_orthotype, "match_orders"= match_orders)
  
  return(assigned_orthotype)
  
}


