
#' @description This function assign the orthotypes from the XGBoost prediction matrix
#' @param X Ordered normalized prediction matrix with columns are the training labels in the classifier
#' @param ann_data The annotation data
#' @param col_name_ref The cell type reference also appeared in X
#' @param percent_cutoff The cutoff used to call the orthotypes. All passed elements will be used.
#' @author Junqiang Wang
#' @export
AssignOrthotype<-function(X,
                          ann_data,
                          col_name_ref=NA,
                          col_name_add="orthotype",
                          percent_cutoff=50
                          ){
  
  require(reshape2)
  require(hash)
  
  
  if(is.na(col_name_ref)){stop("Error: please assign the col_name_ref parameter to select the column in the ann_data as the reference cell type!")}
  
  X0<-X
  
  # remove the NA column
  tmp<-match("NA", colnames(X))
  
  if(length(tmp)>0){
    
    X<-X[, -tmp]
    
  }
  
  X<-melt(X)
  
  # chose the biggest values for each test_predlable group
  
  X<-X %>% 
    dplyr::group_by(test_predlabels) %>%
    dplyr::filter(value == max(value)) %>%
    dplyr::filter(value>percent_cutoff) %>%
    as.data.frame()
  

  # https://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters
  X[] <- lapply(X, as.character)
  
  # create the hash
  keys<-c(X$test_id, X$test_predlabels)
  values<-c(paste0(X$test_predlabels, ":", X$test_id), paste0(X$test_predlabels, ":", X$test_id))
  
  hash_orthotype<-hash::hash(keys=keys, values=values)
  
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
  
  assigned_orthotype<-list("ann_data"=ann_data, "hash_orthotype"=hash_orthotype)
  
  return(assigned_orthotype)
  
}


