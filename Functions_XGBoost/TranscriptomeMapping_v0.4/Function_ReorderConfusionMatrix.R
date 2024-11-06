#' @description This function reorders the rows and columns of the confusion matrix
#' @param row_scale Scale the confusion matrix by row
#' @param make_diagonal Place the best match on the diagonal line for optimal positioning.
#' @author Junqiang Wang
#' @export 
ReorderConfusionMatrix <- function(X,
                                   row_scale=TRUE, 
                                   reorder_row=TRUE,
                                   reorder_column=TRUE,
                                   make_diagonal=FALSE,
                                   percent_matrix_partation=0.5){
  
  message("Only performing row scaling by default, please transpose the confusion matrix if column scale should be performed!")
  
  require(reshape2)
  
  if (row_scale){ 
    
  if(max(X)<1){stop("Error: the confusion matrix is already scaled!")}
    
   X = t(scale(t(X), center=FALSE, scale=colSums(t(X))))
  
   X = X*100 
   
   X<-round(X, digits=2)
   
   }
  
  
  # get the sub-matrix without NA column
  
  X0<-X
  
  tmp<-match("NA", colnames(X))
  tmp<-setdiff(tmp, NA)
  
  if(length(tmp)>0){
    
    X<-X[, -tmp]
  }
  
  
  if(make_diagonal){ 
    reorder_row<-TRUE
    reorder_column<-TRUE
    percent_matrix_partation<-NA
                    
  }
  
  
  #-------------------------------Make diagonal
  
    if (reorder_column==TRUE){
      
   # 1. arrange by column (decreasing order)
      tmp<-apply(X, 2, max)
      
      factor_levels_col<-colnames(X)[order(tmp, decreasing=TRUE)]
      
  
      X<-X[ , match(factor_levels_col, colnames(X))]
      
      
    }
  
  
  if (reorder_row==TRUE){
    # 2. arrange by row (decreasing order)
      tmp<-apply(X, 1, max)
      
      factor_levels_row<-rownames(X)[order(tmp, decreasing=TRUE)]
      
      # factor_levels <- c(factor_levels, setdiff(colnames(X), factor_levels)) 
      
      X<-X[match(factor_levels_row, rownames(X)), ]
      
  }
  
  
  
  if (make_diagonal==TRUE){
    
      # 3. rearrange to make diagonal
    
      tmp<-apply(X, 2, which.max)
      
      factor_levels_row <-unique(rownames(X)[tmp])
    
      factor_levels_row <- c(factor_levels_row, setdiff(rownames(X), factor_levels_row))
      
  }
  
  
  

# if cutoff is set
  if(!is.na(percent_matrix_partation)){
    
  tmp<-apply(X, 2, FUN = function(x){which(x>percent_matrix_partation)})
  
  tmp<-unique(unlist(tmp))
  
  factor_levels_row <-unique(rownames(X)[tmp])
  
  factor_levels_row <- c(factor_levels_row, setdiff(rownames(X), factor_levels_row))
  
  }
  

   factor_levels_col <- c(factor_levels_col, setdiff(colnames(X0), factor_levels_col))
  
   factor_levels_row <- c(factor_levels_row, setdiff(rownames(X0), factor_levels_row))
   

# rearrange X by factor levels

  X<-X0[match(factor_levels_row, rownames(X0)), match(factor_levels_col, colnames(X0))]
  
  X_return<-X
  
}





