# Entropy calculation from the prediction matrix
#' @param x The prediction matrix
#' @author  Junqiang Wang
#' @export
Entropy<-function(x,
                              trained_labels_in_column=TRUE){
  
  if(!isTRUE(trained_labels_in_column)){
    message("Rotate the matrix, assuming the predicted labels are the columns.")
    x<-t(x)
  }
  
  entropy<-vector()
  
  entropy<- apply(x, 1, function(p) -sum(p * log(p+1e-9)))
  
  message("Entropy of predictions:", entropy)
  
  entropy.sum<-sum(entropy)
  
  message("Total entropies:", entropy.sum)
  
  return(entropy.sum)
  
}



# Normalized Entropy calculation from the prediction matrix
#' @param x x is prediction matrix
#' @author  Junqiang Wang
#' @export
NormalizedEntropy<-function(x,
                            trained_labels_in_column=TRUE){
  
  if(!isTRUE(trained_labels_in_column)){
    message("Rotate the matrix, assuming the predicted labels are the columns.")
    x<-t(x)
  }
  
  entropy<-vector()
  
  entropy<- apply(x, 1, function(p) -sum(p * log(p+1e-9)))
  
  message("Entropy of predictions:", entropy)
  
  entropy.max<-apply(x, 1, function(x) log(length(x)))
  
  
  entropy.norm<-sum(entropy)/sum(entropy.max)
  
  entropy.norm<-round(entropy.norm, digits = 2)
  
  message("Normalized Entropy:", entropy.norm)
  
  return(entropy.norm)
  
}


# Run Examples
# m<-c(0.9, 0.05, 0.05, 0.33, 0.33, 0.34,0.25, 0.25, 0.5,0.5, 0.5, 0.0)
# x<-matrix(m, nrow=3)
# x<-t(x)
# x

# entropy<-NormalizedEntropy(x=x, predicted_label_by_row=TRUE); entropy
#
# 


