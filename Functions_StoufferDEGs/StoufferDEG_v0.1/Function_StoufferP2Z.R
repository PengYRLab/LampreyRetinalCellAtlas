#' @description This function converts multiple P values to a Z score using the Stouffer's method
#' @param dat Input 
#' @author Junqiang Wang
#' @export
StoufferP2Z<-function(dat){
  
  minNum<- .Machine$double.xmin
  maxNum<- .Machine$double.xmax
  
  z.minNum<-qnorm(minNum)
  z.maxNum<-qnorm(minNum, lower.tail = FALSE)
  
  
  dat<-as.matrix(dat)
  
  dat<-qnorm(dat, lower.tail = FALSE)
  
  dat[dat=="Inf"]<-z.maxNum
  
  dat[dat=="-Inf"]<-z.minNum
  
  z.stouffer<-apply(dat, 1, function(x){
    x<-sum(x)/sqrt(length(x))
  })
  
  return(z.stouffer)
  
}



