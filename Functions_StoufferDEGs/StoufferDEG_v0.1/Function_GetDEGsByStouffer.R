

GetDEGsByStouffer<-function(seu.marker.list,
                            remove.na=TRUE
                            ){

  
  require(tidyverse)
  require(purrr)
  require(dplyr)
  
  seu.marker.list.merged<-purrr::reduce(seu.marker.list, full_join, by = c("cluster", "gene"))
  
  
  if(remove.na==TRUE){
    
    seu.marker.list.merged<-seu.marker.list.merged %>%
      na.omit()
    
  }
  
  
  dat<-seu.marker.list.merged
  
  
  # 
  dat.p<-dat[, grep("p_val\\.", colnames(dat))]
  
  dat.p[is.na(dat.p)]<-0.5
  
  z.stouffer<-StoufferP2Z(dat.p)
  
  p.stouffer<-pnorm(z.stouffer, lower.tail = FALSE)
  
  p.stouffer.fdr<-p.adjust(p.stouffer, method="fdr")
  

  #
  dat.log2fc<-dat[, grep("avg_log2FC\\.", colnames(dat))]
  
  dat.log2fc[is.na(dat.log2fc)]<-0
  
  
  avg_log2FC.mean<-apply(dat.log2fc, 1, function(x){
    x<-sum(x)/length(x)
    
    
    
  })
  
  
  seu.marker.list.merged.updated<-bind_cols(list(seu.marker.list.merged, z.stouffer, p.stouffer, p.stouffer.fdr, avg_log2FC.mean))
  
  colnames(seu.marker.list.merged.updated)<-c(colnames(seu.marker.list.merged), "z.stouffer", "p.stouffer", "p.stouffer.fdr", "avg_log2FC.mean")
  
  
  seu.marker.list.merged.updated<-seu.marker.list.merged.updated %>% 
    as.data.frame() %>%
    group_by(cluster) %>%
    arrange(desc(z.stouffer), .by_group = TRUE)
  
  
  
  return(seu.marker.list.merged.updated)
  
}




