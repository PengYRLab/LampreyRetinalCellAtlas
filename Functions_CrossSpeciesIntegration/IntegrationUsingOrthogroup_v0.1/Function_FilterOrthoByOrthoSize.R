
#' @description This function filters the orthogroup by the orthgroup size
#' @param x the orthogroup table
#' @param size_max the maximal size passed the filtering
#' @param pattern_gene_collapse_input the pattern that gene collapse used in the input
#' @author  Junqiang Wang
#' @export
  FilterOrthoByOrthoSize<-function(x, 
                                   size_max=5,
                                   pattern_gene_collapse_input=","
                                   ){
    
    
     x<-as.data.frame(x)
     
    count_ortho <- apply(x, MARGIN=c(1,2), FUN=CountOrtho, pattern=pattern_gene_collapse_input)
    
    tmp <- apply(count_ortho, 1, max) 
    
    tmp <- which(tmp>size_max)
    
    if(length(tmp)>0){
      
      x<-x[-tmp, ]
    }
    
    return(x)
    
  }
  
  
  CountOrtho<-function(x, 
                       pattern=NULL){
    
    x<-str_replace_all(x, " ", "")
    
    x<-str_split(x, pattern=pattern, simplify=TRUE) %>% as.vector() %>% unique() %>% length() 
    
    return(x)
    
  }
  
  
