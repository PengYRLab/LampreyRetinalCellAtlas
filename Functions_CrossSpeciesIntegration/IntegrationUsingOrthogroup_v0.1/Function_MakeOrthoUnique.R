#' @description  preprocess the table and remove genes with the same names in each row for each species
#' @param  x the orthogroup table (e.g. orthofinder output)
#' @param pattern_gene_collapse_input the concatenated pattern
#' @param pattern_gene_collapse_output the concatenated pattern
#' @author Junqiang Wang
#' @export
MakeOrthoUnique<-function(x, 
                          pattern_gene_collapse_input=",", 
                          pattern_gene_collapse_output=","){
    
    require(dplyr)
    require(stringr)
    

    x<-str_replace_all(x, " ", "")
    
    x<-str_split(x, pattern=pattern_gene_collapse_input, simplify=TRUE) %>% as.vector() %>% unique()
    
    x<-paste0(as.character(x), collapse=pattern_gene_collapse_output)
    
    return(x)
    
  }
  
  CountOrtho<-function(x, 
                       pattern=","){
    
    x<-str_replace_all(x, " ", "")
    
    x<-str_split(x, pattern=pattern, simplify=TRUE) %>% as.vector() %>% unique() %>% length() 
    
    return(x)
    
  }

  
  
  