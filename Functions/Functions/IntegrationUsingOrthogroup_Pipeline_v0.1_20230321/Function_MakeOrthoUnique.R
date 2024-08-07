
#' @description  preprocess the table and remove genes with the same names in each row for each species
#' @param  x the orthogroup table (e.g. orthofinder output)
#' @param pattern_gene_collapse_input the pattern that gene collapse used in the input
#' @param pattern_gene_collapse_output the pattern that gene collapse used in the output
#' @author Junqiang Wang
#' @export
MakeOrthoUnique<-function(x, 
                          pattern_gene_collapse_input=",", 
                          pattern_gene_collapse_output=","){
    
    require(dplyr)
    require(stringr)
    
   # x<-table[2,2]
   # change in version 3: remove the blank/empty string
    x<-str_replace_all(x, " ", "")
    
    x<-str_split(x, pattern=pattern_gene_collapse_input, simplify=TRUE) %>% as.vector() %>% unique()
    
    x<-paste0(as.character(x), collapse=pattern_gene_collapse_output)
    
    return(x)
    
  }
  
  # table<-apply(table, MARGIN=c(1,2), FUN=MakeOrthoUnique)
  # 
  # 
  # #
  # table<-table %>%
  #   as.data.frame() %>%
  #   filter(Lamprey!="" & Mouse!="") %>%
  #   distinct(Lamprey, .keep_all= TRUE) %>%
  #   distinct(Mouse, .keep_all= TRUE)
  # 
  # table
  
  
  
  # get 1-many matched genes (no bigger than 3 for each species)
  # count the orthologues 
  CountOrtho<-function(x, 
                       pattern=","){
    
    x<-str_replace_all(x, " ", "")
    
    x<-str_split(x, pattern=pattern, simplify=TRUE) %>% as.vector() %>% unique() %>% length() 
    
    return(x)
    
  }

  
  
  