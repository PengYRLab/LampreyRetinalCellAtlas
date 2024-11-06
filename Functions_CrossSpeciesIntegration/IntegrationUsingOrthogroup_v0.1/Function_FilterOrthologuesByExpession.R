#' @description  This function performs filtering
#' @param x a column vector of the input genes
#' @param y a vector used for filtering the input
#' @param pattern_collapse_x_input the pattern that gene collapse used in the input
#' @param pattern_collapse_x_output the pattern that gene collapse used in the output
#' @author Junqiang Wang
#' @export
FilterOrthologuesByExpression<-function(x,
                                        y,
                                        pattern_collapse_x_input=",",
                                        pattern_collapse_x_output=NULL
                                        ){
  
  
  if(is.null(pattern_collapse_x_output)){
    
   pattern_collapse_x_output<-pattern_collapse_x_input}
  
   x<-str_replace_all(x, " ", "")
  
   x<-lapply(x, FUN=str_split, pattern=pattern_collapse_x_input, simplify=TRUE)
  
   x<-lapply(x, intersect, y=y)
   
   x<-lapply(x, function(x) if (length(x) == 0)  "NA" else x)
   
   x<-lapply(x, function(x){paste0(as.character(x), collapse=pattern_collapse_x_output)})
   
   x<-do.call(rbind, x)
   
   return(x)
  
}




