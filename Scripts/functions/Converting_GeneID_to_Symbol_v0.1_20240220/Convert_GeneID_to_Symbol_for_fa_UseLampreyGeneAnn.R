
#input gene_ann
# input fa
# This function converts the GeneID to gene symbol in the fa file if gene symbol exists 

# load the library
#require(dplyr)
#require(seqinr)
#require(stringr)
#require(tibble)

# input arguments:
#filein_gene_annotation_path<-"/Users/mumu/Documents/research/LampreyRetina/analysis/gtf_symbol_update/Lampreygene3_wTargetIDupdated_mumu.xlsx" 
#filein_pep_path<-"/Users/mumu/Documents/research/LampreyRetina/analysis/OrthoFinder/proteinSeq/longestPeptide/Pm_lamprey.fa"


create_hash<-function(gene_annotation){
  
  require(hash)
  require(dplyr)
  
  nrow(gene_annotation)
  
  gene_annotation<- gene_annotation %>%
                        dplyr::filter(Target_gene_id !="")
  
  nrow(gene_annotation)
  
  hash_1<-hash(keys=gene_annotation$ref_gene_id, values=gene_annotation$Target_gene_id)
  
  hash_1
  
}



#function convert the gene to gene symbl for individual item
LOC_2_Symbol<-function(ann_each, hash){
  
  require(dplyr)
  #require(seqinr)
  #require(stringr)
  require(tibble)
  require(hash)
  
  #ann_each<-ann[[1]] ; ann_each

   ann_each<-substring(text=ann_each, first=2)
  
  if(has.key(key=ann_each, hash)==TRUE){
    
   ann_new<-values(hash, keys=ann_each)

  # ann_new<-paste0(">", ann_new)
   
  }else{
    
    ann_new<-ann_each
    
  }
  
  return(ann_new)
  
}




GeneID_2_Symbol_UseGeneAnn<-function(filein_gene_annotation_path,
                          filein_pep_path,
                          fileout="out.converted.fa"
                          ){
  
  require(dplyr)
  require(seqinr)
  require(stringr)
  require(tibble)
  require(hash)
  require(readxl)
  
  
  gene_annotation <- read_excel(filein_gene_annotation_path)
  
  
  gene_annotation<- gene_annotation %>% 
     dplyr::distinct(stringtie_gene_id, .keep_all = TRUE) %>%
     dplyr::distinct(Target_gene_id, .keep_all = TRUE)
  
  dim(gene_annotation)
  
  pepseq<- read.fasta(file = filein_pep_path, seqtype = "AA")
  
  pepseq
  
  ann <- getAnnot(pepseq)
  
  # create_hash
  hash_1<- create_hash(gene_annotation)
  
 #ann update 
ann_new<-lapply(ann, LOC_2_Symbol, hash=hash_1)

seq.pep<-getSequence(pepseq)

seq.pep.ann<-ann_new

write.fasta(sequences=seq.pep, names=seq.pep.ann,  file.out=fileout, open="w")

}

# Run_example
# GeneID_2_Symbol_UseGeneAnn(filein_gene_annotation=filein_gene_annotation, filein_pep=filein_pep)

