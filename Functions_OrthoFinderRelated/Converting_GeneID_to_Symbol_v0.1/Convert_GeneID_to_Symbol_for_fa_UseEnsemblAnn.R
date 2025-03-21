
# load the library
require(dplyr)
require(tibble)


create_hash<-function(ensembl_annotation){
  
  require(hash)
  require(dplyr)
  
  nrow(ensembl_annotation)
  
  ensembl_annotation<- ensembl_annotation %>%
                        dplyr::filter(Gene.name !="")
  
  nrow(ensembl_annotation)
  
  hash_1<-hash(keys=ensembl_annotation$Gene.stable.ID.version, values=ensembl_annotation$Gene.name)
  
  hash_1
  
}



ENSG_2_Symbol<-function(ann_each, hash){
  
  require(dplyr)
  #require(seqinr)
  #require(stringr)
  require(tibble)
  require(hash)


  ann_each<-str_extract(string=ann_each, pattern="ENSG.*")
  
  if(has.key(key=ann_each, hash)==TRUE){
    
   ann_new<-values(hash, keys=ann_each)

  # ann_new<-paste0(">", ann_new)
   
  }else{
    
    ann_new<-ann_each
    
  }
  
  return(ann_new)
  
}




GeneID_2_Symbol_UseEnsemblAnn<-function(filein_ensembl_annotation_path,
                          filein_pep_path,
                          fileout="out.converted.fa"
                          ){
  
  require(dplyr)
  require(seqinr)
  require(stringr)
  require(tibble)
  require(hash)
  
  
  ensembl_annotation<-read.csv(file=filein_ensembl_annotation_path, sep="\t")
  
  ensembl_annotation
  
  pepseq<- read.fasta(file = filein_pep_path, seqtype = "AA")
  
  pepseq
  
  ann <- getAnnot(pepseq)
  
  # create_hash
  hash_1<- create_hash(ensembl_annotation)
  
 #ann update 
ann_new<-lapply(ann, ENSG_2_Symbol, hash=hash_1)

seq.pep<-getSequence(pepseq)

seq.pep.ann<-ann_new

write.fasta(sequences=seq.pep, names=seq.pep.ann,  file.out=fileout, open="w")

}

