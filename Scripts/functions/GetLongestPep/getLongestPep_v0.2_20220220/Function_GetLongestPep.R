#
# Input
# The input files are NCBI protein.faa and cds_from_genomic.fna
# 
# The function
# The function retrieve the pep id info and the transcript length info 
# The function retrieve the pep id info and the transcript id info
# The function merge the info together
# The function get the longest transcript id
# The function write the longest peptides with the transcript ID to a new file

# Read the files 
# 

# This function convert the XP id to Loc id


GetLongestPep<-function(filein_cds,
                        filein_pep,
                        fileout_pep.longest
                        ){
#
require(dplyr)
require(seqinr)
require(stringr)
require(tibble)

#
  
  dnaseq<- read.fasta(file = filein_cds)
  
  ann <- getAnnot(dnaseq)
  
  ann %>% unlist() -> ann
  
  # get the cds id and gene loc id
  ann<-str_extract_all(string=ann, pattern="\\[gene=\\S+\\]|\\[protein_id=\\S+\\]", simplify = TRUE)
  
  #ann<-str_extract_all(string=ann, pattern="\\[gene=[\\w|\\d]+\\]|\\[protein_id=\\S+\\]", simplify = TRUE)
  
  ann<-gsub("\\[gene=|\\]", "", ann)
  
  ann<-gsub("\\[protein_id=|\\]", "", ann) %>% as.data.frame()
  
  colnames(ann)<-c("loc_id", "xp_id")
  
  ann.genes<-ann
  
  
  
  #-----------------------------
  # process the peptide files
  #-----------------------------
  
  pepseq<- read.fasta(file = filein_pep, seqtype = "AA")
  
  ann <- getAnnot(pepseq)
  
  ann<-str_extract_all(string=ann, pattern="XP_\\S+", simplify = TRUE)
  
  ann.pep.id<-ann %>% as.vector()
  
  pep.length<-getLength(pepseq)
  
  length(pep.length)
  
  ann<-cbind(ann, pep.length)
  
  colnames(ann)<-c("xp_id", "pep.length")
  
  ann.pep<-ann %>% as.data.frame()
  
  #-----------------------
  # Update the annotation
  #-----------------------
  
  tmp.c<-intersect(ann.genes$xp_id, ann.pep$xp_id) %>% sort()
  
  tmp<-match(tmp.c, ann.genes$xp_id)
  
  ann.genes<-ann.genes[tmp,]
  
  tmp<-match(tmp.c, ann.pep$xp_id)
  
  ann.pep<-ann.pep[tmp,]
  
  ann<-cbind(ann.genes, ann.pep$pep.length) %>% as.data.frame()
  
  colnames(ann)<-c("loc_id", "xp_id", "pep.length")
  
  message("\nthe number of total trancripts are: ", nrow(ann))
  
  message("\nwrite the annotation file \n")
  
  write.csv(ann, file="ann.file.csv", quote=FALSE)
  
  #--------------------------
  # Get the longest gene's id
  #--------------------------
  
  ann<-as.tibble(ann)
  
  ann %>%
    group_by(loc_id) %>%
    slice_max(n = 2, order_by = pep.length)
  
  ann %>%
    group_by(loc_id) %>%
    top_n(n = 1, wt = pep.length) -> pep.longest
  
  #pep.longest <-pep.longest %>% as.data.frame()
  
  pep.longest %>% distinct(loc_id, .keep_all=TRUE) ->pep.longest
  
  dim(pep.longest)
  
  message("\nthe number of the longest transcripts are: ", dim(pep.longest)[1])
  
  #--------------------------
  # get the longest peptide
  #--------------------------
  
  tmp<-match(pep.longest$xp_id, ann.pep.id)
  
  seq.pep.longest<-getSequence(pepseq)[tmp]
  
  seq.pep.longest.ann<-getAnnot(pepseq)[tmp]
  
  
  
  for(i in 1:length(seq.pep.longest.ann)){
    
    tmp<- pep.longest$loc_id[i]
    
    seq.pep.longest.ann[[i]]<-tmp
    
  }
  
  
  
  
  # for(i in 1:length(pepseq.longest)){
  # 
  #   tmp<-paste0(">", pep.longest$loc_id[i])
  # 
  #   attr(pepseq.longest[[i]], 'Annot')<-tmp
  # 
  # }
  # 
  
  #-------------------------------------
  # write the sequence to the file
  #-------------------------------------
  
  message("\nwrite the output peptide file")
  
  write.fasta(sequences=seq.pep.longest, names=seq.pep.longest.ann,  file.out=fileout_pep.longest)
  
}





