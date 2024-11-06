
library(dplyr)
library(seqinr)
library(stringr)
library(tibble)

#-----------------
# process the cds
#------------------

dnaseq<- read.fasta(file = "cds_from_genomic.fna")

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

pepseq<- read.fasta(file = "protein.faa", seqtype = "AA")

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

ann<-ann %>% group_by(loc_id)

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

message("the number of the longest transcripts are: ", dim(pep.longest)[1])

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


#-------------------------------------
# write the sequence to the file
#-------------------------------------

write.fasta(sequences=seq.pep.longest, names=seq.pep.longest.ann,  file.out="longest.lamprey.pep.faa")








