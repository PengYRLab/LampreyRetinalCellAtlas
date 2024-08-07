

rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(RColorBrewer)

library(SingleR)
library(stringr)
library(biomaRt)

library(MumuscRNASeq)
library(ggpubr)
library(Seurat)
library(patchwork)


library("readxl")

library(hash)

set.seed(777)

#----------------------

load("/Users/junqiangwang/Documents/MumuAtPengLab/Data_And_Analysis/LampreyS9S10_StringTie_MTadded_Analysis/seu.lamprey.filtered.rda")

ann <- read_excel("Lampreygene3_wTargetIDupdated_mumu.xlsx")

dim(ann)

ann<- ann %>% dplyr::distinct(stringtie_gene_id, .keep_all = TRUE)

dim(ann)

ann<- ann %>% dplyr::distinct(Target_gene_id, .keep_all = TRUE)

dim(ann)


#hash.ann<-hash(keys=ann$stringtie_gene_id, values=ann$Target_gene_id)
hash.ann<-hash(keys=ann$gene_name, values=ann$Target_gene_id)

genes.orig<-rownames(seu@assays$RNA@counts)

genes.annotated<-genes.orig


for (i in 1:length(genes.orig)){
  
  gi<-genes.orig[i]

  if(has.key(gi, hash.ann)==TRUE){

    genes.annotated[i]<-values(hash.ann, keys=gi, simplify=TRUE)

    }
}


genes.annotated

rownames(seu@assays$RNA@counts)<-genes.annotated


seu.counts<-seu@assays$RNA@counts %>% as.matrix()

unique.row.names<-rownames(seu.counts) %>% unique() %>% sort()



unique.row.names<-setdiff(unique.row.names, "SDS-1")


seu.counts<-seu.counts[match(unique.row.names, rownames(seu.counts)), ]

seu.ann<-seu@meta.data

seu<-CreateSeuratObject(counts = seu.counts, project = "lamprey")
seu@meta.data$orig.ident<-seu.ann$orig.ident
seu@meta.data$Sample<-seu.ann$Sample


save(seu, file="seu.lamprey.filtered.geneAnnotated.rda")


# Check the data

id1<-grep("LOC", genes.annotated)
id2<-grep("MSTR", genes.annotated)

tmp<-genes.annotated[-c(id1, id2)] %>% length()

tmp







