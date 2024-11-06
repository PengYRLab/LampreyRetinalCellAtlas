


rm(list=ls())


library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(RColorBrewer)

library(stringr)

library(ComplexHeatmap)
library(varhandle)
library(ggpubr)

library(MumuscRNASeq)

library(SingleR)

options(ggrepel.max.overlaps = Inf)

set.seed(777)


#-------------
# Load the data
#-------------

load("seu.lamprey.s9.singleRannotated.r2.rda")


seu@meta.data$Class_SingleR<-seu@meta.data$SingleR



ann<-seu@meta.data

ann$seurat_clusters %>% unique() %>% sort()->clusters


for(i in clusters){
  
  tmp.id<-which(ann$seurat_clusters==i)
  
  tmp.cell_type<-ann$Class_SingleR[tmp.id] %>% 
    table() %>%             
    as.data.frame() %>% 
    arrange(desc(Freq))
  
  print(tmp.cell_type)
  
  ann$Class_SingleR[tmp.id]<- tmp.cell_type[1,1] %>% as.character()
  
}

ann$Class_SingleR <-ann$Class_SingleR %>% as.factor()

seu@meta.data<-ann 


#-----------------
# Plot the annotation
#-----------------


p1 <- DimPlot(seu, reduction = "umap", 
              pt.size = 0.01,
              label=TRUE,
              repel=TRUE,
              group.by = "seurat_clusters")


p2 <- DimPlot(seu, reduction = "umap", 
              pt.size = 0.01,
              label=TRUE,
              repel=TRUE,
              group.by = "Class_SingleR")



p1
ggsave("Clustering_SCT_umap_clusters_LampreySingleR_s9_r2_ClusterLevel.pdf", width=13, height = 10)

p2
ggsave("Clustering_SCT_umap_Class_SingleR_LampreySingleR_s9_r2_ClusterLevel.pdf", width=13, height = 10)



#----------------
# save the data
#----------------

save(seu, file="seu.lamprey.s9.singleRannotated.r2.ClusterLevel.rda")












