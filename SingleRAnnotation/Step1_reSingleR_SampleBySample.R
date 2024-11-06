
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
library(SingleR)

set.seed(777)

#--------------
#cutoff
#--------------

my.nFeature_RNA<-600


#--------------------------
# Load the data 
#--------------------------


load(file="seu.ref.lamprey.s9.rda")

seu.ref<-seu


seu<-readRDS("../Lamprey_New.rds")


seu@meta.data$Class<-NA

seu@meta.data$Species<-"Lamprey"

DefaultAssay(seu)<-"RNA"


Idents(seu)<-"orig.ident"

seu<-subset(seu, idents="s9")


seu<-subset(seu, nFeature_RNA>my.nFeature_RNA)

seu.query<-seu



#-------------------------
# SingleR
#-------------------------

ref.data <- seu.ref@assays$RNA@counts %>% as.matrix() %>% Log2CPM()

query.data<-seu.query@assays$RNA@counts %>% as.matrix()%>% Log2CPM()

ref.labels<-seu.ref@meta.data$Class

seu<-seu.query

rm(seu.query)



#----------------
# single-cell level annotation
#------------------


predictions <- SingleR(test=query.data, 
                       ref=ref.data, 
                       labels=ref.labels)

table(predictions$labels)


save(predictions, file="predictions.s9.rda")




plotScoreHeatmap(predictions)


plotDeltaDistribution(predictions)

seu@meta.data$SingleR<-predictions$labels



seu@meta.data$Class<-seu@meta.data$SingleR %>% as.factor()


#---------------------
# redo clustering
#---------------------


seu<-SCTransform(seu, method = "glmGamPoi", return.only.var.genes = TRUE)

# 
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(seu, resolution = my.resolution)

seu@meta.data$seurat_clusters<-seu@meta.data$seurat_clusters %>% as.numeric() %>% as.factor()

p1 <- DimPlot(seu, reduction = "umap", 
              pt.size = 0.01,
              label=TRUE,
              repel=TRUE,
              group.by = "seurat_clusters")


p2 <- DimPlot(seu, reduction = "umap", 
              pt.size = 0.01,
              label=TRUE,
              repel=TRUE,
              group.by = "SingleR")


p1
ggsave("Clustering_SCT_umap_clusters_LampreySingleR_s9_r2.pdf", width=13, height = 10)

p2
ggsave("Clustering_SCT_umap_Class_LampreySingleR_s9_r2.pdf", width=13, height = 10)



seu@meta.data$Class<-seu@meta.data$SingleR %>% as.factor()


save(seu, file="seu.lamprey.s9.singleRannotated.r2.rda")







