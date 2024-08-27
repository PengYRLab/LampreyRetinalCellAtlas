
#-------------
# RGC clustering analysis
#--------------


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
library(harmony)
library(cowplot)

options(ggrepel.max.overlaps = Inf)

set.seed(777)

#------------
# Load the data
#------------

load("../seu.lamprey.rgc.r2 _finalObject_archived.rda")

#----------------------
# Check the quality
#---------------------

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size =0, ncol = 3)


plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


Idents(seu)<-"seurat_clusters"


#-------------------------
#  Do Harmony
# if error, check as.factor()
#--------------------------

DefaultAssay(seu)<-"RNA"

seu<-SCTransform(seu, method = "glmGamPoi", return.only.var.genes = TRUE)

seu <- RunPCA(seu, verbose = FALSE, features = VariableFeatures(object = seu))


seu@meta.data$orig.ident<-seu@meta.data$orig.ident %>% as.vector() %>% as.factor()


seu <- seu %>% 
  RunHarmony(group.by.vars=c("orig.ident"), plot_convergence = TRUE)


#-----------------
# Plot the data
#-----------------

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = seu, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = seu, features = "harmony_1", group.by = "Class_SingleR", pt.size = .1)
p3 <- VlnPlot(object = seu, features = "harmony_1", group.by = "seurat_clusters", pt.size = .1)


plot_grid(p1)

plot_grid(p2)

plot_grid(p3)



#-----------------
# Dim reduction 
#-----------------

seu <- seu %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 1.33) %>% 
  identity()


seu@meta.data$seurat_clusters<- seu@meta.data$seurat_clusters %>% as.numeric() %>% as.vector() %>% as.factor()

seu@meta.data$samples<-seu@meta.data$orig.ident

#--------------
# The
#--------------

pt.size<-0.5

options(repr.plot.height = 4, repr.plot.width = 6)
p1<-DimPlot(seu, reduction = "umap", group.by="seurat_clusters", label = TRUE, pt.size = .1)
p1
ggsave(file="Harmony_Lamprey_Dimplot_clusters_r4.pdf", width=7, height=7)


p2<-DimPlot(seu, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = .1)
p2

ggsave(file="Harmony_Lamprey_Dimplot_origIdent_r4.pdf", width=7, height=7)


#-------------------
# Plot the quality
#-------------------


VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size =0, ncol = 1, group.by = "seurat_clusters" )

ggsave(file="Harmony_Lamprey_QP_VlnPlot_r4.pdf", width=10, height=8)


#------------------
# Cellular composition
#-------------------

#plot as proportion or percentage of cluster
ggplot(seu@meta.data, aes(x=seurat_clusters, fill=orig.ident)) +
  geom_bar(position = "fill")+ylab("Percentage")+xlab("Clusters") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1, size=6))

ggsave(file="Harmony_Lamprey_CellCompostionBydonor_r4.pdf", width=8.15, height=2.6)





