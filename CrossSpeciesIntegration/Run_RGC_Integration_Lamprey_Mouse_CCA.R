


# rm(list=ls())

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
library(Seurat)
library(cowplot)
library(dplyr)

set.seed(777)


#-----------------------
# Create the directory
#-----------------------

pwd<-getwd()
pwd

dir<-paste0(pwd, "/", "figures_dimplot")

dir

if(dir.exists(dir)==FALSE){
  
  dir.create(dir)
  
}

setwd(dir)

getwd()

#-------------------------------
#  The data lamprey
#-------------------------------

load("seu.lamprey.RGC.types.rda")

Idents(seu)<-"orig.ident"

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu@meta.data$orig.ident %>% table()

seu@meta.data$sample<-seu@meta.data$orig.ident

seu@meta.data$subtype %>% unique() %>% sort()


seu@meta.data$subtype<-paste0("RGC", seu@meta.data$subtype)

p1<-DimPlot(seu, reduction = "umap",
            group.by="subtype",
            label = TRUE,
            repel=TRUE,
            cols=col_vector,
            raster=FALSE,
            label.box=FALSE)
p1
p1<-p1+gg.theme+labs(title="")
p1


seu@meta.data$subtype %>% table()


Idents(seu)<-"subtype"
seu<-subset(seu, downsample=100)
seu@meta.data$subtype %>% table()
Idents(seu)<-"subtype"

#

DefaultAssay(seu)<-"RNA"

seu@meta.data$orig.ident<-seu@meta.data$orig.ident %>% as.vector() %>% as.factor()

seu@meta.data$species<-"lamprey"

seu@meta.data$cell_type<-seu@meta.data$subtype

seu@meta.data$subtype2<-seu@meta.data$subtype

seu.lamprey<-seu

rm(seu)


#----------------------------
#  The data mouse
#----------------------------

load("rgc.published.45celltypes.rda")


Idents(seu)<-"orig.ident"

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu@meta.data$sample<-seu@meta.data$orig.ident

seu@meta.data$subtype %>% unique() %>% sort()

p1<-DimPlot(seu, reduction = "umap",
            group.by="subtype",
            label = TRUE,
            repel=TRUE,
            cols=col_vector,
            raster=FALSE,
            label.box=FALSE)
p1
p1<-p1+gg.theme+labs(title="")
p1


seu@meta.data$orig.ident %>% table()



Idents(seu)<-"orig.ident"

seu<-subset(seu, idents=c("aRGC1", "aRGC2", "aRGC3", "aRGC4"))

#seu<-subset(seu, idents=c("MouseACS6"))

Idents(seu)<-"subtype"

seu<-subset(seu, downsample=100)

seu

seu@meta.data$subtype %>% table()


Idents(seu)<-"subtype"




p1<-DimPlot(seu, reduction = "umap",
            group.by="subtype",
            label = TRUE,
            repel=TRUE,
            cols=col_vector,
            raster=FALSE,
            label.box=FALSE)
p1
p1<-p1+gg.theme+labs(title="")
p1

seu@meta.data$subtype %>% table()


DefaultAssay(seu)<-"RNA"

seu@meta.data$orig.ident<-seu@meta.data$orig.ident %>% as.vector() %>% as.factor()

seu@meta.data$species<-"mouse"

seu@meta.data$cell_type<-seu@meta.data$subtype

seu@meta.data$subtype2<-paste0(seu@meta.data$subtype, "(mm)")


#
seu.meta<-seu@meta.data

seu.count<-seu@assays$RNA@counts


#------------- Gene name convertor

source("GeneNameConvertorLocal_Mat_mouse2Human.R")

seu.count<-GeneNameConvertorLocal_Mat_Mouse2Human(mat=seu.count)

#
seu<-CreateSeuratObject(counts=seu.count)

rm(seu.count)

seu@meta.data<-seu.meta

seu.mouse<-seu

#-----------------------------
# Merge these two objects
#-----------------------------

seu<-merge(x=seu.lamprey, y=seu.mouse)

#
seu@meta.data$species<- seu@meta.data$species %>% as.factor()
seu@meta.data$orig.ident<- seu@meta.data$orig.ident %>% as.factor()

#
Idents(seu)<-"orig.ident"
seu.list<-SplitObject(seu, split.by = "orig.ident")
seu.list

#
#seu.list <- lapply(X = seu.list, FUN = SCTransform)

seu.list <- lapply(X = seu.list, FUN = SCTransform, method = "glmGamPoi", vst.flavor="v1", conserve.memory =TRUE)

features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)

features1 <- SelectIntegrationFeatures(object.list = seu.list[c(1,2)], nfeatures = 3000)

features2 <- SelectIntegrationFeatures(object.list = seu.list[c(3,4,5,6)], nfeatures = 3000)

features12<-intersect(features1, features2)

features

intersect(features, features12)

features<-features12

#


seu.list <- PrepSCTIntegration(object.list = seu.list, anchor.features = features)

seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#
anchors <- FindIntegrationAnchors(object.list = seu.list, 
                                  # reference=ref.id,  
                                  normalization.method = "SCT", 
                                  scale=TRUE,
                                  anchor.features = features,
                                  reduction = "cca",
                                  dims = 1:50)

# rm(seu.list)

seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)

# rm(anchors)

#DefaultAssay(seu)<-"integrated"

seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(seu, resolution = 0.5)

seu@meta.data$seurat_clusters<-seu@meta.data$seurat_clusters %>% as.numeric() %>% as.factor()

seu@meta.data$integrated_clusters<-seu@meta.data$seurat_clusters

#-------------------------------------
# Do plot 
#-------------------------------------


#---
cols2<-c("blue", "green4")

p1<-DimPlot(seu, reduction = "umap", 
            group.by = "species", 
            order=rev(rownames(seu@meta.data)),
            label = FALSE, 
            cols=cols2,
            repel=TRUE)

p1


p1<-p1+gg.theme3

p1

file<- paste0("Dimplot_species", ".pdf")

ggsave(file=file, width=5.5, height=5)



#------

p1<-DimPlot(seu, reduction = "umap", 
            group.by = "species", 
            order=rev(rownames(seu@meta.data)),
            label = FALSE, 
            cols=cols2,
            repel=TRUE)

p1

p1<-p1+gg.theme3+NoLegend()

p1

file<- paste0("Dimplot_species_Nolgend_width5.5_height5", ".pdf")

ggsave(file=file, width=5.5, height=5)


#-----------------------------
# subtype 
#-----------------------------

Idents(seu)<-"species"

seu.mouse<-subset(seu, idents="mouse")

seu.lamprey<-subset(seu, idents="lamprey")


p3<-DimPlot(seu.mouse, reduction = "umap", 
            group.by="subtype", 
            cols=col_vector,
            label = TRUE, repel=TRUE)

p3

p3<-p3+gg.theme3

p3


file<- paste0("Dimplot_subtypes_mouse", ".pdf")

ggsave(file=file, width=6, height=5)


#-----------

p3<-DimPlot(seu.mouse, reduction = "umap", 
            group.by="subtype", 
            cols=col_vector,
            label = TRUE, repel=TRUE)

p3

p3<-p3+gg.theme3+NoLegend()

p3

file<- paste0("Dimplot_subtypes_mouse_Nolgend_width5.5_height5", ".pdf")

ggsave(file=file, width=5.5, height=5)


#---------
# lamprey
#--------

p3<-DimPlot(seu.lamprey, reduction = "umap", 
            group.by="subtype", 
            cols=col_vector,
            label = TRUE, repel=TRUE)

p3

p3<-p3+gg.theme3

p3


file<- paste0("Dimplot_subtypes_lamprey", ".pdf")

ggsave(file=file, width=6, height=5)


#-----------

p3<-DimPlot(seu.lamprey, reduction = "umap", 
            group.by="subtype", 
            cols=col_vector,
            label = TRUE, repel=TRUE)

p3

p3<-p3+gg.theme3+NoLegend()

p3

file<- paste0("Dimplot_subtypes_lamprey_Nolgend_width5.5_height5", ".pdf")

ggsave(file=file, width=5.5, height=5)



#-----------------------
# Plot the quality
#-----------------------


VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size =0, ncol = 1, group.by = "seurat_clusters" )

ggsave(file="QP_VlnPlot.pdf", width=5, height=7)


#-------------------------
# The END!
#-------------------------





