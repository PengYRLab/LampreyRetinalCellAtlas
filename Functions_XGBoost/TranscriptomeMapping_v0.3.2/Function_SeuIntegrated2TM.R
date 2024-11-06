#' This function creates the output data for TM from the seu
#' The output data is a list which includes the mat of the variables genes and degs, the train id, the variable genes, the degs
#' The degs are found by stoufer integration of different batches
#' @author  Junqiang Wang
#' @export

SeuIntegrated2TM<-function(
                 seu,
                 assay="SCT",
                 slot="scale.data",
                 seu.idents="species",
                 seu.train.idents="mouse",
                 seu.test.idents="lamprey",
                 label.train="seurat_clusters",
                 label.test="seurat_clusters",
                 common.features=NULL
){
  
  #
  require(Seurat)
  
  to.return = list()
  
  DefaultAssay(seu)<-assay
  
  Idents(seu)<- seu.idents
  
  seu.list<-SplitObject(seu, split.by =  seu.idents)
  tmp1<-match(seu.train.idents, names(seu.list))
  tmp2<-match(seu.test.idents, names(seu.list))
  
  seu.train<-seu.list[[tmp1]]
  seu.test<-seu.list[[tmp2]]
  
  #
  tmp<-match(label.train, colnames(seu.train@meta.data))
  train_id<-seu.train@meta.data[,tmp] %>% as.vector()
  names(train_id)<-rownames(seu.train@meta.data)
  train_id<-train_id %>% as.factor()
  
  #
  tmp<-match(label.test, colnames(seu.test@meta.data))
  test_id<-seu.test@meta.data[,tmp] %>% as.vector()
  test_id<-test_id %>% as.factor()
  names(test_id)<-rownames(seu.test@meta.data)
  

  if(is.null(common.features)){common.features<-rownames(seu)}
  
  #
  to.return$var.genes_train<-common.features
  
  #
  DefaultAssay(seu.train)<-assay
  DefaultAssay(seu.test)<-assay
  
  # fetch data in a slot
 
  to.return$train_mat<- GetAssayData(object = seu.train, assay = assay, slot = slot) %>% as.matrix()
  
  to.return$train_id<-train_id
  
  to.return$var.genes_test<-common.features
  
  to.return$test_mat<- GetAssayData(object = seu.test, assay = assay, slot = slot) %>% as.matrix()
  
  to.return$test_id<-test_id
  
  xg_data<-to.return
  
  return(xg_data)
  
}

