#' @description This function infers the DEGs for individual batches
#' @param de.method DEG inference method, can be "MAST" or "wilcox"
#' @param idents The idents for individual samples
#' @param batch The idents for batch
#' @author Junqiang Wang
#' @export
GetMarkersByBatches<-function(seu,
                              idents,
                              batch,
                              de.method="MAST",
                              pct.1=NULL,
                              pct.2=NULL
                              ){
  
  require(Seurat)
  require(dplyr)
  require(MAST)
  
  DefaultAssay(seu) <- "RNA"
  
  Idents(seu)<-idents
  
  seu.list <- SplitObject(seu, split.by =batch)
  
  
  marker.list<-vector(mode="list", length=length(seu.list))
  
  names(marker.list)<-names(seu.list)
  
  
  
  #
  for (i in 1:length(seu.list)){
    
    seu.i<-seu.list[[i]]
    
    Idents(seu.i)<-idents
    
    DefaultAssay(seu.i) <- "RNA"
    
  
    # use wilcox
    
    if(de.method=="wilcox"){
      
       seu.i<-SCTransform(seu.i, method = "glmGamPoi", return.only.var.genes = FALSE)
      
        markers <- FindAllMarkers(seu.i,
                                      test.use="wilcox",
                                      only.pos = TRUE )
        
    }

    if(de.method=="poisson"){
      
      markers <- FindAllMarkers(seu.i,
                                test.use="poisson",
                                slot="counts",
                                assay = "RNA",
                                latent.vars = c("orig.ident"),
                                only.pos = TRUE)
      
    }

  # use MAST 

    if(de.method=="MAST"){

        
        tmp<-seu.i@meta.data$orig.ident %>% unique() %>% length()
        
        if(tmp>1) {
        
        markers <- FindAllMarkers(seu.i,
                                      only.pos = TRUE,
                                      test.use="MAST",
                                      assay="RNA",
                                      slot="counts",
                                      latent.vars = c("orig.ident","nCount_RNA"))
        
        }else{
          
          markers <- FindAllMarkers(seu.i,
                                        only.pos = TRUE,
                                        test.use="MAST",
                                        assay="RNA",
                                        slot="counts",
                                        latent.vars = c("nCount_RNA"))
        
      
        }
    }
    
    
    
    # filtering
    
    if(!is.null(pct.1)){
      
      markers<- markers %>%
        dplyr::filter(pct.1>UQ(pct.1))
    }
    
    if(!is.null(pct.2)){
      
      markers<- markers %>%
        dplyr::filter(pct.2<UQ(pct.2))
    }
    
    
    head(markers)
    
    marker.list[[i]]<-markers
    
  }
  
  return(marker.list)
  
}

