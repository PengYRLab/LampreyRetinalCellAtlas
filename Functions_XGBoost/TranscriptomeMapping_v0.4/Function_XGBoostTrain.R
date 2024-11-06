
# The function is modified from a previously publication
# Peng YR, Shekhar K, Yan W, Herrmann D, Sappington A, Bryman GS, van Zyl T, Do MTH, Regev A, Sanes JR. Molecular Classification and Comparative Taxonomics of Foveal and Peripheral Cells in Primate Retina. Cell. 2019 Feb 21;176(5):1222-1237.e22. doi: 10.1016/j.cell.2019.01.004. Epub 2019 Jan 31. PMID: 30712875; PMCID: PMC6424338.


XGBoostTrain <- function(train_Data, 
                          train_labels = NULL, 
                          do.scale=FALSE,
                          scale.mean = NULL, 
                          scale.sd = NULL,
                          max.cells.per.ident = 300, 
                          train.frac = 0.66,
                          width_confusion_figure=10,
                          height_confuison_figure=10
                          ){
  
  if (!is.factor(train_labels)){
    train_labels = factor(train_labels)
  }
  
  if (length(train_labels) != ncol(train_Data)) stop("Error: The size of the training labels is not equal to sample size of the training data!")
  
  if (do.scale){
    if (is.null(scale.mean) & is.null(scale.sd)){
      scale.mean = rowMeans(train_Data)
      scale.sd = apply(train_Data, 1, sd)
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.sd))
    } else {
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.sd))
    }
  }
  
  # Training set vs. validation set
  training.set = c(); 
  validation.set=c()
  training.label = c(); 
  validation.label=c();
  print(paste0("Using either ", train.frac*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training, whichever is smaller"))
  
  for (cc in c(1:length(levels(train_labels)))){
    i = levels(train_labels)[cc]
    cells.in.clust = names(train_labels)[train_labels == i];
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp)
    validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(cc-1,length(train.temp)))
    validation.label = c(validation.label, rep(cc-1, length(validation.temp)))
  }
  
  train_matrix <- xgb.DMatrix(data = t(train_Data[,training.set]), label=training.label)
  
  validation_matrix <- xgb.DMatrix(data = t(train_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,
                     "max_depth"= 6, 
                     subsample = 0.66)
  
  nround  <- 200 # number of XGBoost rounds
  
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  
  validation_prediction <- matrix(validation_pred, 
                                  nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)
  
  A = table(validation.label, valid_predlabels)
  
 rownames(A) = levels(train_labels)
 
 tmp<-colnames(A) %>% as.numeric()
 
 colnames(A) = levels(train_labels)[tmp]; 
  
 PlotConfusionMatrix(A, order="Row", xlab.use = "True", ylab.use = "Predicted")
  
 ggsave("confusion.train.pdf", width=width_confusion_figure, height = height_confuison_figure)

  to.return = list()
  to.return$bst_model = bst_model
  to.return$scale_mean = scale.mean
  to.return$scale_sd = scale.sd
  
  return(to.return)
}
