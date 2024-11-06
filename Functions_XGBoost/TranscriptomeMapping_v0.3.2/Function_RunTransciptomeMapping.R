#' @param xg_data xgboost input data
#' @param p_cutoff p cutoff used for assign the NA labels 
#' @author Junqiang wang
#' @export
RunTransciptomeMapping<-function(xg_data,
                                 p_cutoff=0.5
                                 ){

require(ggpubr)
require(xgboost)
require(reshape2)
require(xgboost)

set.seed(777)

# 
train_data <- xg_data$train_mat
test_data <- xg_data$test_mat
train_id <- xg_data$train_id
test_id <- xg_data$test_id

bst_model <- XGBoostTrain(train_Data=train_data, 
                          train_labels = train_id, 
                          do.scale = FALSE)

save(bst_model, file="bst_model.rda")

load("bst_model.rda")

# Use trained model to predict on test data
numberOfClasses <- length(levels(train_id))

test_data <- xgb.DMatrix(t(test_data))

test_pred <- predict(bst_model$bst_model, 
                     newdata = test_data)

test_prediction <- matrix(test_pred, 
                          nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses)


# Find best label for each cell

test_pred_margins <- apply(test_prediction,2,max)

test_predlabels <- apply(test_prediction,2,which.max)

test_predlabels <- levels(train_id)[test_predlabels]


#-----------------------------------------------------------------------------------------
# set the cutoff and assign the label to NA if p.softmax < p.cutoff
#-----------------------------------------------------------------------------------------

id.na<-which(test_pred_margins<p_cutoff)


if(length(id.na)<1){
  
  
  test_predlabels <- factor(test_predlabels, levels=levels(train_id))
  
  C <- table(test_id, test_predlabels)
  
  C
  
  }else{

    test_predlabels[id.na]<-"NA"

    test_predlabels <- factor(test_predlabels, levels=c(levels(train_id), "NA"))

    C <- table(test_id, test_predlabels)
 
    C
 
  }

#---------------------------------
# add non-used label columns
#----------------------------------

id_unmathed<-setdiff(unique(train_id), colnames(C))

#
if(length(id_unmathed)>0){

Cc<-matrix(0, nrow=nrow(C), ncol=length(id_unmathed))

colnames(Cc)<-id_unmathed

C2<-cbind(C, Cc)

C2

tmp<-c(levels(train_id), "NA")

tm<-intersect(tmp, colnames(C2))

C2<-C2[,tmp]

C<-C2

}

return(C)

}


