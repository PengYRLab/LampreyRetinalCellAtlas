
# Plot 
PlotOrderedConfusionMatrix <- function(X,
                                   color_low="white", 
                                   # color_high="green4", 
                                   color_high= "#1b3865",
                                   max_size=8, 
                                   ylab_use="Prediction", 
                                   xlab_use="Train", 
                                   x_lab_rot=FALSE, 
                                   max_perc=100, 
                                   title_use = "Prediction matrix"){
  
  
  require(reshape2)
  
  X0<-X
  
  X = melt(X)
  colnames(X) = c("Prediction", "Train", "Percentage")
  
  X$Prediction<-factor(X$Prediction, levels=rev(rownames(X0)))
  
  X$Train<-factor(X$Train, levels=colnames(X0))
  
  
  p <- ggplot(X, aes(y = Prediction,  x = Train)) + 
    geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =color_low,   high = color_high, limits=c(0, 100 ))+
    scale_size(range = c(1, max_size), limits = c(0,max_perc))+   
    theme_bw() #+nogrid
  
  p = p + xlab(xlab_use) + ylab(ylab_use) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic")) + 
    ggtitle(title_use)+
    labs(color = "Match %") +
    labs(size="")+
    theme(text=element_text(size=12, face="italic"))
  # +
  #labs(size="Match %")
  
  if (x_lab_rot) {
    p <- p + theme(axis_text.x = element_text(angle = 90, vjust = 0.5))
  }
  
  print(p)
  
  p
  
}



