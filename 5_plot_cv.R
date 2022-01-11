library(caret)

pdf("plots/cv_results.pdf")

# Confusion Matrix for Trait CV
trait.cv = read.csv("data/results/cv/loocv_traits.csv", stringsAsFactors = T)
cm = confusionMatrix(data = trait.cv$Prediction, reference = trait.cv$Treatment)

# Snippet taken from https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package
draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'D', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'WW', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'D', cex=1.2, srt=90)
  text(140, 335, 'WW', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

draw_confusion_matrix(cm)

# Draw gene CVs
for(trait in list.files("data/results/cv/", pattern = "loocv_genes_*")) {
  d = read.csv(paste0("data/results/cv/", trait))
  plot(d$Prediction ~ d$Value, main = trait)
  abline(lm(d$Prediction ~ d$Value))
  cortest = cor.test(d$Value, d$Prediction)
  
  rmse = RMSE(pred = d$Prediction, obs = d$Value)
  nrmse = rmse/mean(d$Value) * 100
  mae = MAE(pred = d$Prediction, obs = d$Value)
  nmae = mae/mean(d$Value) * 100
  
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, 'RMSE', cex=1.2, font=2)
  text(10, 70, round(rmse, 3), cex=1.2)
  text(30, 85, 'NRMSE', cex=1.2, font=2)
  text(30, 70, paste0(round(nrmse, 1), " %"), cex=1.2)
  text(50, 85, 'MAE', cex=1.2, font=2)
  text(50, 70, round(mae, 3), cex=1.2)
  text(70, 85, 'NMAE', cex=1.2, font=2)
  text(70, 70, paste0(round(nmae, 1), " %"), cex=1.2)
  
  text(50, 35, 'R^2', cex=1.2, font=2)
  text(50, 20, round(cortest$estimate^2, 3), cex=1.2)
  
}

dev.off()