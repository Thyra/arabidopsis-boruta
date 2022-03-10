# This script performs a Boruta trait selection on the wide-shaped table in intermediary_data/1745AJ_Phenotyping_nafixed.csv
# Results are saved unfixed and roughFixed (with or without Tentative decisions)

library('Boruta')
library('randomForest')
library(caret)
library(ggplot2)
set.seed(17022019)

data.cleaned = read.csv("data/intermediary/1745AJ_Phenotyping_cleaned.csv")
# Cut off recovery phase
# data.cleaned = data.cleaned[data.cleaned$DAS < 31,]
# Reshape data to wide format
data = reshape(data.cleaned, idvar='Plant.ID', timevar='DAS', v.names=setdiff(names(data.cleaned), c("Plant.ID", "DAS", "Treatment", "Time")), direction='wide')
# Drop the plant ID
data = data[,setdiff(names(data), c("Plant.ID"))]
data$Treatment = as.factor(data$Treatment)
# Drop the plants with missing values (harvested I assume)
data = data[rowSums(is.na(data)) == 0,]

boruta_result = Boruta(Treatment ~ ., data=data, maxRuns=1500, doTrace=2)
trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
trait_importance_roughfixed = cbind(trait = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
write.csv(trait_importance_roughfixed, "data/results/trait_selection_result_roughFixed.csv", row.names=F)

# ---- LOOCV-CV -------
write("i,Treatment,Prediction",file="data/results/cv/loocv_traits.csv")
for(i in 1:nrow(data)) {
  # Split data into training set and test row
  data.train   = data[-i,]
  data.train.x = subset(data.train, select = -c(Treatment) )
  data.train.y = data.train$Treatment
  data.test    = data[i,]
  
  # Select features on training set
  boruta_result = Boruta(data.train.x, data.train.y, maxRuns=1500)
  trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
  trait_importance_roughfixed = cbind(trait = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
  selected.features = trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed",]$trait
  
  # Only keep selected features
  data.train.x = subset(data.train.x, select = selected.features)
  
  # Train the model and predict
  model = randomForest(data.train.x, data.train.y)
  write(paste(i, as.character(data.test$Treatment), as.character(predict(model, newdata = subset(data.test, select = selected.features))), sep=","),file="data/results/cv/loocv_traits.csv",append = T)
}

# Plot it!
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
pdf("plots/trait_cv_confusion_matrix.pdf")
draw_confusion_matrix(cm)
dev.off()

# ---- Importance plots ----
# Plot importance over time per category
trait.categories = read.csv("data/input/trait_selection_list.csv")
trait.categories$Variable.ID = make.names(trait.categories$Variable.ID)
trait.categories$Category = trait.categories$imaging.modality
trait.categories[trait.categories$trait.type %in% c("architectural", "biomass-related"),]$Category = "architectural"

selected.points = trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed", c("trait", "medianImp")]
selected.points$trait.name = sub("\\.\\d+$", "", selected.points$trait)
selected.points$DAS = as.numeric(gsub("^.*\\.", "", selected.points$trait))

merged = merge(selected.points[c("trait.name", "DAS", "medianImp")], trait.categories[c("Variable.ID", "Category")], by.x="trait.name", by.y="Variable.ID")
agg = aggregate(medianImp ~ Category + DAS, data=merged, FUN=sum)
pdf("plots/aggregated_importances_per_category_and_day.pdf")
p1 <- ggplot(agg, aes(DAS, medianImp, fill = Category)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  geom_segment(x = 30.5, xend=30.5, y=0, yend=45, linetype="dashed", color = "black", size=.7)
print(p1)
dev.off()
