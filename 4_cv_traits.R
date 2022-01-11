# Perform leave-one-out cross-validation for the important trait selection

library('Boruta')
library('randomForest')
set.seed(17022019)

data.nafixed = read.csv("data/intermediary/1745AJ_Phenotyping_nafixed.csv")
data.nafixed$Treatment = as.factor(data.nafixed$Treatment)

write("i,Treatment,Prediction",file="data/results/cv/loocv_traits.csv")
for(i in 1:nrow(data.nafixed)) {
  # Split data into training set and test row
  data.train   = data.nafixed[-i,]
  data.train.x = subset(data.train, select = -c(Treatment) )
  data.train.y = data.train$Treatment
  data.test    = data.nafixed[i,]
  
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