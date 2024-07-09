# This script performs a Boruta trait selection on cleaned table in intermediary_data/1745AJ_Phenotyping_nafixed.csv
# Results are saved unfixed and roughFixed (with or without Tentative decisions)
# The script then performs a leave-one-out cross-validation to see how well the selected traits can predict the treatment class

library('Boruta')
library('randomForest')
set.seed(17022019)

data.nafixed = read.csv("data/intermediary/1745AJ_Phenotyping_nafixed.csv")

# Transform to wide format
data.nafixed = reshape(data.nafixed, idvar='Plant.ID', timevar='DAS', v.names=setdiff(names(data.nafixed), c("Plant.ID", "DAS", "Treatment")), direction='wide', sep="_")
data.nafixed = data.nafixed[-c(1)] # Drop Plant.ID
data.nafixed$Treatment = as.factor(data.nafixed$Treatment)

boruta_result = Boruta(Treatment ~ ., data=data.nafixed, maxRuns=1500, doTrace=2)
trait_importance = attStats(boruta_result)
trait_importance = cbind(trait = rownames(trait_importance), trait_importance)
# Some results from playing with this:
#  - order of most important traits varies a bit each time this is run but median trait importance stays quite stable (< 1% difference between runs)
write.csv(trait_importance, "data/results/trait_selection_result.csv", row.names=F)

trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
trait_importance_roughfixed = cbind(trait = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
write.csv(trait_importance_roughfixed, "data/results/trait_selection_result_roughFixed.csv", row.names=F)

# Perform Leave-one-out cross-validation
selected.features = trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed",]$trait
data.nafixed = data.nafixed[,c("Treatment", selected.features)]
write("i,Treatment,Prediction",file="data/results/cv/loocv_traits.csv")
for(i in 1:nrow(data.nafixed)) {
  # Split data into training set and test row
  data.train   = data.nafixed[-i,]
  data.train.x = subset(data.train, select = -c(Treatment) )
  data.train.y = data.train$Treatment
  data.test    = data.nafixed[i,]
  
  # Train the model and predict
  model = randomForest(data.train.x, data.train.y)
  write(paste(i, as.character(data.test$Treatment), as.character(predict(model, newdata = subset(data.test, select = selected.features))), sep=","),file="data/results/cv/loocv_traits.csv",append = T)
}
