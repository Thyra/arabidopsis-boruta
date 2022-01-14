# NOTE: This needs to be run with --max-ppsize=500000 because otherwise the pointer stack is too large and R will throw an error when running Boruta
library('Boruta')
library(randomForest)
library(stringr)
# You can use this seed to reproduce our exact results. If you change it, you might get a different set of genes
# or gene-trait connections but they should be just as valid candidates for drought stress response/recovery (see discusison)
set.seed(582477)

mothertable = read.csv("data/input/curated_mothertable20200515.lfs.txt.gz", sep="\t")
# Remove any genes are not differentially expressed between drought and WW treatment
mothertable = mothertable[mothertable$differentialinatleastoneconstrast == "TRUE",]
# We're only using the TPM values (Stephan used counts but Andrea suggested to go with TPM)
mothertable.tpm = mothertable[43:81]
colnames(mothertable.tpm) <- sub("_tpm$", "", colnames(mothertable.tpm)) # remove _count from the end of the columns
row.names(mothertable.tpm) <- mothertable$target_id
# Transpose it so that we have "HarvestX_repY" --> "Gene 1", "Gene 2", ...
mothertable.transposed = as.data.frame(t(mothertable.tpm))
mothertable.transposed$Mothertable.ID = row.names(mothertable.transposed)

mothertable.transposed$Treatment = as.factor(str_split_fixed(row.names(mothertable.transposed), "_", 3)[,2])
mothertable.transposed$DAS = as.integer(str_split_fixed(str_split_fixed(row.names(mothertable.transposed), "_", 3)[,1], "d", 2)[,2])

# Remove all DAS >= 30 then drop the column
mothertable.transposed = mothertable.transposed[mothertable.transposed$DAS < 30,]
mothertable.transposed = subset(mothertable.transposed, select = -c(DAS))

# Runs Boruta of specified trait ~ all genes and returns only the confirmed genes
runBoruta <- function(trait, dataframe) {
  
  # For speed/performance optimization, we're creating a y vector (trait) and an x dataframe containing all the
  # possible predictors (dataframe-trait)
  y = dataframe[, trait]
  x = dataframe[, !names(dataframe) %in% c(trait)]
  boruta_result = Boruta(x, y, pValue=0.01, maxRuns=800, ntree=1000)
  
  trait_importance = attStats(boruta_result)
  trait_importance = cbind(gene = rownames(trait_importance), trait_importance)
  trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
  trait_importance_roughfixed = cbind(gene = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
  return(trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed",])
}

write("i,Treatment,Prediction",file="data/results/cv/loocv_treatment_from_genes.csv")
for(i in 1:nrow(mothertable.transposed)) {
  data.train = mothertable.transposed[-i,]
  data.train.x = subset(data.train, select = -c(Treatment) )
  data.train.y = data.train$Treatment
  data.test  = mothertable.transposed[i,]
  
  boruta.result = runBoruta("Treatment", data.train)
  write.csv(boruta.result, paste0("data/results/cv/selected_genes/treatment_from_genes_", row.names(mothertable.transposed)[i], ".csv"), row.names=F)
  selected.features = boruta.result$gene
  
  # Only keep selected features
  data.train.x = subset(data.train.x, select = selected.features)
  
  # Train the model and predict
  model = randomForest(data.train.x, data.train.y)
  write(paste(row.names(mothertable.transposed)[i], as.character(data.test$Treatment), as.character(predict(model, newdata = subset(data.test, select = selected.features))), sep=","),file="data/results/cv/loocv_treatment_from_genes.csv",append = T)
}

