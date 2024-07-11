library("Boruta")
library(randomForest)
set.seed(17022019)

mothertable = read.csv("./data/input/curated_mothertable20200515.lfs.txt.gz", sep = "\t")
# Remove any genes are not differentially expressed
mothertable = mothertable[mothertable$differentialinatleastoneconstrast == "TRUE", ]
row.names(mothertable) = mothertable$target_id
# We're only using the TPMs
mothertable = mothertable[43:81]
colnames(mothertable) <- sub("_tpm$", "", colnames(mothertable))
# Transpose it so that we have "HarvestX_repY" --> "Gene 1", "Gene 2", ...
mothertable.transposed = as.data.frame(t(mothertable))
mothertable.transposed$Mothertable.ID = row.names(mothertable.transposed)

# Load the traits we want to look at
important.traits = read.csv("data/results/trait_selection_result_roughFixed.csv")
important.traits = sub("\\_\\d+$", "", important.traits[important.traits$decision == "Confirmed", ]$trait)
important.traits = unique(important.traits)

# This tells us which plants belong to which harvest and rep
mothertable.mapping = read.csv("data/input/MothertableToPlantID_mapping.csv")

# Load Phenotpying data
pheno = read.csv("data/intermediary/1745AJ_Phenotyping_nafixed.csv")

# We only need the important traits, plant name, treatment and DAS
pheno = pheno[c(names(pheno)[c(1, 2, 3)], important.traits)]
pheno$Day = paste("d", pheno$DAS, sep = "")

# Merge by Plant.ID and Day (and Treatment, but that doesn't matter)
# We lose a number of samples here because we only have comparable phenotyping
# data between DAS 20 and 35. Even some samples at DAS 20 are lost because
# they were apparently not phenotyped with the new configuration before sampling.
pheno.merged = merge(mothertable.mapping, pheno)

options("expressions" = 50000)

# Runs Boruta of specified trait ~ all genes and returns only the confirmed genes
runBoruta <- function(trait, dataframe) {
  # For speed/performance optimization, we're creating a y vector (trait) and an x dataframe containing all the
  # possible predictors (dataframe-trait)
  y = dataframe[, trait]
  x = dataframe[, !names(dataframe) %in% c(trait)]
  boruta_result = Boruta(x, y, pValue = 0.01, maxRuns = 800, doTrace = 2, ntree = 1000)

  trait_importance = attStats(boruta_result)
  trait_importance = cbind(gene = rownames(trait_importance), trait_importance)
  trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
  trait_importance_roughfixed = cbind(gene = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
  return(trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed", ])
}

important.genes = data.frame(trait = character(), gene = character(), medianImp = double())
for (trait in important.traits) {
  print(paste("Analyzing trait", trait))

  # Samples were taken from two plants so we need to average the phenotypic data
  plant.means = aggregate(as.formula(paste(trait, "~ Mothertable.ID")), data = pheno.merged, FUN = mean)
  # Join Trait and transposed mothertable
  mothertable.forBoruta = merge(plant.means, mothertable.transposed)[, -1] # remove Mothertable.ID so it's all numbers

  # Make one Boruta run and another one and keep only the genes that were confirmed in both runs.
  # Keep going making more runs until no more genes are being thrown out.
  confirmed.genes.allRuns = runBoruta(trait, mothertable.forBoruta)
  # Throw out all genes from Boruta input that have been rejected
  mothertable.forBoruta = mothertable.forBoruta[c(trait, as.character(confirmed.genes.allRuns$gene))]
  confirmed.genes.currentRun = runBoruta(trait, mothertable.forBoruta)

  # While there are genes in the allRuns dataframe that are not confirmed in the current run, do another run!
  while(length(setdiff(unique(as.character(confirmed.genes.allRuns$gene)), unique(as.character(confirmed.genes.currentRun$gene)))) > 0) {
    genes.intersect = intersect(unique(as.character(confirmed.genes.allRuns$gene)), unique(as.character(confirmed.genes.currentRun$gene)))
    # Remove genes from both datasets that were not present in both
    confirmed.genes.currentRun = confirmed.genes.currentRun[confirmed.genes.currentRun$gene %in% genes.intersect, ]
    confirmed.genes.allRuns = confirmed.genes.allRuns[confirmed.genes.allRuns$gene %in% genes.intersect, ]
    confirmed.genes.allRuns = rbind(confirmed.genes.allRuns, confirmed.genes.currentRun)
    # Do another run!
    print(paste(" Doing another run now (", length(genes.intersect), " genes left)", sep = ""))
    # Throw out all genes from Boruta input that have been rejected again
    mothertable.forBoruta = mothertable.forBoruta[c(trait, as.character(confirmed.genes.allRuns$gene))]
    confirmed.genes.currentRun = runBoruta(trait, mothertable.forBoruta)
  }

  selected_genes = data.frame(aggregate(medianImp ~ gene, data = confirmed.genes.allRuns, FUN = median))
  selected_genes$trait = trait
  important.genes = rbind(important.genes, selected_genes)

  # Do Leave-one-out cross-validation
  mothertable.forBoruta = mothertable.forBoruta[c(trait, as.character(selected_genes$gene))]
  write("i,Value,Prediction",file=paste0("data/results/cv/loocv_genes_",trait,".csv"))
  for (i in 1:nrow(mothertable.forBoruta)) {
    data.train = mothertable.forBoruta[-i, ]
    data.test  = mothertable.forBoruta[i, ]

    data.train.x = subset(data.train, select = selected_genes$gene)
    data.train.y = data.train[, trait]
    model = randomForest(data.train.x, data.train.y, ntree = 1000)
    write(paste(i, data.test[, trait], predict(model, newdata = subset(data.test, select = selected_genes$gene)), sep=","),file=paste0("data/results/cv/loocv_genes_",trait,".csv"),append = T)
  }
}

important.genes = important.genes[order(important.genes$trait, important.genes$gene), ]
write.csv(important.genes, "data/results/gene_selection_result_roughFixed.csv", row.names = FALSE)
