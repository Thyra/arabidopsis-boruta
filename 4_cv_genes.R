# NOTE: This needs to be run with --max-ppsize=500000 because otherwise the pointer stack is too large and R will throw an error when running Boruta
library('Boruta')
library(randomForest)
# You can use this seed to reproduce our exact results. If you change it, you might get a different set of genes
# or gene-trait connections but they should be just as valid candidates for drought stress response/recovery (see discusison)
set.seed(582477)

mothertable = read.csv("data/input/curated_mothertable20200515.lfs.txt.gz", sep="\t")
# Remove any genes are not differentially expressed between drought and WW treatment
mothertable = mothertable[mothertable$differentialinatleastoneconstrast == "TRUE",]
# We're only using the count values for now, not the TPMs (originally we used TPMs but it appears that Stephan used counts)
mothertable.counts = mothertable[4:42]
colnames(mothertable.counts) <- sub("_count$", "", colnames(mothertable.counts)) # remove _count from the end of the columns
row.names(mothertable.counts) <- mothertable$target_id
# Transpose it so that we have "HarvestX_repY" --> "Gene 1", "Gene 2", ...
mothertable.transposed = as.data.frame(t(mothertable.counts))
mothertable.transposed$Mothertable.ID = row.names(mothertable.transposed)

write.csv(mothertable.transposed, "data/intermediary/mothertable_transposed.csv")

# Load the traits we want to look at
important.traits = read.csv("data/results/trait_selection_result_roughFixed.csv")
important.traits = sub("\\.\\d+$", "", important.traits[important.traits$decision == "Confirmed",]$trait)
important.traits = unique(important.traits)

# This tells us which plants belong to which harvest and rep
mothertable.mapping = read.csv("data/input/MothertableToPlantID_mapping.csv")

# Ignore the first 18 plants/lines because that's where Stephan weirdly had NA's in his script.
# We'll see if we can reproduce his results this way.
mothertable.mapping = mothertable.mapping[19:nrow(mothertable.mapping),]

# Load Phenotpying data
pheno = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")
pheno$Plant.ID = as.factor(pheno$Plant.ID)
pheno$Treatment = as.factor(pheno$Treatment)

# -------- Do some outlier correction --------------
# We only need the important traits, plant name, treatment and DAS
pheno = pheno[c(names(pheno)[c(1,6,10)], important.traits)]

# Reshape data to wide format
data = reshape(pheno, idvar='Plant.ID', timevar='DAS', v.names=setdiff(names(pheno), c("Plant.ID", "DAS", "Treatment")), direction='wide', sep="_")
# Drop the plant ID
# data = data[,setdiff(names(data), c("Plant.ID"))]
# Get rid of columns that only contain NAs (again)
# data = data[colSums(!is.na(data)) > 0] # We can't do that because then reshape -> long gets confused

data.drought = data[data$Treatment == "D",]
data.well_watered = data[data$Treatment == "WW",]

# Outlier detection -> replace with NA's which are then replaced with medians in following step.
medians = apply(data.drought[,setdiff(names(data.drought), c("Plant.ID", "Treatment"))], 2, FUN=median, na.rm=T)
sds = apply(data.drought[,setdiff(names(data.drought), c("Plant.ID", "Treatment"))], 2, FUN=sd, na.rm=T)
# This is how you would work with MADs instead of SDs:
# mads = apply(data.well_watered[,setdiff(names(data.well_watered), c("Plant.ID", "Treatment"))], 2, FUN=mad, na.rm=T)
for(traitname in names(data.drought[,setdiff(names(data.drought), c("Plant.ID", "Treatment"))])) {
  print(traitname)
  if(nrow(data.drought[which(data.drought[traitname] > medians[traitname]+2*sds[traitname] | data.drought[traitname] < medians[traitname]-2*sds[traitname]), ]) > 0)
    data.drought[which(data.drought[traitname] > medians[traitname]+2*sds[traitname] | data.drought[traitname] < medians[traitname]-2*sds[traitname]), ][traitname] <- NA
}

# Remove NA's, currently just by replacing with column median but might want to use a more sophisticated
# method like rfImpute that weighs by proximity.
data.drought = na.roughfix(data.drought)

# Outlier detection -> replace with NA's which are then replaced with medians in following step.
medians = apply(data.well_watered[,setdiff(names(data.well_watered), c("Plant.ID", "Treatment"))], 2, FUN=median, na.rm=T)
sds = apply(data.well_watered[,setdiff(names(data.well_watered), c("Plant.ID", "Treatment"))], 2, FUN=sd, na.rm=T)
for(traitname in names(data.well_watered[,setdiff(names(data.well_watered), c("Plant.ID", "Treatment"))])) {
  print(traitname)
  if(nrow(data.well_watered[which(data.well_watered[traitname] > medians[traitname]+2*sds[traitname] | data.well_watered[traitname] < medians[traitname]-2*sds[traitname]), ]) > 0)
    data.well_watered[which(data.well_watered[traitname] > medians[traitname]+2*sds[traitname] | data.well_watered[traitname] < medians[traitname]-2*sds[traitname]), ][traitname] <- NA
}

# Remove NA's, currently just by replacing with column median but might want to use a more sophisticated
# method like rfImpute that weighs by proximity.
data.well_watered = na.roughfix(data.well_watered)

# Combine drought and well-watered data again
data.nafixed = rbind(data.drought, data.well_watered)
pheno.fixed = reshape(data.nafixed, direction="long", varying = names(data.nafixed)[3:722], timevar="DAS", times=9:44, sep="_", idvar="Plant.ID")
# Get rid of columns that contain NAs (because they were dropped for one treatment)
# data.nafixed = data.nafixed[colSums(!is.na(data.nafixed)) == nrow(data.nafixed)]

# --------------------

pheno.fixed$Day = paste("d", pheno.fixed$DAS, sep="")

# Merge by Plant.ID and Day (and Treatment, but that doesn't matter)
pheno.merged = merge(mothertable.mapping, pheno.fixed)

options("expressions"=50000)

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

for(trait in important.traits) {
  print(paste("Analyzing trait", trait))
  write("i,Value,Prediction",file=paste0("data/results/cv/loocv_genes_",trait,".csv"))
  
  plant.means = aggregate(as.formula(paste(trait, "~ Mothertable.ID")), data=pheno.merged, FUN=mean)
  
  # Join Trait and transposed mothertable
  mothertable.forBoruta = merge(plant.means, mothertable.transposed)[,-1] # remove Mothertable.ID so it's all numbers
  
  for(i in 1:nrow(mothertable.forBoruta)) {
    data.train = mothertable.forBoruta[-i,]
    data.test  = mothertable.forBoruta[i,]
    
    boruta.result = runBoruta(trait, data.train)
    write.csv(boruta.result, paste0("data/results/cv/selected_genes/", trait, "_", i, ".csv"), row.names=F)
    selected.features = boruta.result$gene
    
    # Build model and predict
    data.train.x = subset(data.train, select = selected.features)
    data.train.y = data.train[, trait]
    model = randomForest(data.train.x, data.train.y)
    write(paste(i, data.test[, trait], predict(model, newdata = subset(data.test, select = selected.features)), sep=","),file=paste0("data/results/cv/loocv_genes_",trait,".csv"),append = T)
    
  }
}
