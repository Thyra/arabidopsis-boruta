# This script prepares the input data from 1745AJ for important trait selection by:
# - filtering columns to only contain top view traits
# - filter rows to start at DAS 21
# - Geting rid of empty columns
# - Filtering traits by the trait selection list
# - VIF filtering the traits
# - Reshaping the data to wide format (Trait X DAS)
# - detect outliers (2 SDs below or above median) and replace with the median of that specific treatment and day
# - remove any NA columns again
# The result is then saved to data/intermediary/1745AJ_Phenotyping_nafixed.csv

require('fmsb')
require('randomForest')

data = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")
# Filter columns to only contain top view traits
data = data[,c(1,6,10, grep("top.*", names(data)))]
# Filter rows to only contain rows starting at DAS 21
# this reduces our plant count to 113, sadly.
# Also limit data to DAS 35
data = data[data$DAS >= 21 & data$DAS <= 35,]
# Get rid of columns that only contain NAs
data = data[colSums(!is.na(data)) > 0]
# Filter out any traits that don't have a YES in the trait selection list (either because they have a NO or because they're not present)
trait_info = read.csv("data/input/trait_selection_list.csv")
keep_columns = intersect(
  names(data),
  c("Plant.ID", "Treatment", "DAS", make.names(trait_info[trait_info$Selection == "YES",]$Variable.ID)))
data = data[keep_columns]

# Reshape data to wide format
data = reshape(data, idvar='Plant.ID', timevar='DAS', v.names=setdiff(names(data), c("Plant.ID", "DAS", "Treatment")), direction='wide')
# Drop the plant ID
data = data[,setdiff(names(data), c("Plant.ID"))]
# Get rid of columns that only contain NAs (again)
data = data[colSums(!is.na(data)) > 0]

#### ------  OUTLIER DETECTION  --------- ###
data$Treatment = as.factor(data$Treatment)
data.drought = data[data$Treatment == "D",]
data.well_watered = data[data$Treatment == "WW",]

# Outlier detection -> replace with NA's which are then replaced with medians in following step.
medians = apply(data.drought[,setdiff(names(data.drought), c("Treatment"))], 2, FUN=median, na.rm=T)
sds = apply(data.drought[,setdiff(names(data.drought), c("Treatment"))], 2, FUN=sd, na.rm=T)
for(traitname in names(data.drought[,setdiff(names(data.drought), c("Treatment"))])) {
  print(traitname)
  if(nrow(data.drought[which(data.drought[traitname] > medians[traitname]+2*sds[traitname] | data.drought[traitname] < medians[traitname]-2*sds[traitname]), ]) > 0)
    data.drought[which(data.drought[traitname] > medians[traitname]+2*sds[traitname] | data.drought[traitname] < medians[traitname]-2*sds[traitname]), ][traitname] <- NA
}

# Remove NA's, currently just by replacing with column median but might want to use a more sophisticated
# method like rfImpute that weighs by proximity.
data.drought = na.roughfix(data.drought)

# Outlier detection -> replace with NA's which are then replaced with medians in following step.
medians = apply(data.well_watered[,setdiff(names(data.well_watered), c("Treatment"))], 2, FUN=median, na.rm=T)
sds = apply(data.well_watered[,setdiff(names(data.well_watered), c("Treatment"))], 2, FUN=sd, na.rm=T)
for(traitname in names(data.well_watered[,setdiff(names(data.well_watered), c("Treatment"))])) {
  print(traitname)
  if(nrow(data.well_watered[which(data.well_watered[traitname] > medians[traitname]+2*sds[traitname] | data.well_watered[traitname] < medians[traitname]-2*sds[traitname]), ]) > 0)
    data.well_watered[which(data.well_watered[traitname] > medians[traitname]+2*sds[traitname] | data.well_watered[traitname] < medians[traitname]-2*sds[traitname]), ][traitname] <- NA
}

# Remove NA's, currently just by replacing with column median but might want to use a more sophisticated
# method like rfImpute that weighs by proximity.
data.well_watered = na.roughfix(data.well_watered)

# Combine drought and well-watered data again
data.nafixed = rbind(data.drought, data.well_watered)
# Get rid of columns that contain NAs (because they were dropped for one treatment)
data.nafixed = data.nafixed[colSums(!is.na(data.nafixed)) == nrow(data.nafixed)]

write.csv(data.nafixed, "data/intermediary/1745AJ_Phenotyping_nafixed.csv", row.names=F)
