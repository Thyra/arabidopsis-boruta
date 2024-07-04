# This script prepares the input data from 1745AJ for important trait selection by:
# - filtering columns to only contain top view traits
# - filter rows to start at DAS 20 after the camera configuration switch
# - Geting rid of empty columns
# - Filtering traits by the trait selection list
# - Reshaping the data to wide format (Trait X DAS)
# - detect outliers (2 SDs below or above median) and replace with the median of that specific treatment and day
# - remove any NA columns again
# The result is then saved to data/intermediary/1745AJ_Phenotyping_nafixed.csv

require('fmsb')
require('randomForest')

data = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")
# Filter columns to only contain top view traits
data = data[,c(1,6,10,11, grep("top.*", names(data)))]
# Filter rows to only contain rows starting at DAS 20 after the camera configuration change
# this reduces our plant count to 113, sadly.
# Also limit data to DAS 35 because they become very noisy afterwards because of the baskets.
data = data[data$DAS.FLOAT >= 20.45 & data$DAS <= 35,]
# Get rid of columns that only contain NAs
data = data[colSums(!is.na(data)) > 0]
# Filter out any traits that don't have a YES in the trait selection list (either because they have a NO or because they're not present)
trait_info = read.csv("data/input/trait_selection_list.csv")
keep_columns = intersect(
  names(data),
  c("Plant.ID", "Treatment", "DAS", make.names(trait_info[trait_info$Selection == "YES",]$Variable.ID)))
data = data[keep_columns]

# Reshape data to wide format
data = reshape(data, idvar='Plant.ID', timevar='DAS', v.names=setdiff(names(data), c("Plant.ID", "DAS", "Treatment")), direction='wide', sep="_")

# Get rid of columns that only contain NAs (again)
data = data[colSums(!is.na(data)) > 0]

#### ------  OUTLIER DETECTION  --------- ###
data$Treatment = as.factor(data$Treatment)

replace_outliers <- function(df, clnames) {
  medians = apply(df[,clnames], 2, FUN=median, na.rm=T)
  sds = apply(df[,clnames], 2, FUN=sd, na.rm=T)
  for(traitname in clnames) {
    print(traitname)
    if(nrow(df[which(df[traitname] > medians[traitname]+2*sds[traitname] | df[traitname] < medians[traitname]-2*sds[traitname]), ]) > 0)
      df[which(df[traitname] > medians[traitname]+2*sds[traitname] | df[traitname] < medians[traitname]-2*sds[traitname]), ][traitname] <- NA
  }
  df[clnames] = na.roughfix(df[clnames])
  return(df)
}

data.drought = data[data$Treatment == "D",]
data.well_watered = data[data$Treatment == "WW",]

columns_to_fix = setdiff(names(data.drought), c("Plant.ID", "Treatment"))

# Combine drought and well-watered data again
data.nafixed = rbind(replace_outliers(data.drought, columns_to_fix), 
                     replace_outliers(data.well_watered, columns_to_fix))
# Get rid of columns that contain NAs (because they were dropped for one treatment)
data.nafixed = data.nafixed[colSums(!is.na(data.nafixed)) == nrow(data.nafixed)]

# Reshape back to long format
pheno.fixed = reshape(data.nafixed, direction="long", varying = columns_to_fix, timevar="DAS", times=20:35, sep="_", idvar="Plant.ID")

write.csv(pheno.fixed, "data/intermediary/1745AJ_Phenotyping_nafixed.csv", row.names=F)
