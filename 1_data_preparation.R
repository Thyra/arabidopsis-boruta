require('statgenHTP')
require('stringr')
# ----- Phenotyping Data ----------------------
data = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")
# Filter columns to only contain top view traits
data = data[,c(1,6,9,10,11, grep("top.*", names(data)))]
# Filter rows to only contain those after the camera configuration shift and before baskets were installed
data = data[data$DAS.FLOAT >= 20.43 & data$DAS <= 31,]
# Get rid of columns that only contain NAs
data = data[colSums(!is.na(data)) > 0]
# Filter out any traits that don't have a YES in the trait selection list (either because they have a NO or because they're not present)
trait_info = read.csv("data/input/trait_selection_list.csv")
keep_columns = intersect(
  names(data),
  c("Plant.ID", "Treatment", "Time", "DAS", make.names(trait_info[trait_info$Selection == "YES",]$Variable.ID)))
data = data[keep_columns]

## Outlier Detection with statgenHTP

# First we'll add some positional columns (row and column)
data$plant.id = as.numeric(substr(data$Plant.ID, 7,9))
data$print.pot_row = floor((data$plant.id -1) / 32)+1
data$print.pot_column = floor((data$plant.id -1) %% 32)+1
# For every second row the trays should start from the right because the lanes are connected in one big curve
data[data$print.pot_row %% 2 == 0,]$print.pot_column = 34 - data[data$print.pot_row %% 2 == 0,]$print.pot_column
# Make a lane between row 6 and 7
data[data$print.pot_row > 6,]$print.pot_row = data[data$print.pot_row > 6,]$print.pot_row + 1

phenoTP <- createTimePoints(dat = data,
                            experimentName = "1745AJ",
                            genotype = "Treatment",
                            timePoint = "Time",
                            timeFormat = "%d.%m.%y",
                            repId = "plant.id",
                            plotId = "Plant.ID",
                            rowNum = "print.pot_row", colNum = "print.pot_column"
                            )
for(trait in c("top.geometry.fluo.area..px.2.", "top.geometry.fluo.leaf.count",
               "top.geometry.fluo.hull.length", "top.intensity.nir.mean", "top.intensity.vis.hsv.h.mean", "top.intensity.vis.hsv.s.mean",
               "top.intensity.vis.hsv.v.mean", "top.intensity.vis.hsv.h.yellow2green", "top.intensity.fluo.intensity.mean")) {
  plot(phenoTP, 
       traits = trait,
       plotType = "raw", plotLine=T)
  
  # Single outlier points
  # plot(phenoTP, 
  #      traits = "top.intensity.fluo.hsv.v.mean",
  #      plotType = "raw", plotLine=T)
  # 
  # singleOut <- detectSingleOut(TP = phenoTP,
  #                              trait = "top.intensity.fluo.hsv.v.mean",
  #                              confIntSize = 3,
  #                              nnLocfit = 0.3)
  # plot(singleOut, outOnly = FALSE)
  
  
  # Whole Series
  fit.spline <- fitSpline(inDat = as.data.frame(phenoTP),
                          trait = trait,
                          genotypes = "WW",
                          knots = 50,
                          useTimeNumber = TRUE,
                          timeNumber = "DAS")
  predDat <- fit.spline$predDat
  coefDat <- fit.spline$coefDat
  
  outVator <- detectSerieOut(corrDat = as.data.frame(phenoTP),
                             predDat = predDat,
                             coefDat = coefDat,
                             trait = trait,
                             genotypes = "WW",
                             thrCor = 0.8,
                             thrPca = 30,
                             thrSlope = 0.7)
  plot(outVator)
  
  fit.spline <- fitSpline(inDat = as.data.frame(phenoTP),
                          trait = trait,
                          genotypes = "D",
                          knots = 50,
                          useTimeNumber = TRUE,
                          timeNumber = "DAS")
  predDat <- fit.spline$predDat
  coefDat <- fit.spline$coefDat
  
  outVator <- detectSerieOut(corrDat = as.data.frame(phenoTP),
                             predDat = predDat,
                             coefDat = coefDat,
                             trait = trait,
                             genotypes = "D",
                             thrCor = 0.8,
                             thrPca = 30,
                             thrSlope = 0.7)
  plot(outVator)
  
  readline(prompt="Next trait?")
}

# Based on the plots above, I'm taking out plant 188 and 085
data = data[!data$Plant.ID %in% c("1745AJ188", "1745AJ085"),]

out = data[,!(names(data) %in% c("Time", "plant.id", "print.pot_row", "print.pot_column"))]
write.csv(out, "data/intermediary/1745AJ_Phenotyping_cleaned.csv", row.names=F)

# ------- RNASeq-Data --------
mothertable = read.csv("data/input/curated_mothertable20200515.lfs.txt.gz", sep="\t")
# Remove any genes are not differentially expressed between drought and WW treatment
mothertable = mothertable[mothertable$differentialinatleastoneconstrast == "TRUE",]
# We're only using the TPM values (Stephan used counts but Andrea suggested to go with TPM)
mothertable.tpm = mothertable[43:81]
colnames(mothertable.tpm) <- sub("_tpm$", "", colnames(mothertable.tpm)) # remove _tpm from the end of the columns
row.names(mothertable.tpm) <- mothertable$target_id
# Transpose it so that we have "HarvestX_repY" --> "Gene 1", "Gene 2", ...
mothertable.transposed = as.data.frame(t(mothertable.tpm))
mothertable.transposed$Mothertable.ID = row.names(mothertable.transposed)
# Add day and treatment info
mothertable.transposed$day = as.numeric(str_split_fixed(str_split_fixed(mothertable.transposed$Mothertable.ID, "_", 3)[,1], "d" ,2)[,2])
mothertable.transposed$treatment = str_split_fixed(mothertable.transposed$Mothertable.ID, "_", 3)[,2]

# Copy the WW samples from 15 DAS as D samples as well to give a good start point (there is no drought stress there for sure)
add = mothertable.transposed[mothertable.transposed$day == 15,]
add$treatment = "D"
mothertable.transposed = rbind(mothertable.transposed, add)

mothertable.transposed = mothertable.transposed[c((ncol(mothertable.transposed)-2):(ncol(mothertable.transposed)), 1:(ncol(mothertable.transposed)-3))]
write.csv(mothertable.transposed, "data/intermediary/mothertable_transposed.csv", row.names=F)


