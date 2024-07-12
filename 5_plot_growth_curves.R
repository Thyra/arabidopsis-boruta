data = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")

field.capacities = read.csv("./data/input/FieldCapacity.csv", sep=";")
field.capacities$yCoordinate = field.capacities$Field.Capacity * 20 - 1000
field.capacities.drought = field.capacities[field.capacities$Treatment == "D",]
field.capacities.drought.medians = aggregate(yCoordinate ~ DAS, data = field.capacities.drought, median)

col.orange = "#e69f00"
col.darkblue = "#0072b2"

# [in DAYS.FLOAT]
drought.stress.start = 15
drought.stress.end   = 30

# Sort Dataframe by Plant.ID / DAS
data = data[
  with(data, order(Plant.ID, DAS.FLOAT)),
  ]

# Calculate Diff to previous Day
data <- within(data,difference <- ave(ave(top.geometry.vis.area.norm..mm.2., Plant.ID, DAS, FUN = mean),Plant.ID, FUN = function(x) c(NA,diff(x))))

# Remove unneeded columns
data = data[,c("Plant.ID", "DAS", "DAS.FLOAT", "difference", "Treatment")]

# Remove the second measurement of Day 20.
data = data[!(data$DAS == 20 & data$difference == 0),]

source("./colors.R")
data$color = stress_day_colors[1]
data[data$Treatment == "D", ]$color = stress_day_colors[data[data$Treatment == "D", ]$DAS + 1]

plants.without.flowers = c(
  "1745AJ140",
  "1745AJ175", # in this one also 1 leaf was removed
  "1745AJ226",
  "1745AJ254",
  "1745AJ269",
  "1745AJ273",
  "1745AJ277",
  "1745AJ287",
  "1745AJ299",
  "1745AJ316",
  "1745AJ384"
)
# Remove the datapoint at day 32 where the baskets where added except for the plants without flowers because they had no baskets.
data[data$DAS == 32 & !(data$Plant.ID %in% plants.without.flowers),]$difference = NaN

pdf("plots/growth_curves.pdf", 6.692, 4.486, colormodel="cmyk")
par(family="serif", mar=c(4,3.3,0,0))
plot(data$difference ~ data$DAS.FLOAT, col=as.factor(data$Treatment), axes=F, type="n", xlab="",ylab="", xlim=c(9,44), cex.lab=0.75)
title(xlab="Days after Sowing", cex.lab=0.75, line = 2.2)

# Drought stress period
rect(drought.stress.start, -1000, drought.stress.end, 800, col = "#F0F0F0", border = NA)
text(23.5, 700, "Drought Stress", col="#707070", cex=1.05)
lines(field.capacities.drought.medians$yCoordinate ~ field.capacities.drought.medians$DAS, lty=3, col="#707070")
text(12, field.capacities.drought.medians[field.capacities.drought.medians$DAS == 12,]$yCoordinate, "Median Field Capacity\nof drought stressed plants",col="#707070", pos = 1, xpd=T, cex=0.6)
text(29, min(field.capacities.drought.medians$yCoordinate), 
    paste(round((min(field.capacities.drought.medians$yCoordinate)+1000)/20), "%"),
    col="#707070", pos=4, cex=0.6)
text(15, field.capacities.drought.medians[field.capacities.drought.medians$DAS == drought.stress.start,]$yCoordinate, 
     paste(round((field.capacities.drought.medians[field.capacities.drought.medians$DAS == drought.stress.start,]$yCoordinate+1000)/20), "%"),
     col="#707070", pos=4, cex=0.6)

for(p in unique(data$Plant.ID)) {
  plant_data = data[data$Plant.ID == p,]
  #lines(plant_data$difference ~ plant_data$DAS, col=plant_data$color, lwd=0.6)
  segments(x0 = head(plant_data$DAS,-1),
         y0 = head(plant_data$difference,-1),
         x1 = tail(plant_data$DAS,-1),
         y1 = tail(plant_data$difference,-1),
         lwd = 1,
         col = plant_data$color)
}
axis(1, at=c(seq(10, 44, 5), 44), tick=F, cex.axis=0.75)
axis(1, at=seq(10,44,1), col = 0, col.ticks = 1, labels=F, cex.axis=0.75)
axis(2,tick=F,las=2, at=seq(-1000,1000,200), labels=parse(text=paste(seq(-1000, 1000, 200)/100, "~cm^2")), line=-0.5, cex.axis=0.75)
text(9,1000, "Plant Growth per Day", pos=4, xpd=T, cex=1.25)
text(9,900, "Difference to previous day in total area measured from top view", pos=4, xpd=T, cex=0.75)

text(28, -50, "Drought Stressed Plants", pos=2, col=stress_day_colors[30], cex=0.6)
text(28, 570, "Well Watered Plants", pos=2, col=stress_day_colors[1], cex=0.6)

text(32,140,"Installed\nBaskets",cex=0.55)

#text(29, data[data$DAS == 29 & data$Plant.ID == "1745AJ254",]$difference-10, "Plant 1745AJ254", pos=4, cex=0.35, col = col.orange)

abline(h = 0)

dev.off()
