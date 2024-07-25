# This script plots the results of the boruta Trait selection:
# - plot all traits that are important at any day over the whole time course
# - plot cumulative importance of all confirmed traitXDAS combinations by trait category over time

traits = read.csv("data/results/trait_selection_result_roughFixed.csv")
traits$trait.name = sub("_\\d+$", "", traits$trait)
traits$DAS = gsub("^.*\\_", "", traits$trait)

# Filter out traits which don't have Confirmed at any DAS
# --> i.e. traits that are not important at any point in time
keep_traits = unique(traits[traits$decision == "Confirmed",]$trait.name)
traits = traits[traits$trait.name %in% keep_traits,]

# Add Trait Category
trait_categories = read.csv("data/input/trait_selection_list.csv")
trait_categories$trait.name = make.names(trait_categories$Variable.ID)
merged = merge(traits, trait_categories)

# Load Phenotpying data
pheno = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")
# Limit to what I want to show
pheno = pheno[pheno$DAS.FLOAT >= 10 & pheno$DAS <= 35,]

# Plot all the important traits
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)

source("./colors.R")
color_stops = ((10:36) - 10) / 26
colors = stress_day_colors[10:36]

pdf("plots/selected_traits.pdf", height=4)
for(t in keep_traits) {
  cis = group.CI(get(t) ~ DAS + Treatment, pheno, ci=0.84)
  names(cis) = c("DAS", "Treatment", "upper", "mean", "lower")
  merged2 = merge(cis, traits[traits$trait.name == t,], by="DAS", all.x=T)
  merged2$impY = (merged2$medianImp/max(traits$medianImp, na.rm=T)) * (max(merged2$upper, na.rm=T) - min(merged2$lower, na.rm=T)) + min(merged2$lower, na.rm=T)
  merged2$color = merged2$DAS
  merged2[merged2$Treatment == "WW",]$color = 10
  p1 <- ggplot(merged2, aes(x = DAS, y = mean, group = Treatment, color=color)) +
    ggtitle(t) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group=Treatment), fill="grey70", alpha=0.5, colour=NA) +
    geom_line(linewidth = 1.5) +
    scale_colour_gradientn(values = color_stops , colors=colors) +
    geom_point(aes(x = DAS, y = impY), color="black") +
    theme_minimal() +
    theme(legend.position = "none") +
    xlim(10, 35) +
    scale_y_continuous(sec.axis = sec_axis(trans=function(x) {(x - min(merged2$lower, na.rm=T))/(max(merged2$upper, na.rm=T) - min(merged2$lower, na.rm=T)) * max(traits$meanImp, na.rm=T)}, name = "trait importance"))
  print(p1)
}
dev.off()

# Only use traits on the day they are important
merged = merged[merged$decision == "Confirmed",]
# Set intensity to color because they should be plotted together
merged[merged$trait.type == "intensity",]$trait.type = "color"
# then aggregate(medianImp ~ trait.category + DAS, FUN=SUM)
impCurves.architectural = aggregate(medianImp ~ trait.type + DAS, data=merged[merged$trait.type == "architectural",], FUN=sum)
impCurves.color = aggregate(medianImp ~ trait.type + DAS + imaging.modality, data=merged[merged$trait.type == "color",], FUN=sum)
impCurves.architectural$imaging.modality = "architectural (all)"
impCurves = rbind(impCurves.architectural, impCurves.color)

impCurves[impCurves$imaging.modality == "NIR", ]$imaging.modality = "near-infrared"
impCurves[impCurves$imaging.modality == "RGB", ]$imaging.modality = "visible light"
impCurves[impCurves$imaging.modality == "FLUOR", ]$imaging.modality = "fluorescence"

muh_grob <- grid::rectGrob(
  x=1:12, y=0, gp=gpar(
    col=NA, fill=stress_day_colors[24:35], alpha=1))

pdf("plots/cumulated_trait_importance.pdf", width = 10)
print(ggplot(impCurves, aes(DAS, medianImp, fill = imaging.modality)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values = unname(trait_category_colors))) + labs(fill = "Trait Category") +
  geom_segment(x = 7.5, xend=7.5, y=0, yend=60, linetype="dashed", color = "black", linewidth=.7) +
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  geom_segment(x = 8.5, xend=8.5, y=0, yend=26, linetype="dashed", color = "black", linewidth=.7) +
  annotate("text", x=7.8, y=52, label= "Rewatering", hjust=0, size=5) +
  annotate("text", x=8.8, y=22, label= "Baskets", hjust=0, size=5) +
  annotation_custom(
    grob=muh_grob, xmin = 0, xmax = 1, ymin = -2.5, ymax=0
  ) +
  ylab("cumulated median importance")
dev.off()
