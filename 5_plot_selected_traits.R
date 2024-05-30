# This script plots the results of the boruta Trait selection:
# - plot all traits that are important at any day over the whole time course
# - plot cumulative importance of all confirmed traitXDAS combinations by trait category over time

traits = read.csv("data/results/trait_selection_result_roughFixed.csv")
traits$trait.name = sub("\\.\\d+$", "", traits$trait)
traits$DAS = gsub("^.*\\.", "", traits$trait)

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
# Plot all the important traits

library(ggplot2)
library(Rmisc)
library(gridExtra)
pdf("plots/selected_traits.pdf", height=4)
for(t in keep_traits) {
  cis = group.CI(get(t) ~ DAS + Treatment, pheno, ci=0.84)
  names(cis) = c("DAS", "Treatment", "upper", "mean", "lower")
  merged2 = merge(cis, traits[traits$trait.name == t,], by="DAS", all.x=T)
  merged2$impY = (merged2$medianImp/max(traits$medianImp, na.rm=T)) * (max(merged2$upper, na.rm=T) - min(merged2$lower, na.rm=T)) + min(merged2$lower, na.rm=T)
  p1 <- ggplot(merged2, aes(x = DAS, y = mean, group = Treatment, color=Treatment)) +
    ggtitle(t) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group=Treatment), fill="grey70", alpha=0.5, colour=NA) +
    geom_line() +
    geom_point(aes(x = DAS, y = impY), color="black") +
    theme_minimal() +
#    theme(legend.position = "none") +
    scale_y_continuous(sec.axis = sec_axis(trans=function(x) {(x - min(merged2$lower, na.rm=T))/(max(merged2$upper, na.rm=T) - min(merged2$lower, na.rm=T)) * max(traits$meanImp, na.rm=T)}, name = "Trait Importance"))

  # Mean Differences
  meanDiffs = as.data.frame(do.call(rbind, lapply(
    split(pheno, pheno$DAS),
    function(x) {
      n.obs = aggregate(x[[t]], by=list(x$Treatment), FUN=function(y) length(which(!is.na(y))))$x
      if(n.obs[1] < 2 | n.obs[2] < 2) {
        c(NA, NA)
      } else {
        t.test(x[[t]] ~ x[["Treatment"]])$conf.int
      }
    }
  )))
  names(meanDiffs) = c("upper", "lower")
  meanDiffs$DAS = as.integer(rownames(meanDiffs))
  meanDiffs$mean = (meanDiffs$upper + meanDiffs$lower)/2
  p2 <- ggplot(meanDiffs, aes(x = DAS, y = mean)) +
          ggtitle("Difference in means (D - WW)") +
          geom_ribbon(aes(ymin = lower, ymax = upper), fill="grey70", alpha=0.5, colour=NA) +
          geom_line() +
          geom_hline(yintercept=0)
#  grid.arrange(p1, p2)
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

pdf("../results/Cumulative_Weights_per_day.pdf", width = 10)
print(ggplot(impCurves, aes(DAS, medianImp, fill = imaging.modality)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2"))
print(ggplot(impCurves, aes(DAS, medianImp, fill = imaging.modality)) +
        geom_bar(stat="identity", position = "stack") +
        scale_fill_manual(values=c("#009E73","#CC79A7","#0072B2","#D55E00"))) + labs(fill="Trait Category") +
        geom_segment(x = 7.5, xend=7.5, y=0, yend=60, linetype="dashed", color = "black", size=.7) +
  theme_minimal() +
  geom_segment(x = 8.5, xend=8.5, y=0, yend=26, linetype="dashed", color = "black", size=.7) +
  annotate("text", x=7.8, y=52, label= "Rewatering", hjust=0) +
  annotate("text", x=8.8, y=22, label= "Baskets", hjust=0)
dev.off()

