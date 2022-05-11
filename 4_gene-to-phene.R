library(dtw)
library(Boruta)

trait.curves = read.csv("data/results/gaussian-process/trait_curves.lfs.csv.gz")
gene.curves  = read.csv("data/results/gaussian-process/gene_curves.lfs.csv.gz")

# If there are any trait vals earlier than the earliest RNA-Seq val, cut them off (should not be the case)
trait.curves = trait.curves[trait.curves$day >= min(gene.curves$day),]

# ----- DTW ------
# Now we want to know how many RNA-Seq time points we have before the traits start (min 0),
# this is the offset to our window.type function below
window.offset = nrow(gene.curves[gene.curves$day < min(trait.curves$day) & gene.curves$gene == sample(gene.curves$gene, 1),])

dtw.distances = data.frame(trait=character(), gene=character(), distance=numeric())
pdf("plots/dtw-curves.lfs.pdf")
for(trait in unique(trait.curves$trait)) {
  for(gene in unique(gene.curves$gene)) {
    # Z normalize both curves
    t = scale(trait.curves[trait.curves$trait == trait,c("ww.mean", "d.mean")])
    g = scale(gene.curves[gene.curves$gene == gene,c("ww.mean", "d.mean")])
    # You can use a "rigid" step pattern to not allow any warping (i.e. the delay from gene -> phene is constant)
    # if that works well, it might make sense to write my own simple alignment function which weighs the distance by
    # confidence interval produced by the gaussian process
    a = dtw(t, g, open.end = T, open.begin = T, step.pattern = rigid,
            # This window.type function ensures that trait values cannot be mapped to later RNA-Seq values,
            # as that makes no sense biologically (trait change cannot come from RNA-Seq change later on)
            window.type = function(iw, jw, query.size, reference.size) { iw >= (jw - window.offset) } )
    plot(a, t, g, type="twoway", main=paste(trait, gene))
    dtw.distances = rbind(dtw.distances, data.frame(trait=trait, gene=gene, distance=a$distance))
  }
}
dev.off()
write.csv(dtw.distances, "data/results/gene-trait_dtw_distances.csv", row.names = F)
# 
# # ----- Boruta -----
# set.seed(20220310)
# boruta.weights = dtw.distances = data.frame(trait=character(), gene=character(), medianImp=numeric())
# for(trait in names(trait.scores)[2:ncol(trait.scores)]) {
#   boruta.df = merge(trait.scores[c("day", trait)], gene.scores)
#   boruta.df = boruta.df[2:ncol(boruta.df)]
#   names(boruta.df)[1] = "trait"
#   boruta_result = Boruta(trait ~ ., data=boruta.df, maxRuns = 1500, doTrace = 2)
#   trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
#   trait_importance_roughfixed = cbind(gene = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
#   trait_importance_roughfixed$trait = trait
#   trait_importance_roughfixed = trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed", c("trait", "gene", "medianImp")]
#   boruta.weights = rbind(boruta.weights, trait_importance_roughfixed)
# }
# write.csv(boruta.weights, "data/results/gene-trait_boruta_weights.csv", row.names = F)
