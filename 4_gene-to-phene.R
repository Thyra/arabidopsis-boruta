library(dtw)
library(Boruta)

gene.scores = read.csv("data/results/gaussian-scores/gene_scores.lfs.csv")
gene.scores = gene.scores[colSums(is.na(gene.scores)) == 0]
trait.scores = read.csv("data/results/gaussian-scores/trait_scores.lfs.csv")

# If there are any trait vals earlier than the earliest RNA-Seq val, cut them off (should not be the case)
trait.scores = trait.scores[trait.scores$day >= min(gene.scores$day),]

# ----- DTW ------
# Now we want to know how many RNA-Seq time points we have before the traits start (min 0),
# this is the offset to our window.type function below
window.offset = nrow(gene.scores[gene.scores$day < min(trait.scores$day),])

dtw.distances = data.frame(trait=character(), gene=character(), distance=numeric())
pdf("plots/dtw-curves.lfs.pdf")
for(trait in names(trait.scores)[2:ncol(trait.scores)]) {
  for(gene in names(gene.scores)[2:ncol(gene.scores)]) {
    # Z normalize both curves
    t = scale(trait.scores[[trait]])
    g = scale(gene.scores[[gene]])
    a = dtw(t, g, open.end = T, open.begin = T, step.pattern = asymmetric,
            # This window.type function ensures that trait values cannot be mapped to later RNA-Seq values,
            # as that makes no sense biologically (trait change cannot come from RNA-Seq change later on)
            window.type = function(iw, jw, query.size, reference.size) { iw >= (jw - window.offset) } )
    plot(a, t, g, type="twoway", main=paste(trait, gene))
    dtw.distances = rbind(dtw.distances, data.frame(trait=trait, gene=gene, distance=a$distance))
  }
}
dev.off()
write.csv(dtw.distances, "data/results/gene-trait_dtw_distances.csv", row.names = F)

# ----- Boruta -----
set.seed(20220310)
boruta.weights = dtw.distances = data.frame(trait=character(), gene=character(), medianImp=numeric())
for(trait in names(trait.scores)[2:ncol(trait.scores)]) {
  boruta.df = merge(trait.scores[c("day", trait)], gene.scores)
  boruta.df = boruta.df[2:ncol(boruta.df)]
  names(boruta.df)[1] = "trait"
  boruta_result = Boruta(trait ~ ., data=boruta.df, maxRuns = 1500, doTrace = 2)
  trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
  trait_importance_roughfixed = cbind(gene = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
  trait_importance_roughfixed$trait = trait
  trait_importance_roughfixed = trait_importance_roughfixed[trait_importance_roughfixed$decision == "Confirmed", c("trait", "gene", "medianImp")]
  boruta.weights = rbind(boruta.weights, trait_importance_roughfixed)
}
write.csv(boruta.weights, "data/results/gene-trait_boruta_weights.csv", row.names = F)
