library(dtw)

gene.scores = read.csv("data/results/gaussian-scores/gene_scores.lfs.csv")
trait.scores = read.csv("data/results/gaussian-scores/trait_scores.lfs.csv")

for(trait in names(trait.scores)[2:ncol(trait.scores)]) {
  for(gene in names(gene.scores)[2:ncol(gene.scores)]) {
    # Z normalize both curves
    t = scale(trait.scores[[trait]])
    g = scale(gene.scores[[gene]])
    a = dtw(t, g, open.end = T, open.begin = T, step.pattern = asymmetric)
    plot(a, t, g, type="twoway", main=gene)
  }
}
