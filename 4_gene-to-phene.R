library(dtw)
library(Boruta)

trait.curves = read.csv("data/results/gaussian-process/trait_curves.lfs.csv.gz")
gene.curves  = read.csv("data/results/gaussian-process/gene_curves.lfs.csv.gz")

# If there are any trait vals earlier than the earliest RNA-Seq val, cut them off (should not be the case)
trait.curves = trait.curves[trait.curves$day >= min(gene.curves$day),]

# Weights is a vector with the same length as bait
# ref needs to be longer than bait
# WE SHOULD WRITE SOME UNIT TESTS FOR THIS! https://github.com/r-lib/vdiffr can even write tests for plot functions
# Attention: we don't need a max offset limit as long as gene and trait vals stop at the same day.
# Otherwise we do, in order to make sure gene on day 32 is not aligned with trait on day 30. (for example)
# @TODO get the weight thingy right (they are different as they are moving, no?) --> I had weirdly negative distance values??
align_curves <- function(ref, bait, weights.ref, weights.bait) {
  optimal.offset = NA
  min.dist = Inf
  for(o in 0:(nrow(ref)-nrow(bait))) {
    # Essentially Euclidian Distance weighted by weight vector
    current.dist = sum(abs(ref[(1+o):(o+nrow(bait)),]-bait)*cbind(weights.ref[(1+o):(o+nrow(bait))], weights.bait))
    #current.dist = sum(abs(ref[(1+o):(o+nrow(bait)),]-bait))
    if(current.dist < min.dist) {
      min.dist = current.dist
      optimal.offset = o
    }
  }
  return(list(offset = optimal.offset, dist = min.dist, ref=ref, bait=bait))
}

# Calculate coefficients of variation (mean of both treatments) to use for weighting the time points
gene.curves$avg.cv = (gene.curves$ww.std/gene.curves$ww.mean + gene.curves$d.std/gene.curves$d.mean)/2
trait.curves$avg.cv = (trait.curves$ww.std/trait.curves$ww.mean + trait.curves$d.std/trait.curves$d.mean)/2
trait.curves$avg.diff.cv = (trait.curves$ww.diff.std/abs(trait.curves$ww.diff.mean) + trait.curves$d.diff.std/abs(trait.curves$d.diff.mean))/2

curve.distances = data.frame(trait=character(), gene=character(), trait.mode=character(), gene.mode=character(), dist=numeric(), offset=numeric())
# @TODO make this into apply or something
for(trait in unique(trait.curves$trait)) {
  print(trait)
  trait.vals = scale(trait.curves[trait.curves$trait == trait,c("ww.mean", "d.mean")])
  trait.cv = trait.curves[trait.curves$trait == trait,]$avg.cv
  trait.diff.vals = na.omit(scale(trait.curves[trait.curves$trait == trait,c("ww.diff.mean", "d.diff.mean")]))
  trait.diff.cv = na.omit(trait.curves[trait.curves$trait == trait,]$avg.diff.cv)
  for(gene in unique(gene.curves$gene)) {
    gene.vals = scale(gene.curves[gene.curves$gene == gene,c("ww.mean", "d.mean")])
    gene.cv = gene.curves[gene.curves$gene == gene,]$avg.cv
    a = align_curves(gene.vals, trait.vals, 1/gene.cv, 1/trait.cv)
    curve.distances = rbind(curve.distances, list(trait=trait,gene=gene,trait.mode="normal",gene.mode="normal",
                                                  dist=a$dist, offset=a$offset))
    gene.vals = scale(-gene.curves[gene.curves$gene == gene,c("ww.mean", "d.mean")])
    gene.cv = gene.curves[gene.curves$gene == gene,]$avg.cv
    a = align_curves(gene.vals, trait.vals, 1/gene.cv, 1/trait.cv)
    curve.distances = rbind(curve.distances, list(trait=trait,gene=gene,trait.mode="normal",gene.mode="inverted",
                                                  dist=a$dist, offset=a$offset))
    
    gene.vals = scale(gene.curves[gene.curves$gene == gene,c("ww.mean", "d.mean")])
    gene.cv = gene.curves[gene.curves$gene == gene,]$avg.cv
    a = align_curves(gene.vals, trait.diff.vals, 1/gene.cv, 1/trait.diff.cv)
    curve.distances = rbind(curve.distances, list(trait=trait,gene=gene,trait.mode="diff",gene.mode="normal",
                                                  dist=a$dist, offset=a$offset))
    gene.vals = scale(-gene.curves[gene.curves$gene == gene,c("ww.mean", "d.mean")])
    gene.cv = gene.curves[gene.curves$gene == gene,]$avg.cv
    a = align_curves(gene.vals, trait.diff.vals, 1/gene.cv, 1/trait.diff.cv)
    curve.distances = rbind(curve.distances, list(trait=trait,gene=gene,trait.mode="diff",gene.mode="inverted",
                                                  dist=a$dist, offset=a$offset))
  }
}


plot_alignment <- function(a) {
  par(mfrow=2:1, mar=c(0,0,0,0))
  plot(x=seq_along(a$ref[,1]), y = a$ref[,1], type="l")
  lines(x=seq_along(a$bait[,1])+a$offset, y = a$bait[,1], lty=2)
  plot(x=seq_along(a$ref[,2]), y = a$ref[,2], type="l")
  lines(x=seq_along(a$bait[,2])+a$offset, y = a$bait[,2], lty=2)
}


write.csv(curve.distances, "data/results/gaussian-process/curve_distances.csv", row.names=F)
quit()

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
