library(nsgp)
library(stringr)

mothertable.transposed = read.csv("data/intermediary/mothertable_transposed.csv")
pheno = read.csv("data/intermediary/1745AJ_Phenotyping_cleaned.csv")
# Only look at stress phase (not recovery), otherwise the gaussian process interpolates the gene expression change
# already before the rewatering starts (as there are no RNA measurements at day 30)
pheno = pheno[pheno$DAS < 31,]
mothertable.transposed = mothertable.transposed[mothertable.transposed$day < 31 | mothertable.transposed$treatment == "WW",]

# FOR TEST PURPOSES: ONLY DO FIRST 300 genes!
mothertable.transposed = mothertable.transposed[1:303]

pdf("plots/gaussian-process_genes.lfs.pdf")
gene.curves = data.frame(gene=character(),day=numeric(),ww.mean=numeric(),ww.std=numeric(),d.mean=numeric(),d.std=numeric())
for(gene in names(mothertable.transposed)[4:ncol(mothertable.transposed)]) {
  print(gene)
  tryCatch(
    {
      ctrl = mothertable.transposed[mothertable.transposed$treatment == "WW", c("day", gene)]
      names(ctrl) = c("day", "tpm")
      case = mothertable.transposed[mothertable.transposed$treatment == "D", c("day", gene)]
      names(case) = c("day", "tpm")
      res = gpr2sample(ctrl$day, ctrl$tpm, case$day, case$tpm, seq(14,max(mothertable.transposed$day),1/8))
      plot(res, plotratios = "npc", title=gene)
      curves = as.data.frame(cbind(res$ctrlmodel$targets$pmean, res$ctrlmodel$targets$pstd,
                                   res$casemodel$targets$pmean, res$casemodel$targets$pstd))
      names(curves) = c("ww.mean", "ww.std", "d.mean", "d.std")
      curves$day = seq(14,max(mothertable.transposed$day),1/8)
      curves$gene = gene
      gene.curves = rbind(gene.curves, curves[curves$day <= 30,])
    }, error=function(cond) {
      print(" ERROR")
    }
  )
}
dev.off()

write.csv(gene.curves, gzfile("data/results/gaussian-process/gene_curves.lfs.csv.gz"), row.names=F)

trait_importance = read.csv("data/results/trait_selection_result_roughFixed.csv")
important.traits = unique(sub("\\.\\d+$", "", trait_importance[trait_importance$decision == "Confirmed",]$trait))
pheno.curves = data.frame(trait=character(),day=numeric(),ww.mean=numeric(),ww.std=numeric(),
                          d.mean=numeric(),d.std=numeric(),ww.diff.mean=numeric(),ww.diff.std=numeric(),
                          d.diff.mean=numeric(),d.diff.std=numeric())
pdf("plots/gaussian-process_traits.pdf")
for(trait in important.traits) {
  # For each trait we will do two regressions: One on the actual values and one on the change in comparison
  # to the previous day. That is because some genes will probably have a cumulative effect, so you would expect
  # their expression to correlate with the change in phenotype rather than the phenotype itself (e.g. when a growth
  # hormone stops being expressed, the plant does not shrink.)
  print(trait)
  ctrl = pheno[pheno$Treatment == "WW", c("Plant.ID","DAS", trait)]
  ctrl = ctrl[order(ctrl$Plant.ID, ctrl$DAS),]
  names(ctrl) = c("Plant.ID", "day", "val")
  ctrl$diff <- ave(ctrl$val, ctrl$Plant.ID, FUN=function(x) c(NA, diff(x)))
  
  case = pheno[pheno$Treatment == "D", c("Plant.ID", "DAS", trait)]
  case = case[order(case$Plant.ID, case$DAS),]
  names(case) = c("Plant.ID", "day", "val")
  case$diff <- ave(case$val, case$Plant.ID, FUN=function(x) c(NA, diff(x)))

  res.val = gpr2sample(ctrl$day, ctrl$val, case$day, case$val, seq(20,30,1/8), nskernel = F)
  plot(res.val, plotratios = "npc", title=trait)
  res.diff = gpr2sample(ctrl$day, ctrl$diff, case$day, case$diff, seq(20,30,1/8), nskernel = F)
  plot(res.diff, plotratios = "npc", title=paste0(trait, " - diff"))
  
  curves = as.data.frame(cbind(res.val$ctrlmodel$targets$pmean, res.val$ctrlmodel$targets$pstd,
                               res.val$casemodel$targets$pmean, res.val$casemodel$targets$pstd,
                               res.diff$ctrlmodel$targets$pmean, res.diff$ctrlmodel$targets$pstd,
                               res.diff$casemodel$targets$pmean, res.diff$casemodel$targets$pstd))
  names(curves) = c("ww.mean", "ww.std", "d.mean", "d.std", "ww.diff.mean", "ww.diff.std", "d.diff.mean", "d.diff.std")
  curves$day = seq(20,30,1/8)
  curves$trait = trait
  pheno.curves = rbind(pheno.curves, curves)
}
dev.off()
# Remove diff values before day 21, because we don't really have any for day 20.
pheno.curves[pheno.curves$day < 21,c("ww.diff.mean","ww.diff.std","d.diff.mean","d.diff.std")] = NA

write.csv(trait.scores, "data/results/gaussian-process/trait_curves.lfs.csv", row.names=F)
