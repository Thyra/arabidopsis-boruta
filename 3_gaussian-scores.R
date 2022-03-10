library(nsgp)
library(stringr)

mothertable.transposed = read.csv("data/intermediary/mothertable_transposed.csv")
pheno = read.csv("data/intermediary/1745AJ_Phenotyping_cleaned.csv")
# Only look at stress phase (not recovery)
pheno = pheno[pheno$DAS < 31,]
mothertable.transposed = mothertable.transposed[mothertable.transposed$day < 31 | mothertable.transposed$treatment == "WW",]

# FOR TEST PURPOSES: ONLY DO FIRST 20 genes!
mothertable.transposed = mothertable.transposed[1:23]

gene.scores = as.data.frame(sapply(names(mothertable.transposed)[4:ncol(mothertable.transposed)], function(gene) {
  print(gene)
  tryCatch(
    {
      ctrl = mothertable.transposed[mothertable.transposed$treatment == "WW", c("day", gene)]
      names(ctrl) = c("day", "tpm")
      case = mothertable.transposed[mothertable.transposed$treatment == "D", c("day", gene)]
      names(case) = c("day", "tpm")
      res = gpr2sample(ctrl$day, ctrl$tpm, case$day, case$tpm, seq(14,max(mothertable.transposed$day),1/8))
      plot(res, plotratios = "npc")
      (res$casemodel$targets$pmean - res$ctrlmodel$targets$pmean) * exp(res$ratios$npc)
    }, error=function(cond) {
      seq(14,max(mothertable.transposed$day),1/8) * NA
    }
  )
}))
gene.scores$day = seq(14,max(mothertable.transposed$day),1/8)
gene.scores = gene.scores[gene.scores$day <= 30,]
gene.scores = gene.scores[c(ncol(gene.scores), (1:ncol(gene.scores)-1))]

write.csv(gene.scores, "data/results/gaussian-scores/gene_scores.lfs.csv", row.names=F)

trait_importance = read.csv("data/results/trait_selection_result_roughFixed.csv")
important.traits = unique(sub("\\.\\d+$", "", trait_importance[trait_importance$decision == "Confirmed",]$trait))
trait.scores = as.data.frame(sapply(important.traits, function(trait) {
  print(trait)
  ctrl = pheno[pheno$Treatment == "WW", c("DAS", trait)]
  names(ctrl) = c("day", "val")
  case = pheno[pheno$Treatment == "D", c("DAS", trait)]
  names(case) = c("day", "val")
  res = gpr2sample(ctrl$day, ctrl$val, case$day, case$val, seq(20,30,1/8), nskernel = F)
  plot(res, plotratios = "npc")
  (res$casemodel$targets$pmean - res$ctrlmodel$targets$pmean) * exp(res$ratios$npc)
}))
trait.scores$day = seq(20,30,1/8)
trait.scores = trait.scores[c(ncol(trait.scores), (1:ncol(trait.scores)-1))]

write.csv(trait.scores, "data/results/gaussian-scores/trait_scores.lfs.csv", row.names=F)
