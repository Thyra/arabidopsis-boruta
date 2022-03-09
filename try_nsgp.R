# You can use this seed to reproduce our exact results. If you change it, you might get a different set of genes
# or gene-trait connections but they should be just as valid candidates for drought stress response/recovery (see discusison)
set.seed(582477)
library(nsgp)
library(stringr)

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

write.csv(mothertable.transposed, "data/intermediary/mothertable_transposed.csv")
mothertable.transposed$day = as.numeric(str_split_fixed(str_split_fixed(mothertable.transposed$Mothertable.ID, "_", 3)[,1], "d" ,2)[,2])
mothertable.transposed$treatment = str_split_fixed(mothertable.transposed$Mothertable.ID, "_", 3)[,2]

# FOR THE GENES, we shall use these scores:
gene = sample(mothertable$target_id, 1)
ctrl = mothertable.transposed[mothertable.transposed$treatment == "WW", c("day", gene)]
names(ctrl) = c("day", "tpm")
case = mothertable.transposed[mothertable.transposed$treatment == "D", c("day", gene)]
names(case) = c("day", "tpm")
case = rbind(case, ctrl[ctrl$day == 15,])
res = gpr2sample(ctrl$day, ctrl$tpm, case$day, case$tpm, seq(14,37,1/8))
score.gene = (res$casemodel$targets$pmean - res$ctrlmodel$targets$pmean) * exp(res$ratios$npc)
score.gene = data.frame(x = res$ctrlmodel$targets$x, score=score.gene)
plot(res, plotratios = "npc")

pheno = read.csv("data/input/1745AJ_Phenotyping_formatted.lfs.csv.gz", sep=";")

pheno = pheno[pheno$DAS.FLOAT >= 20.43 & pheno$DAS.FLOAT < 31,]

trait = "top.intensity.vis.lab.a.mean"
# top.intensity.vis.lab.a.mean
# 
ctrl = pheno[pheno$Treatment == "WW", c("DAS", trait)]
names(ctrl) = c("day", "val")
case = pheno[pheno$Treatment == "D", c("DAS", trait)]
names(case) = c("day", "val")

res = gpr2sample(ctrl$day, ctrl$val, case$day, case$val, seq(20,30,1/8), nskernel = F)
score = (res$casemodel$targets$pmean - res$ctrlmodel$targets$pmean) * exp(res$ratios$npc)
score = data.frame(x = res$ctrlmodel$targets$x, score=score)
plot(res, plotratios = "npc")
