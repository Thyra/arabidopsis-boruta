library(UpSetR)
library(grid)

pdf("plots/upset_plots_selected_genes.pdf", 20, 15)

# Plot gene overlaps from treatment by genexp classification
treatment_from_genes.selected_genes = list()
for(run in list.files("data/results/cv/selected_genes/", pattern = "treatment_from_genes_*")) {
  treatment_from_genes.selected_genes[[sub("treatment_from_genes_", "", run)]] = read.csv(paste0("data/results/cv/selected_genes/", run))$gene
}
upset(fromList(treatment_from_genes.selected_genes), order.by = "freq", nsets=27, nintersects=NA)
grid.text("Treatment from Gene Expression",x = 0.65, y=0.95, gp=gpar(fontsize=20))

dev.off()
