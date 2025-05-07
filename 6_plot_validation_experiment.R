library(Rmisc)
library(ggplot2)
library(mratios)
library(RColorBrewer)
library(gridExtra)
library(grid)

plot.curve <- function(selected.genotypes, trait) {
  cis.list <- list()
  
  # Calculate CIs for each genotype
  for (genotype in selected.genotypes) {
    pheno.genotype = pheno[pheno$Genotype == genotype,]
    cis = group.CI(get(trait) ~ DAS + Treatment, pheno.genotype, ci=0.84)
    names(cis) = c("DAS", "Treatment", "upper", "mean", "lower")
    cis$genotype = genotype
    cis.list[[genotype]] <- cis
  }
  
  # Calculate CIs for WT (control)
  cis.control = group.CI(get(trait) ~ DAS + Treatment, control, ci=0.84)
  names(cis.control) = c("DAS", "Treatment", "upper", "mean", "lower")
  cis.control$genotype = "WT"
  cis.list[["WT"]] <- cis.control
  
  # Merge all CIs together
  merged = do.call(rbind, cis.list)
  
  p1 <- ggplot(merged, aes(x = DAS, y = mean, group = interaction(Treatment, genotype), color = genotype, linetype = Treatment)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group = interaction(Treatment, genotype), fill = genotype), alpha = 0.2, colour = NA) +
    geom_line(size = .8) +
    theme_minimal() + 
    scale_color_manual(values = genotype.colors) +
    scale_fill_manual(values = genotype.colors) +
    theme(legend.position = "bottom")
  return(p1)
}

pheno = read.csv("data/input/2139DP_IAP_data.csv.lfs.gz")
pheno$DAS = pheno$Day..Int.
# Save WT plants as control to compare to
control = pheno[pheno$Genotype == "WT",]

genotypes = read.csv("data/input/2139DP_mutants.csv")
important.traits = read.csv("data/results/gene_selection_result_roughFixed.csv")
important.traits$gene <- sub("\\.\\d+$", "", important.traits$gene)
important.traits <- important.traits[important.traits$gene %in% genotypes$Targeted.Gene,]

genotype.colors = c("#E41A1C", "#377EB8", "#4DAF4A")  # Red, Blue, Green for distinguishable colors

pdf("figures/SupplementaryFigure9.pdf", 60, 60, onefile = TRUE)
plot.list <- list()

# Add column titles for the genes
# Start with the "column legend" for the top "cell"
p <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Trait / Gene", angle = 0, size = 9, hjust = 0.5)
plot.list[[paste(1, 1, sep = "-")]] <- ggplotGrob(p)
col.index <- 2
for (gene in unique(genotypes$Targeted.Gene)) {
  gene_label <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = gene, size = 9, hjust = 0.5) + 
    theme_void()
  plot.list[[paste("1", col.index, sep = "-")]] <- ggplotGrob(gene_label)
  col.index <- col.index + 1
}

row.index <- 2
for (trait in unique(important.traits$trait)) {
  # Add row title for the trait
  trait_label <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = gsub(" ", ".", paste(strwrap(gsub("\\.", " ", trait), width = 20), collapse = "\n")), angle = 0, size = 9, hjust = 0.5) + 
    theme_void()
  plot.list[[paste(row.index, "1", sep = "-")]] <- ggplotGrob(trait_label)
  col.index <- 2
  for (gene in unique(genotypes$Targeted.Gene)) {
    associated.genotypes <- genotypes$Genotype[genotypes$Targeted.Gene == gene]
    if (gene %in% important.traits$gene[important.traits$trait == trait]) {
      # Generate the plot for the specific trait and genotypes
      p <- ggplotGrob(plot.curve(associated.genotypes, trait))
      plot.list[[paste(row.index, col.index, sep = "-")]] <- p
    } else {
      # Add an empty plot if no association
      p <- ggplot() + theme_void()
      plot.list[[paste(row.index, col.index, sep = "-")]] <- ggplotGrob(p)
    }
    col.index <- col.index + 1
  }
  row.index <- row.index + 1
}

# Arrange the plots in a grid
grid.arrange(
  grobs = plot.list, 
  nrow = length(unique(important.traits$trait)) + 1, 
  ncol = length(unique(genotypes$Targeted.Gene)) + 1
)

# Add a new page and print the genotypes dataframe as a table
grid.newpage()
table_grob <- tableGrob(genotypes, rows = NULL, theme = ttheme_default(base_size = 40))
grid.draw(table_grob)

dev.off()
