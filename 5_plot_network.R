library(stringr)
source("colors.R")

important_genes = read.csv("data/results/gene_selection_result_roughFixed.csv")
# Extract Trait Category
important_genes$trait.category = sapply(as.character(important_genes$trait), function(x) { paste(strsplit(x, "\\.")[[1]][2:3], collapse=".") }, USE.NAMES = F)
# For geometry traits, spectrum doesn't matter
important_genes[grepl("geometry.*", important_genes$trait.category),]$trait.category = "geometry"
important_genes$color <- trait_category_colors[as.character(important_genes$trait.category)]

# Shorten trait names to make them legible
important_genes$label = sapply(as.character(important_genes$trait), function(x) { paste(strsplit(x, "\\.")[[1]][4:length(strsplit(x, "\\.")[[1]])], collapse=".") }, USE.NAMES = F)
important_genes$label = str_remove(important_genes$label, "\\.\\.px(\\.2)?")


sink("data/intermediary/network.dot")
cat("digraph G {\n")
cat('overlap="false"\n')

cat(paste0('"', as.character(unique(important_genes$gene)), '" [label="-" color="#000000bb" fontcolor="#000000bb" fillcolor="#ffffffbb"]\n'))

# cat(paste0('"', ordered.genes, '"', " [fillcolor=orange style=filled]\n"))
cat(paste0('"', important_genes$gene, '" -> ', '"', important_genes$trait, '"[weight=',important_genes$medianImp,', penwidth=',important_genes$medianImp,' color="',trait_category_colors[important_genes$trait.category],'bb" arrowhead=none] \n'))
#cat(paste0('"', important_genes$gene, '" -> ', '"', important_genes$trait.x, '"\n'))

for (category in unique(important_genes$trait.category)) {
  nodes = unique(important_genes[important_genes$trait.category == category,])
  cat(paste0('"', nodes$trait, '" [label="', nodes$label, '" shape=box fontcolor=',ifelse(category == "intensity.fluo", "black", "white"),' fillcolor="', trait_category_colors[category],'" style=filled fontsize=50 width=',nchar(nodes$label)*30/72,']\n'))
}

cat("}")
sink()

# The actual network images I used were done in Cytoscape; just import the dot file using the plugin http://apps.cytoscape.org/apps/dotapp
# and then clicking "prefuse force directed layout"
system("neato -Tpdf data/intermediary/network.dot > plots/network.pdf")

