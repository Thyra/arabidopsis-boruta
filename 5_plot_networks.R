trait.categories = read.csv("data/input/trait_selection_list.csv")
trait.categories$Variable.ID = make.names(trait.categories$Variable.ID)
trait.categories$Category = trait.categories$imaging.modality
trait.categories[trait.categories$trait.type %in% c("architectural", "biomass-related"),]$Category = "architectural"

# ---- DTW ----
dtw.dist = read.csv("data/results/gene-trait_dtw_distances.csv")
# Only keep those with a distance <= 10
dtw.dist = dtw.dist[dtw.dist$distance <= 10,]
dtw.dist$weight = 100/dtw.dist$distance

trait.categories.dtw = trait.categories[trait.categories$Variable.ID %in% unique(dtw.dist$trait),]

sink("data/results/network_dtw.dot")
cat("digraph G {\n")
cat('overlap="false"\n')

cat(paste0('"', as.character(unique(trait.categories.dtw[trait.categories.dtw$Category == "architectural",]$Variable.ID)), '" [shape=box fillcolor=green style=filled]\n'))
cat(paste0('"', as.character(unique(trait.categories.dtw[trait.categories.dtw$Category == "FLUOR",]$Variable.ID)), '" [shape=box fillcolor=red style=filled]\n'))
cat(paste0('"', as.character(unique(trait.categories.dtw[trait.categories.dtw$Category == "RGB",]$Variable.ID)), '" [shape=box fillcolor=lightblue style=filled]\n'))

cat(paste0('"', dtw.dist$gene, '" -> ', '"', dtw.dist$trait, '"[weight=',dtw.dist$weight,'] \n'))
cat("}")
sink()

# The actual network images I used were done in Cytoscape; just import the dot file using the plugin http://apps.cytoscape.org/apps/dotapp
# and then clicking "prefuse force directed layout"
system("neato -Tpdf data/results/network_dtw.dot > plots/network_dtw.pdf")

