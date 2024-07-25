# Load Cluster Mapping
clusters = read.csv("data/input/pam_clustering_10cl_050923.txt.gz", sep="\t")

# Load important genes and categories
important_genes = read.csv("data/results/gene_selection_result_roughFixed.csv")
# Extract Trait Category
important_genes$trait.category = sapply(as.character(important_genes$trait), function(x) { paste(strsplit(x, "\\.")[[1]][2:3], collapse=".") }, USE.NAMES = F)
# For geometry traits, spectrum doesn't matter
important_genes[grepl("geometry.*", important_genes$trait.category),]$trait.category = "geometry"

merged = merge(important_genes, clusters, by.x="gene", by.y="target_id")
ctable = table(merged$trait.category, merged$cluster)
 
# Check if Clusters are enriched in certain subnets
enrichment = data.frame(subnet=character(),cluster=integer(),p.value=double())
for(roow in seq(1, nrow(ctable))) {
  subnet.name = row.names(ctable)[roow]
  for(column in seq(1, ncol(ctable))) {
    cluster.id = colnames(ctable)[column]
    # Matrix is like this:
    #                 Cluster   OtherClusters
    # Category          X           Y
    # OtherCategories   Z           W
    m = matrix(
      c(ctable[roow,column],sum(ctable[,column])-ctable[roow,column],
        sum(ctable[roow,])-ctable[roow,column],sum(ctable)-sum(ctable[,column])-sum(ctable[roow,])),
      nrow=2,ncol=2)
    t = fisher.test(m, alternative="greater")
    enrichment = rbind(enrichment, data.frame(subnet=subnet.name,cluster=as.integer(cluster.id),p.value=t$p.value))

  }
}
enrichment$q.value <- p.adjust(enrichment$p.value, method="BY")
enrichment.sig = enrichment[enrichment$q.value < 0.05,]

# Now let's create a Sankey Diagram of how the clusters in our network are distributed accross categories
# Attention: this has to be saved manually because it can only be exported as HTML (it's interactive)
library(networkD3)
source("colors.R")
nodes = data.frame("name" = c(row.names(ctable), colnames(ctable)))
links.raw = integer()
links.groups = character()
for(roow in seq(1, nrow(ctable))) {
  for(column in seq(1, ncol(ctable))) {
    print(paste(roow, column))
    links.raw = c(links.raw, roow-1, nrow(ctable)+column-1, ctable[roow,column])
    if(nrow(enrichment.sig[enrichment.sig$subnet == row.names(ctable)[roow] & enrichment.sig$cluster == colnames(ctable)[column],]) > 0)
      links.groups = c(links.groups, "sig")
    else
      links.groups = c(links.groups, "nosig")
  }
}
links = as.data.frame(matrix(links.raw,
  byrow = TRUE, ncol = 3))
names(links) = c("source", "target", "value")
links$color = as.factor(links.groups)
my_color <- paste0('d3.scaleOrdinal() .domain(["',paste(nodes$name, collapse='","'),'", "sig", "nosig"]) .range(["',paste(c(unname(trait_category_colors),rep("gray",10)), collapse='","'),'", "#ffd45e", "lightgray"])')
# Rename categories for consistency
nodes$name[1:4] = c("architectural (all)", "fluorescence", "NIR", "VIS")
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name", LinkGroup = "color",
              fontSize= 12, nodeWidth = 30, colourScale=my_color)



