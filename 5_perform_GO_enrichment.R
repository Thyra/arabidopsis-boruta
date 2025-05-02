library(topGO)
set.seed(582477)

mothertable <- read.csv("data/input/curated_mothertable20200515.lfs.txt.gz", sep="\t")

all.genes <- unique(as.character(mothertable$target_id))
all.genes <- gsub("\\..*","",all.genes)

# Load Cluster Mapping
clusters = read.csv("data/input/pam_clustering_10cl_050923.txt.gz", sep="\t")
clusters$target.id = sapply(strsplit(as.character(clusters$target_id),"\\."), `[`, 1)

geneID2GO <- readMappings("data/input/ATH_GOslim_04012020_format.lfs.map.gz")

important.genes = read.csv("data/results/gene_selection_result_roughFixed.csv")
important.genes$gene = sapply(strsplit(as.character(important.genes$gene),"\\."), `[`, 1)
network.genes = unique(as.character(important.genes$gene))
  
# Extract Trait Category
important.genes$trait.category = sapply(as.character(important.genes$trait), function(x) { paste(strsplit(x, "\\.")[[1]][2:3], collapse=".") }, USE.NAMES = F)
# For geometry traits, spectrum doesn't matter
important.genes[grepl("geometry.*", important.genes$trait.category),]$trait.category = "geometry"
  
analyze_enrichment <- function(sig_gene_list) {
  GOdata_sig <- new("topGOdata",
                    description = "name",
                    ontology = "BP",
                    allGenes = sig_gene_list,
                    nodeSize = 10,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
  resultsGOfisher <- runTest(GOdata_sig, algorithm = "classic",
                             statistic = "fisher")
  tableGOresults <- GenTable(GOdata_sig,
                             classicFisher=resultsGOfisher,
                             topNodes=length(resultsGOfisher@score))
  
  tableGOresults$classicFisher <- ifelse(tableGOresults$classicFisher == "< 1e-30", 1e-30, tableGOresults$classicFisher)
  tableGOresults$classicFisher <- as.numeric(tableGOresults$classicFisher)
  
  tableGOresults <- tableGOresults[tableGOresults$classicFisher<0.05,]
  
  ret = NULL
  if(nrow(tableGOresults)>0){
    
    sig_list <- names(sig_gene_list)[sig_gene_list == 1]
    tableGOresults$all_ids <- NULL
    tableGOresults$sig_ids <- NA
    
    for (a in 1:length(tableGOresults$GO.ID)){
      all_ids <- genesInTerm(GOdata_sig,tableGOresults[a,"GO.ID"])
      sig_ids <- intersect(all_ids[[1]],sig_list)
      tableGOresults$sig_ids[a] <- paste(sig_ids, collapse=";")
    }
    
    if(length(levels(sig_gene_list))>1){
      # Sort rows by ascending p value
      tableGOresults = tableGOresults[order(tableGOresults$classicFisher),]
      
      ret = tableGOresults
    }
  }
  return(ret)
}
  
upregulated.clusters = c(1,3,4,5,9)
downregulated.clusters = c(2,6,10,7,8)
allgenes = c(upregulated.clusters, downregulated.clusters)
regulation.types = c("upregulated.clusters", "downregulated.clusters")

enrichment.table = NULL
for(regulation.type in regulation.types) {
  xregulated.genes = unique(as.character(clusters[clusters$cluster %in% get(regulation.type),]$target.id))
  
  # ALL GENES IN THE NETWORK
  intersect_Genes <- intersect(unique(as.character(important.genes$gene)), xregulated.genes)
  sig_gene_list <- factor(as.integer(all.genes %in% intersect_Genes))
  
  names(sig_gene_list) <- all.genes
  
  result <- analyze_enrichment(sig_gene_list)
  result$up.or.downregulated = strsplit(regulation.type,"\\.")[[1]][1]
  result$subnetwork = "whole network"
  if(is.null(enrichment.table))
    enrichment.table <- result
  else
    enrichment.table <- rbind(enrichment.table, result)
  
  # Per subnetwork
  for(category in unique(as.character(important.genes$trait.category))) {
    intersect_Genes <- intersect(unique(as.character(important.genes[important.genes$trait.category == category,]$gene)), xregulated.genes)
    sig_gene_list <- factor(as.integer(all.genes %in% intersect_Genes))
    
    names(sig_gene_list) <- all.genes
    
    result <- analyze_enrichment(sig_gene_list)
    result$up.or.downregulated = strsplit(regulation.type,"\\.")[[1]][1]
    result$subnetwork = category
    enrichment.table <- rbind(enrichment.table, result)
  }
}
enrichment.table = enrichment.table[c(8,9,1:7)]
write.csv(enrichment.table, "data/results/GO_term_enrichments.csv", row.names=F)
