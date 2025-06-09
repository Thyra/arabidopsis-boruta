#May 14 2020
#alternative to kmeans clustering -PAM- is less sentsitive to outliers
getwd()

###from kmeans_051920.R
### change of Nullwerte auf 0.0001 for the genes below
version
#check for Atids (June23/2020)-request from Dennis
#'AT1G05580|AT1G13550|AT1G33840|AT1G52810|AT2G18450|AT2G24610|AT3G07820|AT3G62170|AT5G12000|AT5G42180|AT5G62420'
install.packages("dplyr")
install.packages("scales")
install.packages("rlang")
library(rlang)
library(dplyr)
sessionInfo()

setwd("/prj/ath-dr-rec/PAM_10_GOtermanno_april2723")
dfr <- read.delim("/prj/ath-dr-rec/PAM_10_GOtermanno_april2723/curated_mothertable_08252020.txt")
row.names(dfr) <- dfr$target_id
head(dfr)
colnames(dfr)
dfr2 <- dfr[,c(1,82:94,104:112,114,116,118,121:127)]
row.names(dfr) <- dfr$target_id
colnames(dfr2)
colnames(dfr2)[which(names(dfr2) == "mean_H10d36_D_rec")] <- "36_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H2d17_D")] <- "17_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H3d20_D")] <- "20_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H4d22_D")] <- "22_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H5d24_D")] <- "24_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H6d27_D")] <- "27_DR" 
colnames(dfr2)[which(names(dfr2) == "mean_H7d29_D")] <- "29_DR" 
colnames(dfr2)[which(names(dfr2) == "mean_H8d31_D_rec")] <- "31_DR" 
colnames(dfr2)[which(names(dfr2) == "mean_H9d34_D_rec")] <- "34_DR"
colnames(dfr2)[which(names(dfr2) == "mean_H1d15_WW")] <- "15_WW"
colnames(dfr2)[which(names(dfr2) == "mean_H4d22_WW")] <- "22_WW"
colnames(dfr2)[which(names(dfr2) == "mean_H7d29_WW")] <- "29_WW"
colnames(dfr2)[which(names(dfr2) == "mean_H10d36_WW")] <-"36_WW"
head(dfr2)
colnames(dfr2)

dfr3_to_cluster <-dplyr::filter(dfr2, differentialinatleastoneconstrast == "WAHR")
dim(dfr3_to_cluster)
#17382
colnames(dfr3_to_cluster)
row.names(dfr3_to_cluster) <- dfr3_to_cluster$target_id
head(dfr3_to_cluster)

#check for Atids (June23/2020)-request from Dennis
#check_ats <- dplyr::filter(dfr3_to_cluster, grepl('AT1G05580|AT1G13550|AT1G33840|AT1G52810|AT2G18450|AT2G24610|AT3G07820|AT3G62170|AT5G12000|AT5G42180|AT5G62420',target_id))
#check_ats

#for kmeans we need means columns only
#select for "means" columns

dfr3_to_cluster <- dfr3_to_cluster[,c(2:14)]
head(dfr3_to_cluster)

dfr3_to_cluster[dfr3_to_cluster == 0] <- NA
dfr3_to_cluster <- log2(dfr3_to_cluster)
#remove Nas
#we replace 0 with NA to be able to log transform and remove rows
#with NA values

dfr3_to_cluster <- na.omit(dfr3_to_cluster)
#head(dfr3_to_cluster)

#normalize the data with scale 
dfr3_to_cluster <- t(scale(t(dfr3_to_cluster)))
dim(dfr3_to_cluster)
#[1] 16979    13

#install.packages("cluster")
library(cluster)
library(tidyr)
dim(dfr3_to_cluster)
set.seed(123)
pam<- pam(dfr3_to_cluster,10)
print(pam)
pam$medoids
pam$clusinfo
pam$clustering
pam$objective
pam$isolation
pam$silinfo
pam$diss
pam$id.med
pam$call
#pam$data

# save pam clustering information to a table
write.table(pam$medoids, file="pam10_clustering_medoids_050923_2.txt", quote=F, sep="\t", row.names=T)
write.table(pam$clustering, file="pam_clustering_10cl_050923.txt", quote=F, sep="\t", row.names=T,
            col.names=c("target_id\tcluster"))
write.table(pam$clusinfo, file="pam10_clustinfo_050923.txt", quote=F, sep="\t", row.names=F )
#read.delim("pam_clustering_10cl_050823.txt")

#setwd('/homes/efitzek/Desktop/FDRs/analysis_2020/PAM10_062420_(with 11genes)')
pam <- read.delim("pam_clustering_10cl_050923.txt")
colnames(pam)
head(pam)
head(dfr3_to_cluster)
row.names(dfr3_to_cluster) <- dfr3_to_cluster$target_id

#getwd()
#install.packages("ggtext")
library(ggtext)
library(cowplot)
library(reshape2)
#library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
#install.packages("ggforce")
library(ggforce)

#Orange: #E69F00(Drought stress)
#Hellblau #56B4E9 (well-watered)
setwd('/prj/ath-dr-rec/PAM_10_GOtermanno_april2723/')
medioids <- read.delim("pam10_clustering_medoids_050923.txt")
head(medioids)
colnames(medioids)
#medioids$target_id <- row.names(medioids)
head(medioids)
names(medioids) <- sapply(str_remove_all(colnames(medioids),"X"),"[")
colnames(medioids) <- sub('\\_[^_]+$', '', colnames(medioids))
head(medioids)
med <- melt(medioids)
colnames(med) <- c('sample', 'z.score') 
head(med)
med
pam<- read.delim("pam_clustering_10cl_050923.txt")

cols <- c("15"="#90e8d4",
          "17"="#FED976",
          "20"="#FEB24C",
          "22"="#FD8D3C",
          "24"="#FC4E2A",
          "27"="#E31A1C",
          "29"="#BD0026",
          "31"="#7FCDBB",
          "34"="#41B6C4",
          "36"="#1D91C0")

cluster.plots <- list() 
for (i in 1:10){
  a_cluster <- dfr3_to_cluster[pam$cluster==i,]
  head(a_cluster)
  ww_cluster <- a_cluster[,grep("_WW",colnames(a_cluster))]
  colnames(ww_cluster) <- sub('\\_[^_]+$', '', colnames(ww_cluster))
  head(ww_cluster)
  dr_cluster <- a_cluster[,grep("_DR",colnames(a_cluster))]
  colnames(dr_cluster) <- sub('\\_[^_]+$', '', colnames(dr_cluster))
  ww_cluster <- melt(ww_cluster)
  head(ww_cluster)
  colnames(ww_cluster) <- c('locus', 'sample', 'z.score') 
  head(ww_cluster)
  dr_cluster <- melt(dr_cluster)
  colnames(dr_cluster) <- c('locus', 'sample', 'z.score') 
  head(dr_cluster)
  avg <- dr_cluster%>% group_by(sample) %>% summarise(mean = mean(z.score, na.rm = T))
  head(avg)
    com_cluster.plot <- ggplot() +
    geom_line(alpha=0.1, data=ww_cluster, aes(x=sample, y=z.score,
                                              group=locus), color="#1C67DA")+
    #geom_line(alpha=0.1, data=dr_cluster, aes(x=sample, y=z.score,
     #                                           group=locus), color="#E69F00")+
    #scale_color_gradient(values = c("17"="#FED976","20"="#FEB24C","22"="#FD8D3C","24"="#FC4E2A",
    #                              "27"="#E31A1C","29"="#BD0026","31"="#7FCDBB","34"="#41B6C4","36"="#1D91C0"))+  
    #geom_line(alpha=0.1, data=dr_cluster[ww_cluster$sample >= "15" & ww_cluster$sample <= "22",], 
    #          aes(x=sample, y=z.score, group=locus), color ="#90e8d4")+
    #geom_line(alpha=0.1, data=dr_cluster[ww_cluster$sample >= "29" & ww_cluster$sample <= "36",], 
     #         aes(x=sample, y=z.score, group=locus), color ="#7FCDBB")+
    #geom_line(alpha=0.1, data=dr_cluster[ww_cluster$sample >= "34" & ww_cluster$sample <= "36",], 
    #            aes(x=sample, y=z.score, group=locus), color ="#7FCDBB")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "17" & dr_cluster$sample <= "20",], 
              aes(x=sample, y=z.score, group=locus), color ="#FED976")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "20" & dr_cluster$sample <= "22",], 
                aes(x=sample, y=z.score, group=locus), color ="#FEB24C")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "22" & dr_cluster$sample <= "24",], 
                aes(x=sample, y=z.score, group=locus), color ="#FD8D3C")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "24" & dr_cluster$sample <= "27",], 
                aes(x=sample, y=z.score, group=locus), color ="#FC4E2A")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "27" & dr_cluster$sample <= "29",], 
                aes(x=sample, y=z.score, group=locus), color ="#E31A1C")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "29" & dr_cluster$sample <= "31",], 
                aes(x=sample, y=z.score, group=locus), color ="#7FCDBB")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "31" & dr_cluster$sample <= "34",], 
                aes(x=sample, y=z.score, group=locus), color ="#41B6C4")+
    geom_line(alpha=0.1, data=dr_cluster[dr_cluster$sample >= "34" & dr_cluster$sample <= "36",], 
                aes(x=sample, y=z.score, group=locus), color ="#1D91C0")+
    geom_line(alpha=1, data=avg, aes(x=sample, y=mean),color='black', size = 1.5)+
    #labs(x='days',y='z-score') +
    theme(text = element_text(size=14, face= 'bold'),
          axis.title = element_blank(),
          axis.text = element_text(size = 16, angle=0, vjust=1))
  cluster.plots[[i]] <- com_cluster.plot
}
com_cluster.plot
#save as pdf
pdf("pam10_clusters_060723_version2.pdf", width= 18, height = 10)
print(plot_grid(plotlist=cluster.plots, ncol = 5, nrow = 2))
dev.off()

#GO annotation
source("/prj/ath-dr-rec/PAM_10_GOtermanno_april2723/topGO_alt.GenTable.R")
getwd()

library(topGO)
library(tidyverse)
library(dplyr)
library(readr)

geneID2GO_Ath <- readMappings("ATH_GOslim_2023_04_01_format.map", sep="\t")
head(geneID2GO_Ath)
source("/prj/ath-dr-rec/PAM_10_GOtermanno_april2723/topGO_alt.GenTable.R")
for (i in 1:10) { 
  geneList_pam <- read.delim("pam_clustering_10cl_050823.txt", sep ="\t")
  class(geneList_pam)
  head(geneList_pam)
  geneList_pamc <- factor(as.integer(geneList_pam ["cluster"] == i ))
  head(geneList_pamc)
  total_name_list <- substr(geneList_pam[,"target_id"],1,9)
  class(total_name_list)
  names(geneList_pamc) <- total_name_list
  head(geneList_pamc)
  head(str(geneList_pamc))
  geneList_pamcGO <-new("topGOdata",
                             description = "pam",
                             ontology = "BP",
                             allGenes = geneList_pamc,
                             nodeSize = 10,
                             annot = annFUN.gene2GO,
                             gene2GO = geneID2GO_Ath)
  Fisher_geneList_pamcGO <- runTest(geneList_pamcGO, algorithm = "classic", statistic = "fisher")
  TableGO_geneList_pamc <- topGO_altGenTable(geneList_pamcGO,
                                                  classicFisher = Fisher_geneList_pamcGO,
                                                  topNodes = length(Fisher_geneList_pamcGO@score))
  TableGO_geneList_pamc$q_value <- p.adjust(TableGO_geneList_pamc$classicFisher, method = 'BY')
  sigList <- subset(geneList_pam, cluster == "i")
  sigList <- as.character(sigList [,1])
  sigList <- substr(sigList ,1,9)
  TableGO_geneList_pamc <- TableGO_geneList_pamc[TableGO_geneList_pamc$q_value <0.05,]
  colnames(TableGO_geneList_pamc)
  write.table(TableGO_geneList_pamc, file=paste0("GO_TableGO_geneList_pam10_c", i, "_062420_BP.txt"), sep="\t", row.names=FALSE)
}

