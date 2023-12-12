# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load required packages --------------------------------------------------
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(clusterProfiler)
library(KEGGREST)
library(org.Hs.eg.db)
library(ggplot2)
library(pathview)

# Load the RNA-seq count data and set the first column as row names
OPC_LFC <- read.csv("Export/OPC_LFC.csv", header = TRUE, row.names = 1)
signif_genes <- na.omit(OPC_LFC) %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05)
gene_ids <- row.names(signif_genes)
gene_ids <- as.data.frame(ensembl_ids)%>%
  mutate(entrez_ids = mapIds(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "ENTREZID"))

enrich_result <- enrichKEGG(gene = gene_ids$entrez_ids,
                            organism = "hsa",
                            keyType = "kegg",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

KEGG_result <- as.data.frame(enrich_result@result)
KEGG_sort <- KEGG_result[order(KEGG_result$Count, KEGG_result$p.adjust, decreasing = T), ]
head(KEGG_result)
write.csv(KEGG_result, "Export/OPC_KEGG.csv")

KEGG_topN <- KEGG_sort[1:15, ]
KEGG_topN$GeneRatio_converted <- sapply(strsplit(KEGG_topN$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

KEGG_plot <- ggplot(KEGG_topN, aes(x = GeneRatio_converted, y = reorder(Description, Count), fill = p.adjust)) +
  geom_point(aes(size = Count), shape = 21, stroke = 0.2) +
  scale_fill_viridis_c()+
  scale_size_continuous(range = c(1, 4)) 
KEGG_plot
