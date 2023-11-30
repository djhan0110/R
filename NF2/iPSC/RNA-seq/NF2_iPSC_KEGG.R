# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

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
NF2_iPSC_LFC <- read.csv("Export/NF2_iPSC_LFC.csv", header = TRUE, row.names = 1)
raw_data <- NF2_iPSC_LFC
clean_LFC <- janitor::clean_names(na.omit(raw_data)) %>% 
  mutate(pattern = ifelse(log2fold_change >= 1 & padj < 0.05, "Up",
                          ifelse(log2fold_change <= -1 & padj < 0.05, "Down", "ns")),
         nlogpadj = -log10(padj))
significant_genes <- clean_LFC %>%
  filter(pattern != "ns")
significant_genes <- significant_genes %>%
  mutate(entrez_ids = mapIds(org.Hs.eg.db, keys = rownames(significant_genes), keytype = "ENSEMBL", column = "ENTREZID"))
enrich_result <- enrichKEGG(gene = significant_genes$entrez_ids,
                            organism = "hsa",
                            keyType = "kegg",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

KEGG_result <- as.data.frame(enrich_result@result)
KEGG_sort <- KEGG_result[order(KEGG_result$qvalue, decreasing = F), ]
head(KEGG_result)
KEGG_topN <- KEGG_sort[1:10, ]
KEGG_topN$GeneRatio_converted <- sapply(strsplit(KEGG_topN$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

# KEGG plot ---------------------------------------------------------------
plot_KEGG <- ggplot(KEGG_topN, aes(x = GeneRatio_converted, y = reorder(Description, -qvalue), fill = qvalue)) +
  geom_point(aes(size = Count), shape = 21, stroke = 0.2) +
  scale_fill_viridis_c(direction = -1)+
  scale_size_continuous(range = c(1, 6), breaks = c(3, 4, 5, 7)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "Gene Ratio", y = "", fill = expression(bolditalic("q") ~ bold("-value"))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"),
        legend.position = "top",
        legend.justification = c(1, 0),
        legend.key.height = unit(2, "mm")) 
plot_KEGG

ggsave("Figure/NF2_iPSC_KEGG.tiff", plot_KEGG, units = "in", width = 5, height = 5, device = "tiff", dpi = 300)

# Retrieve gene and plot --------------------------------------------------
separate_topN <- KEGG_topN %>%
  filter(qvalue < 0.05) %>% 
  separate_rows(geneID, sep = "/") %>%
  filter(!is.na(geneID)) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = geneID, keytype = "ENTREZID", column = "SYMBOL"))
split_topN <- split(separate_topN, separate_topN$Description)

bind_kegg <- data.frame()
for (kegg in names(split_topN)) {
  new_df <- paste0("retrive_", make.names(kegg))
  joined_df <- left_join(split_topN[[kegg]], clean_LFC, by = "hgnc_symbol")
  bind_kegg <- bind_rows(bind_kegg, joined_df) %>% 
    mutate(nlogpadj = -log10(padj))
  bind_kegg <- bind_kegg[order(bind_kegg$log2fold_change, decreasing = F), ]
  }

plot_retrieve <- ggplot(bind_kegg, aes(x = log2fold_change, y = reorder(hgnc_symbol, log2fold_change))) +
  geom_col(aes(fill = pattern), color = "black", width = 0.1, linewidth = 0.2) +
  geom_point(aes(fill = pattern, size = nlogpadj), shape = 21, stroke = 0.2) +
  scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "ns" = "black")) +
  scale_size_continuous(range = c(1, 5), breaks = c(0, -1, -2)) +
  labs(title = "KEGG pathway",
       x = expression(bold("-log"["2"]*"FC")), y = element_blank()) +
  theme_classic(base_size = 12) +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"),
        axis.text.y = element_text(color = "black", face = "bold.italic"),
        legend.position = "")

plot_retrieve
ggsave("Figure/NF2_iPSC_KEGG_DEG.tiff", plot_retrieve, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

# Pathview ----------------------------------------------------------------
gene_data <- significant_genes$log2fold_change
names(gene_data) <- significant_genes$entrez_ids
pathview(gene.data = gene_data, species = "hsa", pathway.id = KEGG_topN$ID[1])
pathview(gene.data = gene_data, species = "hsa", pathway.id = KEGG_topN$ID[1], kegg.native=FALSE)
