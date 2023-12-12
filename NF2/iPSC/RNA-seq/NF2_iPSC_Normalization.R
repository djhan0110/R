# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

# Load packages -----------------------------------------------------------
library(readr)
library(tidyverse)
library(DESeq2)
library(PCAtools)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggforce)

# Load count matrix -------------------------------------------------------
NF2_iPSC_RC <- read.csv("Raw_data/NF2_iPSC_RC.csv")
raw_data <- NF2_iPSC_RC

clean_df <- janitor::clean_names(na.omit(raw_data)) %>%
  column_to_rownames(var = "x")

# Assuming ctrl_sample is your control sample
ctrl_sample <- "Wi"
cell_line <- factor(rep(c("Wi", "H11"), each = 3)) %>% 
  relevel(ctrl_sample)

# Create sample information data frame
info_df <- data.frame(cell_line)
rownames(info_df) <- colnames(clean_df)
info_df

# Subset the clean_df for a specific gene (ENSG00000186575)
clean_df[rownames(clean_df) %in% "ENSG00000186575", ]

# DESeq2 -------------------------------------------------------------
ddsm <- DESeqDataSetFromMatrix(countData = clean_df,
                               colData = info_df,
                               design = ~ cell_line)

keep_rows <- rowSums(counts(ddsm)) >= ncol(ddsm)
keep_ddsm <- ddsm[keep_rows, ]
dds <- DESeq(keep_ddsm)
vst_dds <- vst(dds, blind = F)
annotation_vst <- as.data.frame(assay(vst_dds)) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(vst_dds), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotation_vst, "Export/NF2_iPSC_VST.csv")

experiment_df <- resultsNames(dds)[2]
result_dds <- results(dds)
annotation_result <- as.data.frame(result_dds) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(result_dds), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_result[annotation_result$hgnc_symbol %in% "NF2", ]
write.csv(annotation_result, "Export/NF2_iPSC_Result.csv")

LFC_dds <- lfcShrink(dds, coef = experiment_df, type = "apeglm")
annotation_LFC <- as.data.frame(LFC_dds) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(LFC_dds), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_LFC[annotation_LFC$hgnc_symbol %in% "NF2", ]
write.csv(annotation_LFC, "Export/NF2_iPSC_LFC.csv")

# Principal component analysis --------------------------------------------
mapping_vst <- annotation_vst

annotation_symbol <- mapIds(org.Hs.eg.db, keys = rownames(mapping_vst),
                            column = "SYMBOL", keytype = "ENSEMBL")
rownames(mapping_vst) <- make.unique(ifelse(is.na(annotation_symbol), rownames(mapping_vst), annotation_symbol))
mapping_vst <- mapping_vst %>% 
  dplyr::select(-hgnc_symbol)

pca_vst <- pca(mapping_vst, metadata = info_df, removeVar = 0.1)
scree_plot <- screeplot(pca_vst, axisLabSize = 10, titleLabSize = 10)
bi_plot <- biplot(pca_vst, showLoadings = TRUE,
                  labSize = 5, pointSize = 5, sizeLoadingsNames = 5,
                  lab = c("wi_i_psc_1"= "Ctrl-1",
                          "wi_i_psc_2"= "Ctrl-2",
                          "wi_i_psc_3"= "Ctrl-3",
                          "h11_i_psc_1"= "KO-1",
                          "h11_i_psc_2"= "KO-2",
                          "h11_i_psc_3"= "KO-3"))
pairs_plot <- pairsplot(pca_vst)
pairs_plot
loadings_plot <- plotloadings(pca_vst, labSize = 2)
loadings_plot

precomp_vst <- prcomp(t(as.matrix(assay(vst_dds))))
pc_eigenvalues <- precomp_vst$sdev^2
eigenvalues_df <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

# Scree plot --------------------------------------------------------------
eigenvalues_df %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# Loadings plot -----------------------------------------------------------
annotation_loadings <- data.frame(precomp_vst$rotation) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(precomp_vst$rotation), keytype = "ENSEMBL", column = "SYMBOL"))
topN = 10
topN_genes <- na.omit(annotation_loadings) %>% 
  dplyr::select(PC1, PC2, hgnc_symbol) %>% 
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>%
  group_by(PC) %>%
  arrange(desc(abs(loading))) %>% 
  dplyr::slice(1:topN) %>% 
  pull(hgnc_symbol) %>%
  unique() 
topN_loadings <- annotation_loadings %>% 
  filter(hgnc_symbol %in% topN_genes)

loadings_plot <- ggplot(topN_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "cm")),
               color = "brown", linewidth = 0.2, alpha = 0.5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = hgnc_symbol),
                  min.segment.length = unit(100, "lines"),
                  force = 2, box.padding = 0.1, max.overlaps = 100,
                  size = 2, fontface = "bold.italic", color = "black") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_size = 8) +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))
loadings_plot
ggsave("Figure/NF2_iPSC_loadings.tiff", loadings_plot, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

# PCA plot ----------------------------------------------------------------
PCA_values <- plotPCA(vst_dds, intgroup = "cell_line", returnData = T)

percentVar <- round(100 * attr(PCA_values, "percentVar"))
pca_plot <- ggplot(PCA_values, aes(PC1, PC2, fill = cell_line)) +
  geom_point(size = 1.5, stroke = 0.3, shape = 21) +
  geom_mark_ellipse(alpha = 0.1, linewidth = NA, expand = unit(1, "mm")) +
  theme_classic(base_size = 8, base_line_size = 0.2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_blank(),
    legend.margin = margin(-2, 0, -5, 0),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(vjust = 2),
    legend.position = "top") +
  scale_fill_manual(labels = c("Wi" = expression(bold(Ctrl^phantom("/"))),
                               "H11" = expression(bold("NF2"^"-/-"))), values = c("#fde725", "#440154")) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  labs(
    x = paste0("PC1: ",percentVar[1],"% variance"),
    y = paste0("PC2: ",percentVar[2],"% variance")) 
pca_plot

ggsave("Figure/PCA_NF2_iPSC.tiff", pca_plot, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)
