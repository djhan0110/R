# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

# library(sva)
# library(DESeq2)
# library(PCAtools)
# library(ggforce)
# library(scales)
# library(ggh4x)

# Load count matrix and sample information --------------------------------
NF2_iPSC_VST <- read_csv("Export/NF2_iPSC_VST.csv")
head(NF2_iPSC_VST)
clean_df <- janitor::clean_names(na.omit(NF2_iPSC_VST)) %>% 
  column_to_rownames(var = "x1") %>% 
  dplyr::select(-hgnc_symbol)

# Assuming ctrl_sample is your control sample
ctrl_sample <- "Wi"
cell_line <- factor(rep(c("Wi", "H11"), each = 3))

# Create sample information data frame
info_df <- data.frame(
  sample = colnames(clean_df),
  cell_line = factor(cell_line, levels = unique(c(ctrl_sample, setdiff(cell_line, ctrl_sample))))
)
info_df

# Principal component analysis --------------------------------------------
precomp_vst <- prcomp(t(as.matrix(clean_df)))
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
pca_vst <- plotPCA(clean_df, intgroup = c("cell_type", "cell_line"), returnData = T)
percentVar_vst <- round(100 * attr(pca_vst, "percentVar"))
pca_vst$group


pc_scores <- data.frame(precomp_vst$x)

pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()



RCs_NF2_iPSC <- read_csv("Raw_data/RCs_NF2_iPSC.csv")
head(RCs_NF2_iPSC)
rc_data <- janitor::clean_names(RCs_NF2_iPSC) %>% 
  column_to_rownames(var = "x1")

cell_line <- c(rep(c("Wi", "H11"), each = 3))
cell_type <- c(rep(c("Ctrl", "KO"), each = 3))
sample_info <- data.frame(sample = colnames(rc_data), cell_line, cell_type)
sample_info$cell_type <- factor(sample_info$cell_type, levels = c("Ctrl", "KO"))
sample_info

# DESeq2 ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = rc_data,
                              colData = sample_info,
                              design = ~ cell_type)
keep <- rowSums(counts(dds)) >=  ncol(dds)
dds <- dds[keep, ]
ddsDE <- DESeq(dds)
vst_counts <- vst(dds, blind = F)


# Define custom fill colors for cell_type
fill_colors <- c("Ctrl" = "#fde725", "KO" = "#440154")

# Create the plot
pca_plot <- ggplot(pca_vst, aes(PC1, PC2, shape = cell_type, fill = cell_type)) +
  geom_point(size = 1.5, stroke = 0.3) +
  geom_mark_ellipse(alpha = 0.25, linewidth = NA, expand = unit(1, "mm")) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  theme_classic(base_size = 8, base_line_size = 0.2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    legend.position = "top",
    legend.margin = margin(-2, 0, -5, 0),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.
  ) +
  scale_fill_manual(values = fill_colors, labels = c("Ctrl" = "Control", "KO" = expression(bold("NF2"^"-/-")))) +
  scale_shape_manual(values = c(21, 24), labels = c("Ctrl" = "Control", "KO" = expression(bold("NF2"^"-/-")))) +
  
  #   
  #   
  # #  guides(fill = guide_legend(override.aes = list(shape = c(21, 24)))) +
  #   scale_shape_manual(labels = c("Ctrl" = "Control", "KO" =)),
  #                      values = c("Ctrl" = 21, "KO" = 24)) +
  #scale_fill_manual(values = fill_colors) +
  labs(
    fill = "",
    shape = "",
    x = "PC1: % variance",
    y = "PC2: % variance") 

pca_plot
ggsave("Figure/NF2_iPSC_PCA.tiff", pca_plot, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)
