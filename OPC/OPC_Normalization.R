# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load required packages --------------------------------------------------
library(readxl)
library(tidyverse)
library(sva)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(PCAtools)
library(ggforce)
library(scales)
library(ggh4x)

# Load raw counts and sample information ----------------------------------
CountMatrix <- read_excel("Raw_data/CountMatrix.xlsx")
rc_data <- as.data.frame(janitor::clean_names(na.omit(CountMatrix)))
head(rc_data)

rownames(rc_data) <- rc_data[["geneid"]]
cln_data <- rc_data[, -which(names(rc_data) == "geneid")]
colnames(cln_data)

column_data <- c("iPS8A-1", "iPS8A-2", "iPS8A-3",
                 "iPSWi-1", "iPSWi-2", "iPSWi-3",
                 "OPC8A-1", "OPC8A-2", "OPC8A-3",
                 "OPCWi-1", "OPCWi-2", "OPCWi-3")
colnames(cln_data) <- column_data
head(cln_data)

cell_line <- c(rep(c("SK-8A", "Wi"), 2, each = 3))
cell_type <- c(rep(c("iPSC", "OPC"), each = 6))
sample_info <- data.frame(sample = colnames(cln_data), cell_line, cell_type)
sample_info$cell_type <- factor(sample_info$cell_type, levels = c("iPSC", "OPC"))
sample_info$cell_line <- factor(sample_info$cell_line, levels = c("Wi", "SK-8A"))
sample_info
rows_with_na <- cln_data[!complete.cases(cln_data), ]
print(rows_with_na)

# DESeq2 ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = cln_data,
                              colData = sample_info,
                              design = ~ cell_type)

keep <- rowSums(counts(dds)) >=  ncol(dds)
dds <- dds[keep, ]
ddsDE <- DESeq(dds)

# Normalization -----------------------------------------------------------
normalized_counts <- as.data.frame(counts(ddsDE, normalized = T))
annotated_counts <- normalized_counts %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(normalized_counts), keytype = "ENSEMBL", column = "SYMBOL"))
head(annotated_counts)
write.csv(annotated_counts, "Export/OPC_NCs.csv")

# Variance Stabilized Data ------------------------------------------------
vsd_counts <- vst(dds, blind = F)
vsd_df <- as.data.frame(assay(vsd_counts))
annotated_vsd <- vsd_df %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(vsd_df), keytype = "ENSEMBL", column = "SYMBOL"))
head(annotated_vsd)
write.csv(annotated_vsd, "Export/OPC_VSD.csv")

# Log2 fold change --------------------------------------------------------
resultsNames(ddsDE)
result_ddsDE <- results(ddsDE, name = resultsNames(ddsDE)[2])
annotated_result <- as.data.frame(result_ddsDE) %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(result_ddsDE), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotated_result, "Export/OPC_Result.csv")

LFC_result <- lfcShrink(ddsDE, coef = resultsNames(ddsDE)[2], type = 'apeglm')
annotation_LFC <- as.data.frame(LFC_result) %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(LFC_result), keytype = "ENSEMBL", column = "SYMBOL"))
write.csv(annotation_LFC, "Export/OPC_LFC.csv")

# PCA plot ----------------------------------------------------------------
vsd_matrix <- t(as.matrix(vsd_df))
precomp_vsd <- prcomp(vsd_matrix)
ev_vsd <- precomp_vsd$sdev^2
ev_vsd <- tibble(PC = factor(1:length(ev_vsd)),
                 variance = ev_vsd) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
ev_vsd %>% ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

loading_data <- data.frame(precomp_vsd$rotation)
annotation_loading <- loading_data %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(loading_data), keytype = "ENSEMBL", column = "SYMBOL"))

convert_tibble <- annotation_loading %>% 
  as_tibble(rownames = 'Geneid')

topN_PCA <- na.omit(convert_tibble) %>% 
  dplyr::select(PC1, PC2, Geneid) %>%   # select only the PCs we are interested in
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% # convert to a "long" format
  group_by(PC) %>% # for each PC
  arrange(desc(abs(loading))) %>% # arrange by descending order of loading
  dplyr::slice(1:10) %>% # take the 10 top rows
  pull(Geneid) %>% # pull the gene column as a vector
  unique() # ensure only unique genes are retained

topN_loadings <- convert_tibble %>% 
  filter(Geneid %in% topN_PCA)

pca_vsd <- plotPCA(vsd_counts, intgroup = c("cell_type", "cell_line"), returnData = T)
percentVar_vsd <- round(100 * attr(pca_vsd, "percentVar"))
pca_vsd$group

loading_pl <- ggplot(data = topN_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "cm")),
               color = "brown", linewidth = 0.2, alpha = 0.8) +
  geom_text_repel(aes(x = PC1, y = PC2, label = hgnc_symbol),
                  min.segment.length = unit(100, 'lines'),
                  force = 10, box.padding = 0.1, max.overlaps = 100,
                  size = 2, fontface = "bold.italic", color = "black") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_size = 8) +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"))
loading_pl
ggsave("Figure/OPC_loading.tiff", loading_pl, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

# Define custom fill colors for cell_type
fill_colors <- c("iPSC" = "#fde725", "OPC" = "#440154")

# Create the plot
pca_plot <- ggplot(pca_vsd, aes(PC1, PC2, shape = cell_line, fill = cell_type)) +
  geom_point(size = 1.5, stroke = 0.3) +
  geom_mark_ellipse(alpha = 0.25, linewidth = NA, expand = unit(1, 'mm')) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21) ),
    shape = guide_legend(override.aes = list(fill = NA))) +
  theme_classic(base_size = 8, base_line_size = 0.2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    legend.position = "top",
    legend.margin = margin(-2, 0, -5, 0),
    legend.key.size = unit(0.1, "cm"),
    legend.spacing.x = unit(0.1, "cm")
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_shape_manual(values = c("Wi" = 21, "SK-8A" = 24)) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  labs(
    fill = "",
    shape = "",
    x = "PC1: % variance",
    y = "PC2: % variance") 

pca_plot

ggsave("Figure/OPC_PCA.tiff", pca_plot, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)
