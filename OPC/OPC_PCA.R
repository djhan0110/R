# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load library ------------------------------------------------------------
library(readr)
library(tidyverse)

# Load VST ----------------------------------------------------------------
VST_OPC <- read_csv("Export/VST_OPC.csv")
# Set the "...1" column as row names
rownames(VST_OPC) <- VST_OPC$`...1`

# Remove the "hgnc_symbol" and "...1" columns
vsd_df <- VST_OPC %>%
  select(-hgnc_symbol, -`...1`)
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

topN_PCA <-   na.omit(convert_tibble) %>% 
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
               color = "brown", size = 0.2, alpha = 0.8) +
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
#ggsave("Figure/Loading_OPC.tiff", loading_pl, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

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

ggsave("Figure/PCA_OPC.tiff", pca_plot, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)