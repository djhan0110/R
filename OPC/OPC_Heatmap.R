# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(readxl)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(colorRamp2)

# Load and process data ---------------------------------------------------
OPC_VST <- read_csv("Export/OPC_VST.csv")
OPC_LFC <- read_csv("Export/OPC_LFC.csv")
merged_df <- janitor::clean_names(OPC_VST) %>% 
  left_join(janitor::clean_names(OPC_LFC), by = "x1")

HL_genes <- read_excel("Raw_data/HL_genes.xlsx")
#HL_genes$SYMBOL <- gsub(" ", "", HL_genes$SYMBOL)
annotated_genes <- HL_genes %>%
  mutate(ENSEMBL = mapIds(org.Hs.eg.db, keys = SYMBOL, keytype = "SYMBOL", column = "ENSEMBL"))
rows_with_na_ensembl <- annotated_genes[is.na(annotated_genes$ENSEMBL), ]
print(rows_with_na_ensembl)
syndromic_genes <- annotated_genes %>%
  filter(TYPE == "Syndromic") %>% 
  left_join(merged_df, by = c("ENSEMBL" = "x1"))
syndromic_rows_with_na <- syndromic_genes[is.na(syndromic_genes$padj), ]
print(syndromic_rows_with_na)

nonsyndromic_genes <- annotated_genes %>%
  filter(TYPE == "Non-syndromic") %>% 
  left_join(merged_df, by = c("ENSEMBL" = "x1"))
nonsyndromic_rows_with_na <- nonsyndromic_genes[is.na(nonsyndromic_genes$padj), ]
print(nonsyndromic_rows_with_na)

# Function to process genes and generate heatmap ---------------------------
process_and_generate_heatmap <- function(genes, title, output_file) {
  filtered_df <- merged_df %>% 
    filter(x1 %in% genes$ENSEMBL) %>%
    column_to_rownames("hgnc_symbol.x") %>%
    mutate(attr = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                         ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")))
  
  attr_table <- table(filtered_df$attr)
  selected_matrix <- filtered_df %>% 
    dplyr::select(matches("*.wi_.*|*.8a_.*")) %>% 
    as.matrix()
  
  scaled_matrix <- t(apply(selected_matrix, 1, scale))
  colnames(scaled_matrix) = colnames(selected_matrix)
  
  annotation <- rowAnnotation(
    FC = anno_barplot(filtered_df$log2fold_change,
                      axis_param = list(side = "top",
                                        labels_rot = 0,
                                        gp = gpar(fontface = "bold")),
                      gp = gpar(fill = ifelse(filtered_df$log2fold_change > 1 & filtered_df$padj < 0.05, "firebrick",
                                              ifelse(filtered_df$log2fold_change < -1 & filtered_df$padj < 0.05, "steelblue",
                                                     "white")))), width = unit(2, "cm"),
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontface = "bold"))
  
  heatmap <- Heatmap(selected_matrix, name = "VST", border = TRUE, show_row_dend = FALSE,
                     col = colorRamp2(c(min(selected_matrix), max(selected_matrix)), c("white", "firebrick")), 
                     row_title = title, row_title_gp = gpar(fontface = "bold"),
                     row_names_gp = gpar(fontface = "bold.italic"),
                     column_title_gp = gpar(fontface = "bold"),
                     show_column_names = TRUE, 
                     column_names_side = "bottom",
                     column_names_gp = gpar(fontface = "bold"),
                     column_split = 2, column_title =  c("iPSC", "OPC"),
                     right_annotation = annotation)
  
  tiff(file = output_file, res = 300, width = 6, height = 0.2*nrow(filtered_df), units = "in")
  draw(heatmap)
  dev.off()
}

# Generate heatmaps for non-syndromic and syndromic genes -------------------
process_and_generate_heatmap(nonsyndromic_genes, "Non-syndromic", "Figure/OPC_Heatmap_Non-syndromic.tiff")
process_and_generate_heatmap(syndromic_genes, "Syndromic", "Figure/OPC_Heatmap_Syndromic.tiff")
