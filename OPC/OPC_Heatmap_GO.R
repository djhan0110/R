# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(ComplexHeatmap)
library(colorRamp2)

# Load and process data ---------------------------------------------------
OPC_VST <- read_csv("Export/OPC_VST.csv")
OPC_LFC <- read_csv("Export/OPC_LFC.csv")
merged_df <- janitor::clean_names(OPC_VST) %>% 
  left_join(janitor::clean_names(OPC_LFC), by = "x1") %>% 
  mutate(attr = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                       ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns"))) %>% 
  na.omit()
colnames(merged_df) <- gsub("i_ps", "ips_", colnames(merged_df))
colnames(merged_df) <- gsub("8a", "_8a", colnames(merged_df))
colnames(merged_df) <- gsub("__", "_", colnames(merged_df))
merged_df

# GOs ---------------------------------------------------------------------
GO_ids <- data.frame(GOID = c("GO:0015695", # GO:0015695 organic cation transport
                              "GO:0015867", # GO:0015867 ATP transport
                               "GO:0030647", # GO:0030647 aminoglycoside antibiotic metabolic process
                               "GO:0042490")) %>% # GO:0006897 endocytosis
  mutate(AnnotationDbi::select(GO.db, keys = GOID, keytype = "GOID",  columns = "TERM")) %>%
  mutate(TERM = paste(toupper(substr(TERM, 1, 1)), substr(TERM, 2, nchar(TERM)), sep = ""))

retrieve_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = GO_ids$GOID, keytype = "GOALL",  columns=c("ENSEMBL", "SYMBOL")) %>% 
  mutate(AnnotationDbi::select(GO.db, keys = GOALL, keytype = "GOID",  columns = "TERM")) %>%
  mutate(TERM = paste(toupper(substr(TERM, 1, 1)), substr(TERM, 2, nchar(TERM)), sep = "")) %>% 
  left_join(merged_df, by = c("ENSEMBL" = "x1")) %>% 
  na.omit()

retrieve_genes$TERM <- factor(retrieve_genes$TERM, levels = GO_ids$TERM)
GO_list <- split(retrieve_genes, retrieve_genes$TERM)

keep_genes <- retrieve_genes %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>% 
  column_to_rownames("SYMBOL") %>% 
  na.omit()

mutated_df <- keep_genes %>% 
  filter(keep_genes$attr != "ns")

for (i in 1:nrow(GO_ids)) {
  mutated_df <- mutated_df %>%
    mutate(!!as.name(as.character(GO_ids[i, 1])) := ifelse(ENSEMBL %in% GO_list[[i]]$ENSEMBL, "Y", "N"))
}

mat <- as.matrix(mutated_df[, grep("wi|8a", colnames(mutated_df))])
mat_GO <- as.matrix(mutated_df[, grep("GO:", colnames(mutated_df))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))
type <- gsub("_\\d+$", "", colnames(mat))

ht_annotation = HeatmapAnnotation(Type = type, show_annotation_name = FALSE,
                       col = list(Type = c("ips_8a" = "#fde725", "ips_wi" = "#35b779", "opc_8a" = "#31688e", "opc_wi" = "#440154")),
                       annotation_legend_param = list(Type = list(direction = "horizontal"),
                                                      labels = c("8A iPSC", "Wi iPSC", "8A OPC", "Wi OPC"),
                                                      labels_gp = gpar(fontface = "bold"),
                                                      nrow = 2))
ht_list <- list()
ht_list <- Heatmap(mat_scaled, name = "Expression", column_km = 2, border = TRUE,
                   col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick")),
                   column_title = c("iPSC", "OPC"), column_title_gp = gpar(fontface = "bold"),
                   top_annotation = ht_annotation, 
                   show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE,
                   heatmap_legend_param = list(direction = "horizontal",
                                               labels_gp = gpar(fontface = "bold")),
                   heatmap_width = unit(4, "cm")) +
  
  rowAnnotation(FC = anno_barplot(mutated_df$log2fold_change,
                                  axis_param = list(side = "top",
                                                    gp = gpar(fontface = "bold", fontsize = 8)),
                                  gp = gpar(
                                    fill = ifelse(mutated_df$log2fold_change > 1 & mutated_df$padj < 0.05, "firebrick",
                                                          ifelse(mutated_df$log2fold_change < -1 & mutated_df$padj < 0.05, "steelblue", "white")))),
                annotation_name_side = "top",
                annotation_name_rot = 90,
                annotation_name_gp = gpar(fontface = "bold")) + 
  
  Heatmap(mat_GO, name = "GO", show_row_dend = FALSE, show_column_dend = FALSE, border = TRUE,
          row_names_gp = gpar(fontface = "bold.italic", fontsize = 10),
          column_names_gp = gpar(fontface = "bold"),
          column_labels = c("OCT", "AT", "AMP", "MD"), column_names_side = "top", column_names_rot = 90,
          rect_gp = gpar(col = "white", lty = 1), col = c("Y" = "gold", "N" = "gray90"),
          heatmap_legend_param = list(direction = "horizontal",
                                      labels_gp = gpar(fontface = "bold")),
          heatmap_width = unit(4, "cm"))

tiff("Figure/OPC_Heatmap_GO.tiff", width = 4, height = 12, units = "in", res = 300)
draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
