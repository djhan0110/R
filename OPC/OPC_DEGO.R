# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(org.Hs.eg.db)
library(GO.db)
library(ggplot2)
library(ComplexHeatmap)
library(colorRamp2)
library(seriation)
library(gridExtra)

# Load result -------------------------------------------------------------
OPC_LFC <- read_csv("Export/OPC_LFC.csv") 
clean_LFC <- as.data.frame(janitor::clean_names(na.omit(OPC_LFC))) %>% 
  mutate(nlogpadj = -log10(padj),
         attr = ifelse(padj < 0.05 & log2fold_change > 1, "Up",
                       ifelse(padj < 0.05 & log2fold_change < -1, "Down", "ns")))
filtered_LFC <- clean_LFC %>% 
  filter(attr != "ns")
table(clean_LFC$attr)
table(filtered_LFC$attr)

OPC_VSD <- read_csv("Export/OPC_VSD.csv") 
clean_VSD <- as.data.frame(janitor::clean_names(na.omit(OPC_VSD)))

merged_df <- clean_LFC %>% 
  left_join(clean_VSD, by = "x1")

# GO ----------------------------------------------------------------------
GO_ids <- c(
  "GO:0006897",
  "GO:0015695",
  "GO:0015867",
  "GO:0030647",
  "GO:0042490")

retrieve_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = GO_ids, keytype = "GOALL",  columns=c("ENSEMBL", "SYMBOL")) %>% 
  mutate(AnnotationDbi::select(GO.db, keys = GOALL, keytype = "GOID",  columns = "TERM")) %>% 
  mutate(TERM = paste(toupper(substr(TERM, 1, 1)), substr(TERM, 2, nchar(TERM)), sep = "")) %>% 
  left_join(merged_df, by = c("ENSEMBL" = "x1")) %>% 
  filter(!grepl("ns", attr)) %>% 
  na.omit()

retrieve_genes$TERM <- factor(retrieve_genes$TERM)

distinct_genes <- retrieve_genes %>%
  group_by(TERM) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  ungroup()

heatmap_column <- c(grep("i_ps.*", colnames(distinct_genes), value = TRUE),
                 grep("opc.*", colnames(distinct_genes), value = TRUE))

# Create a list to store separate dataframes and heatmaps
GOID_list <- list()

for (id in GO_ids) {
  # Create a new dataframe and set row names
  mat_df <- as.data.frame(distinct_genes[distinct_genes$GOID == id, c(5, 16:27)])
  rownames(mat_df) <- mat_df$SYMBOL
  
  # Remove the "SYMBOL" column
  mat_df <- mat_df[, -1]
  
  # Assign the dataframe to a variable with a dynamic name
  assign(paste0("mat_", id), as.matrix(mat_df))
}

for (id in GO_ids) {
  seriated_df <- seriate(dist(get(paste0("mat_", id))), method = "TSP")
  assign(paste0("seriated_", id), seriated_df)
}

# Create a list to store heatmaps
Heatmap_list <- list()

for (id in GO_ids) {
  # Generate a heatmap for each GO term and store it in the Heatmap_list
  Heatmap_list[[id]] <- Heatmap(get(paste0("mat_", id)), name = "VST",
                                show_row_dend = FALSE, row_order = get_order(get(paste0("seriated_", id))),
                                row_names_side = "left", row_names_gp = gpar(fontface = 'bold.italic'), 
                                column_split = factor(rep(c("iPSC", "OPC"), each = 6), levels = c("iPSC", "OPC")),
                                column_order = heatmap_column,
                                column_title_gp = gpar(fontface = 'bold'),
                                #left_annotation = annotation_GOs,
                                #row_title = GOs_title, row_title_gp = gpar(col = "black"),
                                border = TRUE)
  
  
  
  
}



# Display or save the heatmaps in a loop
for (id in GO_ids) {
  tiff(
    file = paste0("Figure/OPC_Heatmap_", id, ".tiff"), 
    res = 300, 
    width = 4.3, 
    height = 4.3, 
    units = "in"
  )
  draw(Heatmap_list[[id]])
  dev.off()  # Close the TIFF device after drawing
  
}
