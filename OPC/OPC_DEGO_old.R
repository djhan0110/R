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

# GOIs --------------------------------------------------------------------
# GOIs <- merged_df %>%
#   filter(grepl("^SLC22", hgnc_symbol.x) | grepl("^TMC", hgnc_symbol.x)) %>% 
#   column_to_rownames("hgnc_symbol.x")
# # 
# # SLC22A1, SLC22A2, SLC22A3, SLC22A4, SLC22A5, SLC22A1, 
# # MET channels LHFPL5 TMIE TMC1 TMC2 PIEZO2
# # TRP channels TRP
# # Endocytosis 
# # 
# # transient receptor potential (TRP) and adenosine triphosphate (ATP) receptor channels
# 
# joined_GOIs <- na.omit(left_join(GOIs, merged_df, by = "hgnc_symbol.y")) 
# selected_GOIs <- GOIs %>%
#   dplyr::select(matches(".*wi.*"), matches(".*8a.*"))
# matrix_GOIs <- as.matrix(selected_GOIs)
# GOIs_column <- c(grep("wi.*", colnames(GOIs), value = TRUE),
#                  grep("8a.*", colnames(GOIs), value = TRUE))
# 
# GOIs_row <- seriate(dist(matrix_GOIs), method = "TSP")
# GOIs_column <- c("i_ps_wi_1", "i_ps_wi_2", "i_ps_wi_3", "i_ps8a_1", "i_ps8a_2", "i_ps8a_3",
#                  "opc_wi_1", "opc_wi_2", "opc_wi_3", "opc8a_1", "opc8a_2", "opc8a_3")
# GOIs_title = ""
# #target_GOIs <- c("SLC22A1", "SLC22A2", "SLC22A3", "SLC22A4", "SLC22A5", "SLC22A16")
# 
# annotation_GOIs <- rowAnnotation(FC = anno_barplot(joined_GOIs$log2fold_change.x,
#                                                    gp = gpar(fill = ifelse(joined_GOIs$log2fold_change.x > 1 & joined_GOIs$padj.x < 0.05, "red",
#                                                                            ifelse(joined_GOIs$log2fold_change.x < -1 & joined_GOIs$padj.x < 0.05, "blue", "gray20")))),
#                                  width = unit(2, "cm")) 
# 
# Heatmap_GOIs <- Heatmap(matrix_GOIs, name = "VST",
#                         show_row_dend = FALSE, row_order = get_order(GOIs_row),
#                         row_names_side = "left", row_names_gp = gpar(fontface = 'bold.italic'),
#                         column_split = factor(rep(c("iPSC", "OPC"), 2, each = 3), levels = c("iPSC", "OPC")),
#                         column_order = GOIs_column,
#                         column_title_gp = gpar(fontface = 'bold'),
#                         left_annotation = annotation_GOIs,
#                         border = TRUE)
# Heatmap_GOIs
# tiff(file = "Figure/OPC_Heatmap_GOIs.tiff", res = 300, width = 5, height = 5, units = "in")
# draw(Heatmap_GOIs)
# dev.off()

# GO ----------------------------------------------------------------------
GO_ids <- c(
  "GO:0006897",
  "GO:0015695",
  "GO:0015867",
  "GO:0030647",
  "GO:0042490")
# GO:0006897 endocytosis
# GO:0015695 organic cation transport
# GO:0015867 ATP transport
# GO:0030647 aminoglycoside antibiotic metabolic process
# GO:0042490 mechanoreceptor differentiation
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

# Split the retrieved_genes data frame by the TERM column
split_GO <- split(distinct_genes, distinct_genes$TERM)
#colname <- split_GO[[i]][5] 
split_GO[[1]][16:27]
generate_heatmap <- function(data) {
  # Create a heatmap for the given data frame
  heatmap(data)
}

heatmaps <- lapply(split_GO, function(df) generate_heatmap(as.matrix(df[, 16:27])))  # Exclude the "TERM" column

















# Initialize an empty list to store the selected matrices
selected_matrix_list <- list()

# Patterns to match for column selection
pattern_wi <- ".*wi.*"
pattern_8a <- ".*8a.*"

# Iterate through the matrices in matrix_list
for (i in seq_along(matrix_list)) {
  # Get the current matrix
  current_matrix <- matrix_list[[i]]
  
  # Find the column indices that match the patterns
  columns_to_select <- which(grepl(pattern_wi, colnames(current_matrix)) | grepl(pattern_8a, colnames(current_matrix)))
  
  # Select the specific columns
  selected_matrix <- current_matrix[, columns_to_select]
  
  # Assign row names based on the "SYMBOL.x" column
  rownames(selected_matrix) <- current_matrix$SYMBOL.x
  
  # Add the selected matrix to the selected_matrix_list
  selected_matrix_list[[i]] <- selected_matrix
}









generate_heatmap <- function(data) {
  # Create a heatmap for the given dataset
  heatmap(data)
}
heatmaps <- lapply(heatmap_data_list, generate_heatmap)
heatmaps <- list()  # Initialize an empty list to store the heatmaps

for (i in seq_along(heatmap_data_list)) {
  heatmap <- generate_heatmap(heatmap_data_list[[i]])
  heatmaps[[i]] <- heatmap
}


matrix_list[[1]][1:18]

for i in x 
matrix_GO <- matrix_list[[i]][dplyr::select(matches(".*wi.*"), matches(".*8a.*")] %>% 
  row.names(matrix_GO, matrix_list[[i]]$SYMBOL.x)
# Initialize an empty list to store the named data frames
GO_list <- list()

matrix_list <- lapply(split_GO, function(df) {
  selected_columns <- df %>%
    dplyr::select(matches(".*wi.*"), matches(".*8a.*"), matches(".*SYMBOL.*"))
  matrix_data <- as.matrix(selected_columns)
  return(matrix_data)
})




# Allocate the subsets to separate data frames
gene_term_dataframes <- lapply(names(gene_term_list), function(term) {
  df <- gene_term_list[[term]]
  df <- df[!duplicated(df$SYMBOL), ]
  df <- df %>% 
    filter(!grepl("ns", attr))
  assign(paste0("GO_", term), df, envir = .GlobalEnv)
  return(df)
})

merged_list <- lapply(names(gene_term_list), function(term) {
  df <- gene_term_list[[term]]
  df <- df[!duplicated(df$SYMBOL), ]
  df <- df %>% 
    filter(!grepl("ns", attr))
  return(df)
})




GOs_row <- seriate(dist(as.matrix(matrix_GOs)), method = "TSP")
GOs_title = "GO:0030647\nAminoglycoside antibiotic metabolic process"

annotation_GOs <- rowAnnotation(FC = anno_barplot(joined_GOs$log2fold_change,
                                                 gp = gpar(fill = ifelse(joined_GOs$log2fold_change > 1 & joined_GOs$padj < 0.05, "red",
                                                                         ifelse(joined_GOs$log2fold_change < -1 & joined_GOs$padj < 0.05, "blue",
                                                                                "gray")))),
                               width = unit(2, "cm")) 

Heatmap_GOs <- Heatmap(as.matrix(matrix_GOs), name = "VSD",
                       show_row_dend = FALSE, row_order = get_order(GOs_row),
                       row_names_side = "left", row_names_gp = gpar(fontface = 'bold.italic'), 
                       column_split = factor(rep(c("iPSC", "OPC"), 2, each = 3), levels = c("iPSC", "OPC")),
                       column_order = GOIs_column,
                       column_title_gp = gpar(fontface = 'bold'),
                       left_annotation = annotation_GOs,
                       row_title = GOs_title, row_title_gp = gpar(col = "black"),
                       border = TRUE)

Heatmap_GOs



tiff(file = "Figure/OPC_Heatmap_GOs.tiff", res = 300, width = 4.3, height = 4.3, units = "in")
draw(Heatmap_GOs)
dev.off()

