# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(ggplot2)
library(ggVennDiagram)

# Load and process data ---------------------------------------------------
OPC_VST <- read_csv("Export/OPC_VST.csv")
OPC_LFC <- read_csv("Export/OPC_LFC.csv")
merged_df <- janitor::clean_names(OPC_VST) %>% 
  left_join(janitor::clean_names(OPC_LFC), by = "x1") %>% 
  mutate(attr = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                       ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")))

HL_genes <- read_excel("Raw_data/HL_genes.xlsx")
HL_genes$SYMBOL <- gsub(" ", "", HL_genes$SYMBOL)
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
vd_list <- list(
  syndromic_up = syndromic_genes %>%
    filter(attr != "Down") %>% 
    pull(SYMBOL),
  syndromic_down = syndromic_genes %>%
    filter(attr != "Up") %>% 
    pull(SYMBOL),
  nonsyndromic_up = nonsyndromic_genes %>%
    filter(attr != "Down") %>% 
    pull(SYMBOL),
  nonsyndromic_down = nonsyndromic_genes %>%
    filter(attr != "Up") %>% 
    pull(SYMBOL)
)

# Venn diagram ------------------------------------------------------------
vd_label <- c("Up", "Down")
venn_syndromic <- Venn(vd_list[1:2])
syndromic_pd <- process_data(venn_syndromic, shape_id == '201')

vd_syndromic <- ggplot() +
  
  #scale_fill_gradient(low = 'firebrick', high = 'white') +
  geom_sf(data = venn_region(syndromic_pd), aes(fill = count), show.legend = F) +
  geom_sf(data = venn_setedge(syndromic_pd), color = 'black', linewidth = 0.1, alpha = 1, show.legend = F) +
  geom_sf_label(data = venn_region(syndromic_pd),
                aes(label = paste0(count, '\n(', scales::percent(count/sum(count), accuracy = 2), ')')),
                size = 2, fontface = 'bold', alpha = 0.5) +
  geom_sf_text(data = venn_setlabel(syndromic_pd),
               color = c('black'),
               size = 2, nudge_y = -5, fontface = 'bold', label = vd_label) +
  scale_fill_gradient(low = 'white', high = 'firebrick') +
  labs(title = "Syndromic genes") +
  theme_void() +
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "bold"))
vd_syndromic  
  
venn_setlabel(syndromic_pd)
  
    geom_sf_text(data = venn_setlabel(syndromic_pd),
               color = c('black'),
               size = 2, nudge_y = -5, fontface = 'bold', label = vd_label) 



  scale_fill_gradient(low = 'white', high = 'firebrick') +
  
 
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "bold"))
vdp_up







