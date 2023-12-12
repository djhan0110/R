# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load required packages --------------------------------------------------
library(readr)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(scales)

# Import and clean data set -----------------------------------------------
OPC_LFC <- read_csv("Export/OPC_LFC.csv")
OPC_LFC <- OPC_LFC %>%
  column_to_rownames("...1") 
cln_LFC <- janitor::clean_names(na.omit(OPC_LFC))
signif_LFC <- cln_LFC %>%
  mutate(signif = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                         ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")),
         nlogpadj = -log10(cln_LFC$padj))
GOI_Volcano <- read_excel("Raw_data/GOI_Volcano.xlsx") %>% 
  mutate(ensg_id = mapIds(org.Hs.eg.db, keys = target, keytype = "SYMBOL", column = "ENSEMBL"))
head(GOI_Volcano)
OPC_Volcano <- ggplot(signif_LFC, aes(x = log2fold_change, y = -log10(padj))) +
  geom_point(shape = 21,
             size = ifelse(rownames(signif_LFC) %in% GOI_Volcano$ensg_id, 1, 0.01),
             fill = ifelse(signif_LFC$signif == "Up", "firebrick",
                           ifelse(signif_LFC$signif == "Down", "steelblue", "gray80")),
             alpha = ifelse(rownames(signif_LFC) %in% GOI_Volcano$ensg_id, 1, 0.1),) +
  geom_text_repel(label = ifelse(rownames(signif_LFC) %in% GOI_Volcano$ensg_id, signif_LFC$hgnc_symbol, NA),
                  color = ifelse(signif_LFC$signif == "Up", "firebrick",
                                 ifelse(signif_LFC$signif == "Down", "steelblue", "gray50")),
                  size = 3, fontface = "bold", direction = "both",
                  bg.color = "white", bg.r = 0.1,
                  segment.size = 0.25, segment.color = "black",
                  force = 0.5, box.padding = 0.5, max.overlaps = 200) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2, alpha = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.2, alpha = 0.5) +
  labs(x = expression(bold("Log"["2"]*"FC")),
       y = expression(paste(bold("-Log(adj.p)")))) +
  ylim(-10, 75) +
  theme_classic(base_size = 8) + 
  theme(
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"))

OPC_Volcano  
         
ggsave("Figure/OPC_Volcano.tiff", OPC_Volcano, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)
