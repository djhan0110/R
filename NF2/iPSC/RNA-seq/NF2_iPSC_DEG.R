# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

# Load libraries ----------------------------------------------------------
library(readr)
library(janitor)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(scales)
library(ggh4x)

# Data processing ---------------------------------------------------------
NF2_iPSC_LFC <- read_csv("Export/NF2_iPSC_LFC.csv")
raw_data <- NF2_iPSC_LFC
clean_LFC <- clean_names(na.omit(raw_data)) %>% 
  mutate(pattern = ifelse(log2fold_change >= 1 & padj < 0.05, "Up",
                          ifelse(log2fold_change <= -1 & padj < 0.05, "Down", "ns")),
         nlogpadj = -log10(padj))

# Differentially expressed genes ------------------------------------------
clean_LFC <- clean_LFC[order(clean_LFC$log2fold_change, clean_LFC$nlogpadj, decreasing = T), ]
topN <- 10
top_DEG <- rbind(head(clean_LFC, topN), tail(clean_LFC, topN))
top_DEG$pattern <- factor(top_DEG$pattern, levels = c("Down", "Up"))
top_DEG$hgnc_symbol <- gsub("STON1-GTF2A1L", "STON1-\nGTF2A1L", top_DEG$hgnc_symbol)

plot_DEG <- ggplot(top_DEG, aes(x = log2fold_change, y = reorder(hgnc_symbol, abs(log2fold_change)))) +
  geom_col(aes(fill = pattern), position = position_dodge(0.6), width = 0.1, linewidth = 0.2, color = "black", alpha = 1) +
  geom_point(aes(size = nlogpadj, fill = pattern), position = position_dodge(0.6), shape = 21, color = "black", stroke = 0.3, alpha = 1) +
  scale_size_continuous(breaks = pretty_breaks(n = 5), range = c(1, 4)) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  facet_wrap(~ pattern, scales = "free", ncol = 2, nrow = 1) +
  facetted_pos_scales(x = list(
    pattern == "Up" ~ scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                                      breaks = pretty_breaks(n = 3),
                                      labels = label_number(accuracy = 1)),
    pattern == "Down" ~ scale_x_continuous(expand = expansion(mult = c(0.2, 0)),
                                        breaks = pretty_breaks(n = 3),
                                        labels = label_number(accuracy = 1)))) +
  labs(x = expression(bold("-log"["2"]*"FC")), y = element_blank(), 
       fill = "", size = expression(bold("-Log(Adj."*bolditalic(p)*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 8, base_line_size = 0.5) +
  theme(
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black", size = 6),
    axis.text.y = element_text(face = "bold.italic", color = "black"),
    legend.position = "top",
    legend.text = element_text(color = "black", size = 6, face = "bold"),
    legend.title = element_text(color = "black", size = 6, face = "bold"),
    legend.margin = margin(-6, 40, -5, 0),
    legend.spacing.x = unit(-1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) 
plot_DEG

ggsave("Figure/DEG_NF2_iPSC.tiff", plot_DEG, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

# -------------------------------------------------------------------------
separate_up <- GO_up %>%
  separate_rows(geneID, sep = "/") %>%
  filter(!is.na(geneID))
unique_up <- separate_up %>%
  distinct(geneID) %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = geneID, keytype = "ENSEMBL", column = "SYMBOL"),
         Family = ifelse(grepl("^[[:alpha:]]+", hgnc_symbol),
                         paste0(regmatches(hgnc_symbol, regexpr("^[[:alpha:]]+", hgnc_symbol)), " Family"),
                         "Unknown Family"))
separate_down <- GO_down %>%
  separate_rows(geneID, sep = "/") %>%
  filter(!is.na(geneID))
unique_down <- separate_down %>%
  distinct(geneID) %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = geneID, keytype = "ENSEMBL", column = "SYMBOL"),
         Family = ifelse(grepl("^[[:alpha:]]+", hgnc_symbol),
                         paste0(regmatches(hgnc_symbol, regexpr("^[[:alpha:]]+", hgnc_symbol)), " Family"),
                         "Unknown Family"))
