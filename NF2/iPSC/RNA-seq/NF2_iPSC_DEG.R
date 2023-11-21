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
LFC <- clean_names(na.omit(NF2_iPSC_LFC)) %>% 
  mutate(nlogpadj = -log10(padj), attr = ifelse(log2fold_change >= 1 & padj < 0.05, "Up",
                                                ifelse(log2fold_change <= -1 & padj < 0.05, "Down", "ns")))

# Differentially expressed genes ------------------------------------------
sig_DEG <- LFC %>% 
  filter(attr != "ns")
sig_up <- sig_DEG %>% 
  filter(attr == "Up")
sig_up <- sig_up[order(sig_up$log2fold_change, sig_up$nlogpadj, decreasing = T), ]

sig_down <- sig_DEG %>% 
  filter(attr == "Down")
sig_down <- sig_down[order(sig_down$log2fold_change, sig_down$nlogpadj, decreasing = F), ]

topN <- 10
top_DEG <- rbind(head(sig_up, topN), head(sig_down, topN))
top_DEG$attr <- factor(top_DEG$attr, levels = c("Down", "Up"))
top_DEG$hgnc_symbol <- gsub("STON1-GTF2A1L", "STON1-\nGTF2A1L", top_DEG$hgnc_symbol)

plot_DEG <- ggplot(top_DEG, aes(x = log2fold_change, y = reorder(hgnc_symbol, abs(log2fold_change)))) +
  geom_col(aes(fill = attr), position = position_dodge(0.6), width = 0.1, linewidth = 0.2, color = "black", alpha = 1) +
  geom_point(aes(size = nlogpadj, fill = attr), position = position_dodge(0.6), shape = 21, color = "black", stroke = 0.3, alpha = 1) +
  scale_size_continuous(breaks = pretty_breaks(n = 5), range = c(1, 4)) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  facet_wrap(~ attr, scales = "free", ncol = 2, nrow = 1) +
  facetted_pos_scales(x = list(
    attr == "Up" ~ scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                                      breaks = pretty_breaks(n = 3),
                                      labels = label_number(accuracy = 1)),
    attr == "Down" ~ scale_x_continuous(expand = expansion(mult = c(0.2, 0)),
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
    #legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(-1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    #legend.key.size = unit(1, "mm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) 
plot_DEG

ggsave("Figure/DEG_NF2_iPSC.tiff", plot_DEG, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)
