# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(scales)
library(ggh4x)

# Load and process data ---------------------------------------------------
OPC_LFC <- read_csv("Export/OPC_LFC.csv")
cln_LFC <- janitor::clean_names(na.omit(OPC_LFC)) %>% 
  mutate(nlogpadj = -log10(padj),
         attr = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                       ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")),
         labeller = gsub("-", "-\n", hgnc_symbol))
cln_LFC <- cln_LFC[order(cln_LFC$log2fold_change, decreasing = T), ]
table(cln_LFC$attr)

topN = 10
topN_rows <- bind_rows(head(cln_LFC, n = topN), tail(cln_LFC, n = topN))

# Plot --------------------------------------------------------------------
plot_deg <- ggplot(topN_rows, aes(x = log2fold_change, y = reorder(labeller, abs(log2fold_change)), fill = attr)) +  
  geom_col(position = position_dodge(0.6), width = 0.1, linewidth = 0.2, color = "black") +
  geom_point(aes(size = nlogpadj), shape = 21, stroke = 0.2) +
  facet_wrap(~ attr, scales = "free") +
  facetted_pos_scales(x = list(
    attr == "Up" ~ scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                                      breaks = pretty_breaks(n = 3),
                                      labels = label_number(accuracy = 1)),
    attr == "Down" ~ scale_x_continuous(expand = expansion(mult = c(0.2, 0)),
                                        breaks = pretty_breaks(n = 3),
                                        labels = label_number(accuracy = 1)))) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = expression(bold("log"["2"]*"FC")), y = element_blank(), 
       fill = "", size = expression(bold("-Log(Adj."*bolditalic(p)*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # plot.margin = margin(2, 2, 2, 2, "mm"),
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(face = "bold.italic", color = "black"),
    legend.position = "top",
    legend.justification = "center",
    legend.spacing.x = unit(-0.1, "cm"),
    legend.margin = margin(-5, 0, -8, 0),
    strip.text = element_blank()
    ) 

plot_deg

ggsave("Figure/OPC_DEG.tiff", plot_deg, units = "in", width = 3, height = 3, device = "tiff", dpi = 300)

