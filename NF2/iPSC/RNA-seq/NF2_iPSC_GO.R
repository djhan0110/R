# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(scales)
library(ggforce)
library(pathview)

# Load result -------------------------------------------------------------
NF2_iPSC_LFC <- read_csv("Export/NF2_iPSC_LFC.csv") 
clean_LFC <- as.data.frame(janitor::clean_names(na.omit(NF2_iPSC_LFC))) %>% 
  mutate(nlogpadj = -log10(padj),
         DEG = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                       ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")))

GO_process <- function(gene_subset, DEG) {
  result <- as.data.frame(enrichGO(gene = gene_subset,
                                   OrgDb = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont = "ALL")) %>%
    mutate(nlogq = -log10(qvalue), DEG = DEG)
  str_sub(result$Description, 1, 1) <- str_sub(result$Description, 1, 1) %>% str_to_upper()
  result$conv_GR <- sapply(result$GeneRatio, function(x) eval(parse(text = x)))
  result <- result[order(result$conv_GR, decreasing = T), ]
  
  return(result)
}

GO_up <- GO_process(clean_LFC$x1[clean_LFC$DEG == "Up"], "Up")
GO_down <- GO_process(clean_LFC$x1[clean_LFC$DEG == "Down"], "Down")

split_up <- split(GO_up, GO_up$ONTOLOGY)
split_down <- split(GO_down, GO_down$ONTOLOGY)

topN <- 3
bind_up <- rbind(head(split_up$BP, topN), head(split_up$CC, topN), head(split_up$MF, topN))
bind_down <- rbind(head(split_down$BP, topN), head(split_down$CC, topN), head(split_down$MF, topN))

bind_GO <- rbind(bind_up, bind_down)
bind_GO$DEG <- factor(bind_GO$DEG, levels = c("Up", "Down"))

# Plotting ----------------------------------------------------------------
plot_GO <- ggplot(bind_GO, aes(x = conv_GR, y = reorder(Description, conv_GR))) +
  geom_col(aes(fill = ONTOLOGY), color = "black", width = 0.1, linewidth = 0.2) +
  geom_point(aes(size = Count, fill = ONTOLOGY), shape = 21, color = "black", stroke = 0.3) +
  facet_col(~ DEG, scales = "free", space = "free", strip.position = "right") +
  scale_fill_ordinal() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  scale_size_continuous(range = c(1, 4), breaks = pretty_breaks(n = 4)) +
  labs(x = "Gene Ratio", y = "", size = expression(bold("-Log("*bolditalic(q)*")")), fill = "Ontology") +
  theme_classic() +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"),
        strip.background = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.position = "top",
        legend.justification = c(0.95, 0),
        legend.background = element_blank())

plot_GO
ggsave("Figure/NF2_iPSC_GO.tiff", plot_GO, units = "in", width = 5, height = 5, device = "tiff", dpi = 300)
