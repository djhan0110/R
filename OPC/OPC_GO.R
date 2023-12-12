# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/OPC")
current_dir <- getwd()

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(stringr)
library(scales)
library(ggplot2)
library(ggforce)

# Load raw data file ------------------------------------------------------
Result_OPC <- read_csv("Export/Result_OPC.csv")
LFC <- janitor::clean_names(na.omit(Result_OPC)) %>% 
  mutate(nlogpadj = -log10(padj),
         attr = ifelse(log2fold_change > 1 & padj < 0.05, "Up",
                       ifelse(log2fold_change < -1 & padj < 0.05, "Down", "ns")))
LFC <- LFC[order(LFC$log2fold_change, decreasing = T), ]
table(LFC$attr)

genes_up <- LFC[LFC$attr == "Up", "x1"]
genes_down <- LFC[LFC$attr == "Down", "x1"]

# GO analysis -------------------------------------------------------------
# GO_up <- as.data.frame(enrichGO(gene = genes_up$x1,
#                                 OrgDb = "org.Hs.eg.db",
#                                 keyType = "ENSEMBL",
#                                 ont = "ALL")) %>% 
#   mutate(nlogq = -log10(qvalue), attr = "Up")
# table(GO_up$ONTOLOGY)
# 
# str_sub(GO_up$Description, 1, 1) <- str_sub(GO_up$Description, 1, 1) %>% str_to_upper()
# GO_up$conv_GR <- sapply(GO_up$GeneRatio, function(x) eval(parse(text = x)))
# write.csv(GO_up, "Export/GO_OPC_up.csv")
# 
# GO_down <- as.data.frame(enrichGO(gene = genes_down$x1,
#                                 OrgDb = "org.Hs.eg.db",
#                                 keyType = "ENSEMBL",
#                                 ont = "ALL")) %>%
#   mutate(nlogq = -log10(qvalue), attr = "Down")
# table(GO_down$ONTOLOGY)
# 
# str_sub(GO_down$Description, 1, 1) <- str_sub(GO_down$Description, 1, 1) %>% str_to_upper()
# GO_down$conv_GR <- sapply(GO_down$GeneRatio, function(x) eval(parse(text = x)))
# write.csv(GO_down, "Export/GO_OPC_down.csv")

GO_up <- read_csv("Export/GO_OPC_up.csv")
GO_up <- GO_up[order(GO_up$conv_GR, decreasing = T), ]
GO_down <- read_csv("Export/GO_OPC_down.csv")
GO_down <- GO_down[order(GO_down$conv_GR, decreasing = T), ]
table(GO_up$ONTOLOGY)
table(GO_down$ONTOLOGY)
split_up <- split(GO_up, GO_up$ONTOLOGY)
split_down <- split(GO_down, GO_down$ONTOLOGY)

topN <- 5
bind_up <- rbind(head(split_up$BP, topN), head(split_up$CC, topN), head(split_up$MF, topN))
bind_down <- rbind(head(split_down$BP, topN), head(split_down$CC, topN), head(split_down$MF, topN))

bind_GO <- rbind(bind_up, bind_down)
bind_GO$attr <- factor(bind_GO$attr, levels = c("Up", "Down"))

# Plotting ----------------------------------------------------------------
pl_GO <- ggplot(bind_GO, aes(x = conv_GR, y = reorder(Description, conv_GR))) +
  geom_col(aes(fill = ONTOLOGY), color = "black", width = 0.1, linewidth = 0.2) +
  geom_point(aes(size = Count, fill = ONTOLOGY), shape = 21, color = "black", stroke = 0.3) +
  facet_col(~ attr, scales = "free", space = "free", strip.position = "right") +
#scale_fill_manual(values = c("BP" = "#d95f02", "MF" = "#1b9e77", "CC" = "#7570b3")) +
  scale_fill_ordinal() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                     labels = label_number(accuracy = 0.01),
                     breaks = pretty_breaks(n = 3)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  scale_size_continuous(breaks = pretty_breaks(n = 3), range = c(0.5, 3)) +
  
  labs(x = "Gene Ratio", y = "", size = expression(bold("-Log("*bolditalic(q)*")")), fill = "Ontology") +
  theme_classic(base_size = 6, base_line_size = 0.5) +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text = element_text(size = 6, color = "black", face = "bold"),
        strip.background = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.position = "top",
        legend.justification = c(0.95, 0),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.2, "mm"),
        legend.margin = margin(0, 10, 0, 10))
pl_GO

ggsave("Figure/OPC_GO.tiff", pl_GO, units = "in", width = 3, height = 6, device = "tiff", dpi = 300)

