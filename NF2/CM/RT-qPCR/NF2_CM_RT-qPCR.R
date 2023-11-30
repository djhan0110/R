# Set working environment
rm(list = ls())
cat("\014")
setwd("C:/Workspace/R/NF2/CM/RT-qPCR")

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)

# Set path
current_path <- getwd()
results_path <- file.path(current_path, "Raw_data")
results_files <- list.files(path = results_path, pattern = ".*\\.(csv)$", full.names = T)

function_read_data <- function(file_path) {
  file_extension <- tools::file_ext(file_path)
  if (file_extension %in% c("csv", "txt")) {
    results <- janitor::clean_names(read.csv(file_path, comment.char = "#"))
  } else {
    stop("Unsupported file format. Please provide a CSV file.")
  }
  return(results)
}

# Create a list to store cleaned data
results_list <- lapply(results_files, function_read_data)
results_list <- results_list[!sapply(results_list, is.null)]
housekeeping_gene <- "GAPDH"

for (i in seq_along(results_list)) {
  file_name <- tools::file_path_sans_ext(basename(results_files[i]))
  assign(file_name, results_list[[i]])
  
  assign(paste0("result_", i), results_list[[i]] %>%
           mutate(cq = ifelse(cq == "Undetermined", 40, as.numeric(cq)),
                  sample = as.character(sample)) %>%
           dplyr::select(sample, target, cq, amp_status, tm1, tm2) %>%
           group_by(sample, target) %>%
           mutate(mean = mean(cq), 
                  sd = sd(cq), 
                  n = n(),
                  zscore = (cq - mean)/sd,
                  outlier = ifelse((abs(zscore) > 1) & (sd > 1) & (n != 1), "Outlier", "F")))
  
  assign(paste0("filtered_", i), get(paste0("result_", i)) %>%
           filter(outlier == "F" & target != housekeeping_gene))
  
  assign(paste0("housekeeping_", i), get(paste0("result_", i)) %>%
           filter(outlier == "F" & target == housekeeping_gene) %>% 
           group_by(sample, target) %>%
           summarise(mean_housekeeping = mean(cq), 
                     sd_housekeeping = sd(cq), 
                     n_housekeeping = n()))
  
  assign(paste0("dct_", i), get(paste0("filtered_", i)) %>%
           left_join(get(paste0("housekeeping_", i)), by = "sample") %>%
           mutate(dct = cq - mean_housekeeping,
                  expression = 2^-dct))
  
  assign(paste0("summary_", i), get(paste0("dct_", i)) %>%
           group_by(sample, target.x) %>% 
           summarise(mean_expression = mean(expression),
                     sd_expression = sd(expression)))
  }

sample_1 <- data.frame(sample = as.character(108:129),
                       cell = c("Ctrl", "KO", "Ctrl", "Ctrl", "Ctrl",
                                "KO", "KO", "KO", "KO", "KO", 
                                "KO", "KO", "KO", "Ctrl", "Ctrl",
                                "Ctrl", "Ctrl", "Ctrl", "Ctrl", "KO", 
                                "KO", "KO"),
                       condition = c("-Dox", "-Dox", "-Dox", "-Dox", "-Dox",
                                     "-Dox", "-Dox", "-Dox", "-Dox", "+Dox",
                                     "-Dox", "-Dox", "-Dox", "-Dox", "-Dox",
                                     "-Dox", "-Dox", "-Dox", "-Dox", "-Dox",
                                     "-Dox", "-Dox"))

stat_1 <- left_join(summary_1, sample_1, by = "sample") %>% 
  filter(sample %in% 124:129, target.x != "NF2")

expression_1 <- left_join(dct_1, sample_1, by = "sample") %>% 
  filter(sample %in% 124:129, target.x != "NF2")

plot_qPCR <- ggplot(stat_1, aes(x = sample, y = mean_expression, fill = cell, group = cell)) +
  geom_col(position = position_dodge(), linewidth = 0.5, color = "black", width = 0.8) +
  geom_jitter(data = expression_1, aes(x = sample, y = expression, shape = cell, group = cell),
              position = position_jitterdodge(0.2), color = "black", size = 1, alpha = 0.8, fill = "white") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression),
                position = position_dodge(0.6), width = 0.3, color = "gray30") +
  facet_wrap(~ target.x, scales = "free", ncol = 3) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = expression(bold("Expression (2"^"-\u0394Ct"~")")), fill = "", shape = "") +
  theme_classic() +
  theme(text = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(face = "bold", color = "black", angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        legend.position = "top") 

plot_qPCR

ggsave("Figure/NF2_CM_RT-qPCR.tiff", plot_qPCR, units = "in", width = 5, height = 5, device = "tiff", dpi = 300)
