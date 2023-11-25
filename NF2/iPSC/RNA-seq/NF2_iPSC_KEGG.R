# Set working environment -------------------------------------------------
cat("\014")
rm(list = ls())
setwd("C:/Workspace/R/NF2/iPSC/RNA-seq")

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(scales)
library(ggforce)
library(pathview)