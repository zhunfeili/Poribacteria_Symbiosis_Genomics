#!/usr/bin/env Rscript
# =============================================================================
# KO Analysis: PCoA and Differential Abundance
# Methods: Bray-Curtis PCoA (vegan) & Negative Binomial Model (DESeq2)
# Thresholds: adjusted p < 0.05, |log2FC| > 1
# =============================================================================

library(dplyr)
library(vegan)
library(DESeq2)
library(ggplot2)

ko_counts <- read.csv("ko_raw_counts.csv", row.names = 1, check.names = FALSE)
metadata  <- read.csv("poribacteria_metadata.csv", row.names = 1, stringsAsFactors = FALSE)

common    <- intersect(rownames(ko_counts), rownames(metadata))
ko_counts <- ko_counts[common, ]
metadata  <- metadata[common, , drop = FALSE]
metadata$Group <- as.factor(metadata$Group)

# --- 1. PCoA of KO profiles (Bray-Curtis) ---
ko_rel <- sweep(ko_counts, 1, rowSums(ko_counts), "/")
ko_rel[is.na(ko_rel)] <- 0

ko_dist <- vegdist(ko_rel, method = "bray")
pcoa_res <- cmdscale(ko_dist, k = 2, eig = TRUE)
eig_pct <- pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100

pcoa_df <- data.frame(PC1 = pcoa_res$points[, 1], PC2 = pcoa_res$points[, 2], Group = metadata$Group)

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Symbiotic" = "#E74C3C", "Free-living" = "#3498DB")) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_pct[1]), y = sprintf("PCoA2 (%.1f%%)", eig_pct[2])) +
  theme_minimal()
ggsave("KO_PCoA_BrayCurtis.png", p_pcoa, width = 6, height = 5)

# --- 2. Differential Abundance (DESeq2) ---
count_matrix <- t(ko_counts)
count_matrix <- round(count_matrix) 

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ Group)
dds <- DESeq(dds, test = "Wald", fitType = "local")

res <- results(dds, contrast = c("Group", "Symbiotic", "Free-living"), alpha = 0.05)
res_df <- as.data.frame(res) %>% filter(!is.na(padj))

res_sig <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(res_df, "KO_DESeq2_all.csv")
write.csv(res_sig, "KO_DESeq2_significant.csv")