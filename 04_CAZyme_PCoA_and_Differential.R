#!/usr/bin/env Rscript
# =============================================================================
# CAZyme Analysis: PCoA and Differential Abundance
# Methods: Bray-Curtis PCoA & Wilcoxon rank-sum test with BH correction
# Thresholds: adjusted p < 0.05, |log2FC| > 1
# =============================================================================

library(dplyr)
library(vegan)
library(ggplot2)

caz_data <- read.csv("cazyme_family_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("poribacteria_metadata.csv", row.names = 1, stringsAsFactors = FALSE)

common <- intersect(rownames(caz_data), rownames(metadata))
caz_data <- caz_data[common, ]
metadata <- metadata[common, , drop = FALSE]

# --- 1. PCoA of CAZyme profiles (Bray-Curtis) ---
caz_rel <- sweep(caz_data, 1, rowSums(caz_data), "/")
caz_rel[is.na(caz_rel)] <- 0

caz_dist <- vegdist(caz_rel, method = "bray")
pcoa_res <- cmdscale(caz_dist, k = 2, eig = TRUE)
eig_pct <- pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100

pcoa_df <- data.frame(PC1 = pcoa_res$points[, 1], PC2 = pcoa_res$points[, 2], Group = metadata$Group)

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Symbiotic" = "#E74C3C", "Free-living" = "#3498DB")) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_pct[1]), y = sprintf("PCoA2 (%.1f%%)", eig_pct[2])) +
  theme_minimal()
ggsave("CAZyme_PCoA_BrayCurtis.png", p_pcoa, width = 6, height = 5)

# --- 2. Differential Abundance (Wilcoxon) ---
sym_idx <- which(metadata$Group == "Symbiotic")
fl_idx  <- which(metadata$Group == "Free-living")

results <- lapply(colnames(caz_data), function(fam) {
  sym_vals <- caz_data[sym_idx, fam]
  fl_vals  <- caz_data[fl_idx, fam]
  
  mean_sym <- mean(sym_vals)
  mean_fl  <- mean(fl_vals)
  lfc <- log2((mean_sym + 0.01) / (mean_fl + 0.01))
  
  if (sd(c(sym_vals, fl_vals)) > 0) {
    test <- wilcox.test(sym_vals, fl_vals, exact = FALSE)
    p_val <- test$p.value
  } else {
    p_val <- NA 
  }
  
  tibble(Family = fam, p_value = p_val, log2FC = lfc,
         mean_Symbiotic = mean_sym, mean_Free_living = mean_fl)
})

res_df <- bind_rows(results) %>% filter(!is.na(p_value))
res_df$p_adjusted <- p.adjust(res_df$p_value, method = "BH")
res_df$significant <- res_df$p_adjusted < 0.05 & abs(res_df$log2FC) > 1

write.csv(res_df, "CAZyme_Wilcoxon_all.csv", row.names = FALSE)
write.csv(res_df %>% filter(significant), "CAZyme_Significant.csv", row.names = FALSE)