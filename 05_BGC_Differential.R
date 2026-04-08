#!/usr/bin/env Rscript
# =============================================================================
# Biosynthetic Gene Clusters (BGCs) Analysis
# Method: Wilcoxon rank-sum tests with Benjamini-Hochberg correction
# =============================================================================

library(dplyr)

bgc_data <- read.csv("bgc_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("poribacteria_metadata.csv", row.names = 1, stringsAsFactors = FALSE)

common <- intersect(rownames(bgc_data), rownames(metadata))
bgc_data <- bgc_data[common, ]
metadata <- metadata[common, , drop = FALSE]

sym_idx <- which(metadata$Group == "Symbiotic")
fl_idx  <- which(metadata$Group == "Free-living")
bgc_classes <- colnames(bgc_data)

results <- data.frame()

for (bgc in bgc_classes) {
  sym_vals <- bgc_data[sym_idx, bgc]
  fl_vals  <- bgc_data[fl_idx, bgc]
  
  if (sd(c(sym_vals, fl_vals)) > 0) {
    test <- wilcox.test(sym_vals, fl_vals, exact = FALSE)
    
    results <- rbind(results, data.frame(
      BGC_Class = bgc,
      p_value = test$p.value,
      Mean_Symbiotic = mean(sym_vals),
      Mean_Free_living = mean(fl_vals)
    ))
  }
}

results$p_adjusted <- p.adjust(results$p_value, method = "BH")
results$significant <- results$p_adjusted < 0.05

write.csv(results, "BGC_Differential_Results.csv", row.names = FALSE)