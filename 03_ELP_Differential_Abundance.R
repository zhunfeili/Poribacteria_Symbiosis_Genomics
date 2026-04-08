#!/usr/bin/env Rscript
# =============================================================================
# Eukaryotic-like Proteins (ELPs) Analysis
# Normalization: Relative abundance (raw ELP counts / total Pfam proteins)
# Method: Levene’s test -> Welch’s t-test -> BH correction
# =============================================================================

library(car)
library(dplyr)

elp_counts <- read.csv("elp_raw_counts.csv", row.names = 1, check.names = FALSE)
pfam_totals <- read.csv("total_pfam_counts.csv", row.names = 1)
metadata <- read.csv("poribacteria_metadata.csv", row.names = 1, stringsAsFactors = FALSE)

common <- intersect(intersect(rownames(elp_counts), rownames(pfam_totals)), rownames(metadata))
elp_counts <- elp_counts[common, , drop = FALSE]
pfam_totals <- pfam_totals[common, , drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

# --- Calculate Relative Abundance ---
# Divide raw ELP counts by total Pfam-annotated proteins per genome
elp_rel <- elp_counts / pfam_totals$Total_Pfam

groups <- as.factor(metadata$Group)
domains <- colnames(elp_rel)
results <- data.frame()

# --- Statistical Framework ---
for (domain in domains) {
  vals <- elp_rel[[domain]]
  
  if (sd(vals) > 0) {
    levene_res <- leveneTest(vals ~ groups)
    var_equal <- levene_res$"Pr(>F)"[1] >= 0.05
    
    tt_res <- t.test(vals ~ groups, var.equal = var_equal)
    
    results <- rbind(results, data.frame(
      Domain = domain,
      Levene_p = levene_res$"Pr(>F)"[1],
      T_test_p = tt_res$p.value,
      Mean_Symbiotic = mean(vals[groups == "Symbiotic"]),
      Mean_Free_living = mean(vals[groups == "Free-living"])
    ))
  }
}

results$p_adjusted <- p.adjust(results$T_test_p, method = "BH")
results$significant <- results$p_adjusted < 0.05

write.csv(results, "ELP_Statistical_Results.csv", row.names = FALSE)