#!/usr/bin/env Rscript
# =============================================================================
# KO Pathway Enrichment (Overall Potential)
# Method: Hypergeometric test (clusterProfiler) with BH correction
# Description: Evaluates the overall functional pathway enrichment based on the 
#              complete repertoire of KOs present (>0 count) in each lineage.
# =============================================================================

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Load Data ---
ko_counts <- read.csv("ko_raw_counts.csv", row.names = 1, check.names = FALSE)
metadata  <- read.csv("poribacteria_metadata.csv", row.names = 1, stringsAsFactors = FALSE)
pathway_map <- read.table("ko_pathway_mapping.txt", sep = "\t", stringsAsFactors = FALSE, col.names = c("KO", "Pathway"))

# Clean KO and Pathway IDs
pathway_map$KO <- sub("^ko:", "", pathway_map$KO)
pathway_map$Pathway <- sub("^path:ko", "", pathway_map$Pathway)

# Align matrices
common <- intersect(rownames(ko_counts), rownames(metadata))
ko_counts <- ko_counts[common, , drop = FALSE]
metadata  <- metadata[common, , drop = FALSE]

# Merge for grouping
ko_counts$Strain <- rownames(ko_counts)
metadata$Strain <- rownames(metadata)
data_merged <- merge(ko_counts, metadata[, c("Strain", "Group")], by = "Strain")

# --- 2. Prepare Group-Specific KO Lists ---
prepare_gene_list <- function(df) {
  df %>%
    select(-Strain, -Group) %>%
    pivot_longer(cols = everything(), names_to = "KO", values_to = "Count") %>%
    filter(Count > 0) %>%
    distinct(KO) %>%
    pull(KO)
}

symb_kos <- prepare_gene_list(data_merged %>% filter(Group == "Symbiotic"))
free_kos <- prepare_gene_list(data_merged %>% filter(Group == "Free-living"))

# --- 3. Enrichment Analysis ---
background_kos <- unique(pathway_map$KO)

enrich_kegg <- function(ko_list, group_name) {
  res <- enricher(
    gene = ko_list,
    universe = background_kos,
    TERM2GENE = pathway_map[, c("Pathway", "KO")],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5
  )
  
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    # Save table
    res_df <- as.data.frame(res)
    write.csv(res_df, paste0("KO_Enrichment_Overall_", group_name, ".csv"), row.names = FALSE)
    
    # Save plot
    p <- dotplot(res, showCategory = 20, title = paste(group_name, "Overall KO Enrichment"))
    ggsave(paste0("KO_Enrichment_Overall_", group_name, ".png"), p, width = 8, height = 6, dpi = 300)
    
  } else {
    cat(sprintf("No significant pathways found for %s group.\n", group_name))
  }
}

enrich_kegg(symb_kos, "Symbiotic")
enrich_kegg(free_kos, "Free-living")