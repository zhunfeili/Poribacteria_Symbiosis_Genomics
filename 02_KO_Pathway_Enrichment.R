#!/usr/bin/env Rscript
# =============================================================================
# KO Pathway Enrichment
# Method: Hypergeometric test (clusterProfiler) with BH correction
# =============================================================================

library(clusterProfiler)
library(dplyr)
library(ggplot2)

sig_kos <- read.csv("KO_DESeq2_significant.csv", row.names = 1)
pathway_map <- read.table("ko_pathway_mapping.txt", sep = "\t", col.names = c("KO", "Pathway"))

pathway_map$KO <- sub("^ko:", "", pathway_map$KO)
pathway_map$Pathway <- sub("^path:ko", "", pathway_map$Pathway)

enriched_symb <- rownames(sig_kos %>% filter(log2FoldChange > 1))
enriched_free <- rownames(sig_kos %>% filter(log2FoldChange < -1))
background_kos <- pathway_map$KO 

enrich_kegg <- function(gene_list, group_name) {
  res <- enricher(gene = gene_list, universe = background_kos,
                  TERM2GENE = pathway_map[, c("Pathway", "KO")], 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
  
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    write.csv(as.data.frame(res), paste0("KO_Enrichment_", group_name, ".csv"), row.names = FALSE)
    p <- dotplot(res, showCategory = 15, title = paste(group_name, "Enriched Pathways"))
    ggsave(paste0("KO_Enrichment_", group_name, ".png"), p, width = 7, height = 6)
  }
}

enrich_kegg(enriched_symb, "Symbiotic")
enrich_kegg(enriched_free, "Free-living")