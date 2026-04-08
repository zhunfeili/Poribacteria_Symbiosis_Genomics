Symbioses-driven phylogenomic and genomic divergence of uncultured Candidatus Poribacteria from its free-living conspecifics 

This repository contains the custom R scripts used for the downstream statistical analyses and data visualization in our study on the evolutionary mechanisms of *Candidatus* Poribacteria. 

## Overview
Our study investigates the genomic trade-offs driving the transition of *Ca.* Poribacteria from free-living marine opportunists to specialized benthic symbionts (e.g., sponges and corals). The custom scripts provided here cover four key dimensions of host integration and evolutionary divergence:
1. Metabolic Pathways (KOs)
2. Carbohydrate-Active Enzymes (CAZymes)
3. Biosynthetic Gene Clusters (BGCs)
4. Eukaryotic-Like Proteins (ELPs)

## Repository Structure
This repository focuses strictly on the custom statistical code. Upstream bioinformatics workflows, software versions, and parameters are detailed in the **Materials and Methods** section of the main manuscript. 

* `poribacteria_metadata.txt`: The core metadata mapping all 166 metagenome-assembled genomes (MAGs) to their respective lifestyles (Symbiotic vs. Free-living).
* `01_KO_PCoA_and_Differential.R`: Script for Bray-Curtis PCoA and negative binomial modeling (DESeq2) of KO profiles.
* `02_KO_Pathway_Enrichment.R`: Script for overall functional pathway enrichment using hypergeometric tests.
* `03_ELP_Differential_Abundance.R`: Script for variance and differential testing of ELP relative abundances.
* `04_CAZyme_PCoA_and_Differential.R`: Script for CAZyme family-level comparisons and PCoA visualization.
* `05_BGC_Differential.R`: Script for class-specific BGC comparisons.

## Data Availability & Usage
To maintain a clean and lightweight codebase, the underlying parsed data matrices (e.g., raw KO counts, CAZyme family counts, ELP abundances, and BGC distributions) are provided as Supplementary Datasets associated with the published manuscript. Raw sequencing data and assembled MAGs are deposited in the NCBI database under BioProject PRJNA1450331.

To reproduce the analyses:
1. Download the corresponding parsed matrices from the manuscript's Supplementary Datasets.
2. Place the datasets in the same working directory as these scripts.
3. If necessary, adjust the input filenames in the `read.csv()` functions within the scripts to match the downloaded Supplementary files.
4. Run the R scripts to generate the statistical outputs and visualizations.

## Dependencies
The scripts are written in R (tested on v4.2+) and require the following core packages:
* `tidyverse` (dplyr, tidyr, ggplot2)
* `vegan`
* `DESeq2`
* `clusterProfiler`
* `car`

## License
This project is licensed under the MIT License - see the LICENSE file for details.
