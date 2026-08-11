# HGIM

R scripts for the analyses and visualization presented in:

**Oral-gastric-intestinal microbiome reshapes in gastric carcinogenesis**

## Overview

This repository contains the R scripts used for statistical analyses and
visualization of the gastric cancer microbiome study.

The study integrates high-depth metagenomic data from tongue coating,
gastric fluid, gastric tissue and stool samples from 12,876 individuals
spanning healthy stomach, precancerous lesions and gastric cancer.

The analyses characterize microbiome remodeling across oral, gastric and
intestinal niches during gastric carcinogenesis, including microbial
community structure, differential abundance, biomarker identification,
microbial connectivity, strain-level signatures and clinical associations.

## Repository contents

The scripts are organized according to the corresponding figures in the
manuscript.

### Figure 1

- `Fig1d_Venn_diagram.R`
  - Venn diagram analysis of microbial species across body sites.

- `Fig1e_PCoA_visualization.R`
  - Principal coordinates analysis of microbiome community composition.

- `Fig1f_heatmap_visualization.R`
  - Heatmap visualization of microbial community patterns.

- `Fig1h_phylogenetic_tree.R`
  - Phylogenetic analysis and visualization of microbial species.

### Figure 2

- `Fig2a_Shannon_diversity.R`
  - Analysis and visualization of microbial alpha diversity.

- `Fig2b_upset_biomarker_counts.R`
  - UpSet analysis of shared and body-site-specific GC-associated
    microbial biomarkers.

- `Fig2c_proportion_plot.R`
  - Analysis of the proportion of GC-associated microbial species
    emerging across stages of gastric carcinogenesis.

- `Fig2d_relate_differentially abundant microbial species.R`
  - Differential abundance analysis of microbial species across
    gastric carcinogenesis stages and body sites.

- `Fig2d_relate_species_ranking_by_rfe_model.R`
  - Ranking of microbial species based on the RFE-based biomarker model.

### Figure 3

- `Fig3a_intra_person_similarity_compute.R`
  - Computation of within-individual microbiome similarity across
    anatomical sites.

- `Fig3a_intra_similarity_plot_by_group.R`
  - Visualization of within-individual microbiome similarity across
    disease groups.

- `Fig3b_relate_Co-occurrence Co-abundance in GC and NGC and their linking to Gastric pH.R`
  - Analysis of microbial co-occurrence/co-abundance patterns and
    their association with gastric pH.

- `Fig3d_PH_similarity.R`
  - Analysis of microbiome similarity in relation to gastric pH.

### Figure 4

- `Fig4.R`
  - Statistical analyses and visualization for Figure 4.

### Figure 5

- `Fig5.coxh.risk_index.R`
  - Cox proportional hazards analysis and construction of the risk index.

- `Fig5.Forest_plot_HR.R`
  - Visualization of hazard ratios and associated confidence intervals.

- `Fig5.Gastric tissue communities type.R`
  - Analysis and visualization of gastric tissue microbial community types.

- `Fig5.Survival_KM_curves.R`
  - Kaplan–Meier survival analysis and visualization.

## Data availability

The scripts require processed microbiome profiling data, sample metadata
and intermediate analysis results generated during the study.

Due to privacy and data-sharing restrictions, participant-level clinical
information and raw sequencing data are not included in this repository.

The metagenomic sequencing data and associated databases are described in
the Data Availability section of the manuscript.

The scripts contain references to the required input files and their
corresponding analysis steps.

## Software requirements

The analyses were performed using R.

The R scripts use standard R packages for data manipulation, statistical
analysis and visualization. Required packages are loaded within the
corresponding scripts.

## Reproducibility

The scripts in this repository correspond to the analyses and figures
reported in the manuscript. Input data are processed as described in the
Methods section of the manuscript before being supplied to the R scripts.

Because the underlying participant-level data are not publicly distributed,
the repository does not contain a fully self-contained executable workflow
from raw sequencing data to final figures.

## Citation

If you use these scripts or analyses, please cite the associated manuscript:

> Oral-gastric-intestinal microbiome reshapes in gastric carcinogenesis.

## License

This repository is provided for research and academic use. Please refer to
the repository license for details.
