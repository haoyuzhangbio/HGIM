###############################################
## 1. Load packages
###############################################
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)

###############################################
## 2. Input file paths (modify as needed)
###############################################
#setwd("~/Desktop")
#abun_file <- "2025.11.17.mp4.species.prof.csv"
#meta_file <- "~/Desktop/20260715.total_bracken/ALL11.17.xlsx"
outpath   <- "~/Desktop"

###############################################
## Read 20260715 Bracken relative abundance
###############################################

relative_file <- "~/Desktop/20260715.total_bracken/20260715.total_bracken_profile.relative.txt"

bracken_relative_raw <- read.table(
  relative_file,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

###############################################
## Transpose
###############################################

transpose_bracken <- function(x) {
  sample_ids <- as.character(x[[1]])
  sample_ids <- trimws(sample_ids)
  species_names <- colnames(x)[-1]
  mat <- as.matrix(x[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- sample_ids
  colnames(mat) <- species_names
  mat_t <- t(mat)
  out <- as.data.frame(mat_t, check.names = FALSE)
  return(out)
}

abun <- transpose_bracken(bracken_relative_raw)

###############################################
## Clean sample IDs
###############################################

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", "", x)
  x <- gsub("[[:space:]]+", "", x)
  x <- sub("\\.0$", "", x)
  x <- gsub("\\.", "-", x)
  x <- sub("^W(?=[0-9])", "", x, perl = TRUE)
  x
}

colnames(abun) <- clean_id(colnames(abun))

###############################################
## Read metadata
###############################################

meta_raw <- read_excel(
  "~/Desktop/20260715.total_bracken/ALL11.17.xlsx",
  col_types="text"
)

meta <- meta_raw %>%
  mutate(
    sample = clean_id(sample_ID),
    Tissue = sampletype2
  )

###############################################
## Sample matching
###############################################

common_samples <- intersect(
  colnames(abun),
  meta$sample
)

cat("abun sample count:", ncol(abun), "\n")
cat("meta sample count:", nrow(meta), "\n")
cat("Matched sample count:", length(common_samples), "\n")

###############################################
## Keep matched samples
###############################################

abun <- abun[, common_samples, drop=FALSE]

meta <- meta %>%
  filter(sample %in% common_samples) %>%
  slice(match(common_samples, sample))

stopifnot(
  all(colnames(abun) == meta$sample)
)

###############################################
## 5. Compute mean relative abundance per sampletype2
###############################################

abun_long <- abun %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "abundance") %>%
  left_join(meta %>% select(sample, Tissue), by = "sample")

# Remove samples with no Tissue info
abun_long <- abun_long %>% filter(!is.na(Tissue))

mean_by_type <- abun_long %>%
  group_by(Tissue, gene) %>%
  summarise(mean_abun = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels  = c("Tongue coating", "Gastric fluid", "Gastric tissue", "Stool"),
      ordered = TRUE
    )
  )

###############################################
## 6. Species selection
###############################################

library(readr)

###############################################
## 6.1 Read Gastric tissue specified species list
###############################################

tissue_species <- read.csv(
  "~/Desktop/20260715.total_bracken/braken.tissue.species.list_nopa.csv",
  stringsAsFactors = FALSE
)

# Second column contains species names
gastric_tissue_species <- tissue_species[[2]]

###############################################
## 6.2 Top 15 for other tissues
###############################################

top30_other <- mean_by_type %>%
  filter(
    Tissue != "Gastric tissue"
  ) %>%
  group_by(Tissue) %>%
  arrange(
    desc(mean_abun),
    .by_group = TRUE
  ) %>%
  slice_head(n=15) %>%
  ungroup()

###############################################
## 6.3 keep tissue-associated microbial species
###############################################

gastric_tissue_df <- mean_by_type %>%
  filter(
    Tissue == "Gastric tissue",
    gene %in% gastric_tissue_species
  ) %>%
  arrange(
    desc(mean_abun)
  )

# If more than 30, keep only top 15
gastric_tissue_df <- gastric_tissue_df %>%
  slice_head(n=15)

###############################################
## 6.4 Merge
###############################################

top30_list <- bind_rows(
  top30_other,
  gastric_tissue_df
)

# Fix row order
top30_list$Tissue <- factor(
  top30_list$Tissue,
  levels=c(
    "Tongue coating",
    "Gastric fluid",
    "Gastric tissue",
    "Stool"
  )
)

top30_list <- top30_list %>%
  arrange(Tissue) %>%
  mutate(
    row_id=paste(
      Tissue,
      gene,
      sep="|"
    )
  )

###############################################
## 7. Build heatmap matrix
## Restrict to Gastric tissue only
###############################################

heat.dat_list <- list()

for(i in seq_len(nrow(top30_list))) {
  sp <- top30_list$gene[i]

  # Species raw abundance
  if(sp %in% rownames(abun)) {
    tmp <- abun[
      sp,
      meta$sample
    ]
  } else {
    tmp <- rep(
      0,
      nrow(meta)
    )
    names(tmp) <- meta$sample
  }

  #############################################
  # Process only Gastric tissue region
  #############################################

  tissue_samples <- meta$sample[
    meta$Tissue == "Gastric tissue"
  ]

  # If species is not a gastric tissue specified species,
  # set Gastric tissue region to 0
  if(
    !(sp %in% gastric_tissue_species)
  ) {
    tmp[tissue_samples] <- 0
  }

  heat.dat_list[[i]] <- tmp
}

heat.dat <- do.call(
  rbind,
  heat.dat_list
)

rownames(heat.dat) <- top30_list$row_id

###############################################
## 8. Z-score
###############################################

z.na <- function(x) {
  # Save positions of zeros
  zero_pos <- x == 0

  # Set zeros to NA (not included in mean/sd)
  x[zero_pos] <- NA

  # Log transform (non-zero only)
  x <- log10(x)

  # Compute mean and SD
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)

  # If no valid values or no variation
  if(
    is.na(sd_x) ||
    sd_x == 0
  ) {
    z <- rep(0, length(x))
  } else {
    z <- (x - mean_x) / sd_x
  }

  # Restore original zeros as NA
  z[zero_pos] <- NA

  return(z)
}

mat <- t(
  apply(
    heat.dat,
    1,
    z.na
  )
)

###############################################
## 9. Row split
###############################################

row_split_vec <- top30_list$Tissue

names(row_split_vec) <- top30_list$row_id

###############################################
## 10. Column ordering
###############################################

# Tissue type for each sample
tt <- meta$Tissue[
  match(
    colnames(mat),
    meta$sample
  )
]

tt <- as.character(tt)

# Fixed tissue order
tissue_order <- c(
  "Tongue coating",
  "Gastric fluid",
  "Gastric tissue",
  "Stool"
)

# Order samples
ordered_cols <- unlist(
  lapply(
    tissue_order,
    function(x) {
      cols <- colnames(mat)[
        tt == x
      ]
      cols
    }
  )
)

mat_ordered <- mat[, ordered_cols]

tt_ordered <- tt[
  match(
    colnames(mat_ordered),
    colnames(mat)
  )
]

###############################################
## 11. Column annotation
###############################################

tissue_cols <- c(
  "Tongue coating" = "#299d8f",
  "Gastric fluid"  = "#72bcd5",
  "Gastric tissue" = "#376795",
  "Stool"          = "#e76254"
)

colAnn <- HeatmapAnnotation(
  Tissue = tt_ordered,
  col = list(
    Tissue = tissue_cols
  )
)

###############################################
## 12. Mark top 5
###############################################

top5_list <- top30_list %>%
  group_by(Tissue) %>%
  slice_head(n=5) %>%
  ungroup()

highlight_idx <- match(
  top5_list$row_id,
  rownames(mat_ordered)
)

hi <- rowAnnotation(
  gene = anno_mark(
    at = highlight_idx,
    labels = top5_list$gene,
    labels_gp = gpar(fontsize=7)
  )
)

###############################################
## 13. Color scale
###############################################

heat.color <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#376795", "#BFBFBF", "#e76254")
)

###############################################
## 14. Heatmap
###############################################

ht <- Heatmap(
  mat_ordered,
  name = "Expression",
  col = heat.color,
  na_col = "white",
  cluster_columns = FALSE,
  cluster_rows = FALSE,

  # Split columns by tissue
  column_split = tt_ordered,

  row_split = row_split_vec[rownames(mat_ordered)],

  show_column_names = FALSE,
  show_row_names = FALSE,

  use_raster = TRUE,
  raster_quality = 5,

  border = TRUE,
  top_annotation = colAnn,
  right_annotation = hi
)

###############################################
## 15. Output
###############################################

cairo_pdf(
  "~/Desktop/bodysite_heatmap_Gtissue_selected_species_20260720.pdf",
  width = 10,
  height = 10
)

draw(ht)

dev.off()