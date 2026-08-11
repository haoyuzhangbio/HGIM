############################################################
## Corrected Figure2d biomarker counts
##
## Logic:
## 1. Shared combinations:
##    Species significant in all sources of a combination,
##    with consistent direction
##
## 2. Source only:
##    Species significant in that source,
##    but no other significant source has the same direction
##
## Input:
##   Figure2d.pattern.csv
##
## Output:
##   figure2d_corrected_counts.csv
##   figure2d_corrected_shared_species_list.csv
##   figure2d_corrected_unique_species_list.csv
############################################################

## =========================================================
## 1. Load packages
## =========================================================

library(data.table)
library(dplyr)

## =========================================================
## 2. Read file
## =========================================================

input_file <- "~/Desktop/Fig2_20260720/Figure2d.pattern.kraken.csv"
out_dir <- "~/Desktop"

if (!file.exists(input_file)) {
  stop("Cannot find Figure2d.pattern.csv. Please confirm it is in the same folder as the R script.")
}

pat <- read.csv(
  input_file,
  as.is = TRUE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

## =========================================================
## 3. Check required columns
## =========================================================

required_cols <- c(
  "group2",
  "glmqval_more",
  "source",
  "gene",
  "estimate_more"
)

missing_cols <- setdiff(required_cols, colnames(pat))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Input file missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

## =========================================================
## 4. Extract GC vs nonGC significant differential species
## =========================================================

subdata <- pat[
  pat$group2 == "GC vs nonGC" &
    pat$glmqval_more < 0.05,
]

if (nrow(subdata) == 0) {
  stop("No records with group2 == 'GC vs nonGC' and glmqval_more < 0.05.")
}

subdata$species <- subdata$gene
subdata$group <- subdata$source

## =========================================================
## 5. Fix four body sites
## =========================================================

source_order <- c(
  "Gastric tissue",
  "Gastric fluid",
  "Tongue coating",
  "Stool"
)

source_order <- source_order[source_order %in% unique(subdata$group)]

if (length(source_order) < 2) {
  stop("Fewer than 2 available sources, cannot perform statistics.")
}

## =========================================================
## 6. Keep one direction value per species-source
##    If duplicates exist, take the mean of estimate_more
## =========================================================

species_source <- subdata %>%
  group_by(species, group) %>%
  summarise(
    estimate_more = mean(estimate_more, na.rm = TRUE),
    .groups = "drop"
  )

species_source <- species_source[
  species_source$group %in% source_order,
]

## =========================================================
## 7. Build estimate matrix and sign matrix
## =========================================================

estimate_matrix <- data.table::dcast(
  as.data.table(species_source),
  species ~ group,
  value.var = "estimate_more"
)

estimate_matrix <- as.data.frame(estimate_matrix)
rownames(estimate_matrix) <- estimate_matrix$species
estimate_matrix <- estimate_matrix[, -1, drop = FALSE]

estimate_matrix <- estimate_matrix[
  ,
  source_order,
  drop = FALSE
]

sign_matrix <- sign(estimate_matrix)

present_matrix <- !is.na(estimate_matrix)

## =========================================================
## 8. Function: check if a species has consistent direction
##    across a given set of sources
## =========================================================

is_consistent_for_sources <- function(sign_vec) {
  sign_vec <- sign_vec[!is.na(sign_vec)]
  
  if (length(sign_vec) < 2) {
    return(FALSE)
  }
  
  all(sign_vec > 0) | all(sign_vec < 0)
}

## =========================================================
## 9. Count specified shared combinations
##    Note: inclusive, can be counted multiple times
##    e.g., if GF + TC + Stool is satisfied, it also counts for TC + Stool
## =========================================================

target_intersections <- list(
  c("Gastric fluid", "Tongue coating"),
  c("Gastric tissue", "Tongue coating"),
  c("Stool", "Tongue coating"),
  c("Gastric tissue", "Gastric fluid"),
  c("Gastric fluid", "Stool"),
  c("Gastric tissue", "Stool"),
  c("Gastric fluid", "Tongue coating", "Stool"),
  c("Gastric tissue", "Tongue coating", "Stool")
)

target_intersections <- target_intersections[
  sapply(target_intersections, function(x) all(x %in% colnames(present_matrix)))
]

get_shared_species <- function(sources) {
  keep <- apply(
    present_matrix[, sources, drop = FALSE],
    1,
    function(x) all(x == TRUE)
  )
  
  candidate_species <- rownames(present_matrix)[keep]
  
  candidate_species[
    sapply(candidate_species, function(sp) {
      is_consistent_for_sources(
        as.numeric(sign_matrix[sp, sources])
      )
    })
  ]
}

shared_species_list <- lapply(
  target_intersections,
  get_shared_species
)

names(shared_species_list) <- sapply(
  target_intersections,
  paste,
  collapse = " + "
)

shared_counts <- data.frame(
  item = names(shared_species_list),
  species_count = sapply(shared_species_list, length),
  type = "shared_direction_consistent",
  stringsAsFactors = FALSE
)

## =========================================================
## 10. Count source-specific species
##
## Corrected definition:
## For a given source:
##   The species is significant in this source;
##   AND no other significant source has the same direction;
##   Then count as this source only.
##
## So:
##   If GF and TC are both significant but with opposite directions,
##   the species counts as GF only and TC only separately.
## =========================================================

get_unique_species_for_source <- function(src) {
  all_species <- rownames(present_matrix)
  
  unique_species <- sapply(all_species, function(sp) {
    if (!present_matrix[sp, src]) {
      return(FALSE)
    }
    
    src_sign <- sign_matrix[sp, src]
    
    if (is.na(src_sign) || src_sign == 0) {
      return(FALSE)
    }
    
    other_sources <- setdiff(source_order, src)
    other_sources <- other_sources[present_matrix[sp, other_sources]]
    
    ## Species only significant in this source
    if (length(other_sources) == 0) {
      return(TRUE)
    }
    
    other_signs <- sign_matrix[sp, other_sources]
    
    ## If no other significant source has the same direction, count as source only
    !any(other_signs == src_sign, na.rm = TRUE)
  })
  
  all_species[unique_species]
}

unique_species_list <- lapply(
  source_order,
  get_unique_species_for_source
)

names(unique_species_list) <- paste0(source_order, " only")

unique_counts <- data.frame(
  item = names(unique_species_list),
  species_count = sapply(unique_species_list, length),
  type = "source_specific_or_direction_inconsistent",
  stringsAsFactors = FALSE
)

## =========================================================
## 11. Merge results
## =========================================================

summary_counts <- rbind(
  shared_counts,
  unique_counts
)

## Output in the desired order
desired_order <- c(
  "Gastric fluid + Tongue coating",
  "Gastric tissue + Tongue coating",
  "Stool + Tongue coating",
  "Gastric tissue + Gastric fluid",
  "Gastric fluid + Stool",
  "Gastric tissue + Stool",
  "Gastric fluid + Tongue coating + Stool",
  "Gastric tissue + Tongue coating + Stool",
  "Tongue coating only",
  "Gastric fluid only",
  "Gastric tissue only",
  "Stool only"
)

summary_counts$item <- factor(
  summary_counts$item,
  levels = desired_order
)

summary_counts <- summary_counts[order(summary_counts$item), ]
summary_counts$item <- as.character(summary_counts$item)

print(summary_counts)

write.csv(
  summary_counts,
  file.path(out_dir, "figure2d_corrected_counts.csv"),
  row.names = FALSE
)

## =========================================================
## 12. Export shared combination species lists
## =========================================================

shared_species_long <- do.call(
  rbind,
  lapply(names(shared_species_list), function(nm) {
    data.frame(
      item = nm,
      species = shared_species_list[[nm]],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  shared_species_long,
  file.path(out_dir, "figure2d_corrected_shared_species_list.csv"),
  row.names = FALSE
)

## =========================================================
## 13. Export source-specific species lists
## =========================================================

unique_species_long <- do.call(
  rbind,
  lapply(names(unique_species_list), function(nm) {
    data.frame(
      item = nm,
      species = unique_species_list[[nm]],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  unique_species_long,
  file.path(out_dir, "figure2d_corrected_unique_species_list.csv"),
  row.names = FALSE
)

## =========================================================
## 14. Also output total significant species count per source
## =========================================================

set_size_counts <- data.frame(
  source = source_order,
  significant_species_count = colSums(present_matrix[, source_order, drop = FALSE]),
  stringsAsFactors = FALSE
)

print(set_size_counts)

write.csv(
  set_size_counts,
  file.path(out_dir, "figure2d_all_significant_set_size_counts.csv"),
  row.names = FALSE
)

message("Done!")
message("Output files:")
message("1) figure2d_corrected_counts.csv")
message("2) figure2d_corrected_shared_species_list.csv")
message("3) figure2d_corrected_unique_species_list.csv")
message("4) figure2d_all_significant_set_size_counts.csv")