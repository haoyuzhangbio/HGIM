############################################################
## Generate intra-person pairwise similarity table
##
## Logic:
##   Same uid = same person
##   Different sampletype from same uid = paired sample types
##   similarity = 1 - bray_distance
##
## Input:
##   ~/Desktop/bray_distance.rds
##   ~/Desktop/ALL11.17.xlsx
##
## Output:
##   ~/Desktop/intra_pd1.csv
############################################################

## ================== 0. Basic settings ==================

setwd("~/Desktop")

bray_file <- "~/Desktop/正刊可用图_20260716/bray_distance_output_0715/bray_distance.rds"
meta_file <- "~/Desktop/20260715.total_bracken/ALL11.17.xlsx"
out_file  <- "~/Desktop/intra_pd1.csv"

if (!file.exists(bray_file)) {
  stop(paste0("Cannot find Bray distance file: ", bray_file))
}

if (!file.exists(meta_file)) {
  stop(paste0("Cannot find metadata file: ", meta_file))
}

## ================== 1. Load packages ==================

packages <- c("readxl", "dplyr", "tidyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(readxl)
library(dplyr)
library(tidyr)

## ================== 2. Read data ==================

bray_distance <- readRDS(bray_file)

meta <- read_excel(
  meta_file,
  col_types = "text"
)

meta <- as.data.frame(meta)

## Check column names
print(colnames(meta))

## ================== 3. Check required columns ==================

required_cols <- c(
  "sample_ID",
  "uid",
  "sampletype"
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "ALL11.17.xlsx is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

## ================== 4. Standardize sample type names ==================

## If your sampletype is already Stool / Gastric fluid / Gastric tissue / Tongue coating,
## this step will not change anything.
## If it is fecal / Gjuice / Gtissue / tongue, it will be converted to the unified name.

meta$sampletype_raw <- meta$sampletype

meta$sampletype <- dplyr::recode(
  meta$sampletype,
  "fecal" = "Stool",
  "stool" = "Stool",
  "Stool" = "Stool",

  "Gjuice" = "Gastric fluid",
  "Gastric fluid" = "Gastric fluid",

  "Gtissue" = "Gastric tissue",
  "Gastric tissue" = "Gastric tissue",

  "tongue" = "Tongue coating",
  "Tongue coating" = "Tongue coating",

  .default = meta$sampletype
)

sample_order <- c(
  "Tongue coating",
  "Gastric fluid",
  "Gastric tissue",
  "Stool"
)

meta <- meta %>%
  filter(sampletype %in% sample_order)

## ================== 5. Keep only samples present in Bray distance matrix ==================

bray_samples <- attr(bray_distance, "Labels")

if (is.null(bray_samples)) {
  bray_samples <- rownames(as.matrix(bray_distance))
}

meta <- meta %>%
  filter(sample_ID %in% bray_samples)

if (nrow(meta) == 0) {
  stop("sample_ID in metadata does not match any sample names in bray_distance.rds. Please check sample_ID.")
}

cat("Samples available for analysis in metadata:", nrow(meta), "\n")
cat("Number of uids in metadata:", length(unique(meta$uid)), "\n")

## ================== 6. Construct pairwise combinations across sample types within each uid ==================

## Keep only individuals with at least 2 sample types
uid_sample_count <- meta %>%
  group_by(uid) %>%
  summarise(
    n_sample = n(),
    n_sampletype = n_distinct(sampletype),
    .groups = "drop"
  )

valid_uid <- uid_sample_count %>%
  filter(n_sampletype >= 2) %>%
  pull(uid)

meta2 <- meta %>%
  filter(uid %in% valid_uid)

cat("Number of uids with at least 2 sample types:", length(unique(meta2$uid)), "\n")

## If the same person has multiple samples of the same sample type,
## keep all of them here. All cross-sampletype pairwise combinations will be generated.

pair_list <- list()
counter <- 1
bray_mat <- as.matrix(bray_distance)

for (u in unique(meta2$uid)) {

  tmp <- meta2 %>%
    filter(uid == u)

  if (nrow(tmp) < 2) next

  cmb <- combn(seq_len(nrow(tmp)), 2)

  for (i in seq_len(ncol(cmb))) {

    a <- cmb[1, i]
    b <- cmb[2, i]

    sample1 <- tmp$sample_ID[a]
    sample2 <- tmp$sample_ID[b]

    type1 <- tmp$sampletype[a]
    type2 <- tmp$sampletype[b]

    ## Only keep cross-sampletype combinations
    if (type1 == type2) next

    ## Bray distance
    bd <- as.numeric(bray_mat[sample1, sample2])

    type_pair <- sort(c(type1, type2))

    ## Use naming consistent with previous code
    type_combination <- dplyr::case_when(
      all(type_pair == sort(c("Gastric fluid", "Tongue coating"))) ~ "Gjuice_tongue",
      all(type_pair == sort(c("Gastric tissue", "Tongue coating"))) ~ "Gtissue_tongue",
      all(type_pair == sort(c("Gastric fluid", "Gastric tissue"))) ~ "Gjuice_Gtissue",
      all(type_pair == sort(c("Stool", "Tongue coating"))) ~ "fecal_tongue",
      all(type_pair == sort(c("Stool", "Gastric fluid"))) ~ "fecal_Gjuice",
      all(type_pair == sort(c("Stool", "Gastric tissue"))) ~ "fecal_Gtissue",
      TRUE ~ paste(type_pair, collapse = "_")
    )

    pair_list[[counter]] <- data.frame(
      uid = u,
      sample_ID_1 = sample1,
      sample_ID_2 = sample2,
      sampletype_1 = type1,
      sampletype_2 = type2,
      type_combination = type_combination,
      bray_distance = bd,
      similarity = 1 - bd,
      stringsAsFactors = FALSE
    )

    ## Attach metadata from the first sample
    ## Usually diagnosis/age/gender etc. are the same within the same uid
    meta_cols_to_add <- setdiff(
      colnames(tmp),
      c("sample_ID", "sampletype", "sampletype_raw")
    )

    for (cc in meta_cols_to_add) {
      pair_list[[counter]][[cc]] <- tmp[[cc]][a]
    }

    counter <- counter + 1
  }
}

intra_pd1 <- bind_rows(pair_list)

if (nrow(intra_pd1) == 0) {
  stop("No cross-sampletype pairs generated within any uid. Please check uid and sampletype.")
}

if ("Age" %in% colnames(intra_pd1)) {
  intra_pd1$age <- as.numeric(intra_pd1$Age)
}

if ("Gender" %in% colnames(intra_pd1)) {
  intra_pd1$gender <- intra_pd1$Gender
}

if ("Smoking" %in% colnames(intra_pd1)) {
  intra_pd1$smoke <- intra_pd1$Smoking
}

if ("Drinking" %in% colnames(intra_pd1)) {
  intra_pd1$drink <- intra_pd1$Drinking
}

if ("BMI" %in% colnames(intra_pd1)) {
  intra_pd1$BMI <- as.numeric(intra_pd1$BMI)
}

## ================== 7. Ensure common columns exist and cast types ==================

## If metadata has "group" but not "diagose", copy the column
if (!"diagose" %in% colnames(intra_pd1) && "group" %in% colnames(intra_pd1)) {
  intra_pd1$diagose <- intra_pd1$group
}

## If group is already in abbreviated form, you may map to full diagnosis names here.
## Not forced here to avoid accidentally modifying your actual grouping.

## Convert numeric columns
numeric_cols <- intersect(
  c("Age", "BMI", "PH", "bray_distance", "similarity"),
  colnames(intra_pd1)
)

for (cc in numeric_cols) {
  intra_pd1[[cc]] <- as.numeric(intra_pd1[[cc]])
}

## ================== 8. Check results ==================

cat("Number of generated pairs:", nrow(intra_pd1), "\n")
cat("type_combination distribution:\n")
print(table(intra_pd1$type_combination, useNA = "always"))

cat("diagose distribution:\n")
if ("diagose" %in% colnames(intra_pd1)) {
  print(table(intra_pd1$diagose, useNA = "always"))
}

cat("similarity range:\n")
print(range(intra_pd1$similarity, na.rm = TRUE))

## ================== 9. Export ==================

write.csv(
  intra_pd1,
  out_file,
  row.names = FALSE
)