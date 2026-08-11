#############################
## 1. Basic settings & load packages
#############################

# Input file paths
abun_file <- "~/Desktop/20260715.total_bracken/20260715.total_bracken_profile.relative.txt"
meta_file <- "~/Desktop/20260715.total_bracken/ALL11.17.xlsx"

# Output directory and prefix
setwd("~/Desktop")
out_dir    <- "~/Desktop"
out_prefix <- "Fig.2A.shannon.2026.07.15"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load required R packages
library(vegan)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(Rmisc)
library(ggbreak)

#############################
## 2. ID cleaning function
#############################

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

#############################
## 3. Read abundance table
#############################

# File format:
# First column = samples
# Subsequent columns = species
df_raw <- read.table(
  abun_file,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

cat("Abundance table raw dimensions:\n")
print(dim(df_raw))

cat("First 10 column names:\n")
print(colnames(df_raw)[1:min(10, ncol(df_raw))])

# First column is the sample name column
sample_col <- colnames(df_raw)[1]

sample_ids <- clean_id(df_raw[[sample_col]])
species_names <- colnames(df_raw)[-1]

# Extract numeric matrix: rows = samples, columns = species
mat <- as.matrix(df_raw[, -1, drop = FALSE])
mode(mat) <- "numeric"

rownames(mat) <- sample_ids
colnames(mat) <- species_names

# Transpose to format needed by subsequent code:
# rows = species, columns = samples
mat_t <- t(mat)

df <- as.data.frame(
  mat_t,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Transposed abundance matrix dimensions:\n")
print(dim(df))

cat("First 10 sample names after transpose:\n")
print(colnames(df)[1:min(10, ncol(df))])

#############################
## 4. Read metadata
#############################

meta_raw <- read_excel(
  meta_file,
  col_types = "text"
)

meta_raw <- as.data.frame(meta_raw)

if (!all(c("sample_ID", "sampletype2", "group") %in% colnames(meta_raw))) {
  stop("meta_file must contain columns: sample_ID, sampletype2, group")
}

# Clean metadata sample_ID
meta_raw$sample_ID_raw <- as.character(meta_raw$sample_ID)
meta_raw$sample_ID <- clean_id(meta_raw$sample_ID)

#############################
## 5. Match samples
#############################

common_samples <- intersect(colnames(df), meta_raw$sample_ID)

cat("Abundance table sample count:", ncol(df), "\n")
cat("Metadata sample count:", nrow(meta_raw), "\n")
cat("Common sample count:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  
  cat("\nFirst 20 sample names in abundance table:\n")
  print(head(colnames(df), 20))
  
  cat("\nFirst 20 sample_ID in metadata:\n")
  print(head(meta_raw$sample_ID, 20))
  
  cat("\nFirst 20 raw sample_ID in metadata:\n")
  print(head(meta_raw$sample_ID_raw, 20))
  
  stop("No overlapping sample_ID between abundance table and metadata. Check ID formats above.")
}

# Align by common_samples
df <- df[, common_samples, drop = FALSE]

meta <- meta_raw[
  match(common_samples, meta_raw$sample_ID),
  ,
  drop = FALSE
]

stopifnot(all(colnames(df) == meta$sample_ID))

cat("\nSample matching successful!\n")
cat("Aligned abundance matrix dimensions:\n")
print(dim(df))

# Subset and reorder based on common samples
df_sub   <- df[, common_samples, drop = FALSE]
meta_sub <- meta_raw[match(common_samples, meta_raw$sample_ID), ]

# Confirm order consistency
stopifnot(identical(colnames(df_sub), meta_sub$sample_ID))

#############################
## 5. Calculate Shannon and richness
#############################

# Note: vegan::diversity requires rows = samples, columns = species, so transpose
shannon_vec <- diversity(t(df_sub), index = "shannon")

# Richness: count of species > 0
richness_vec <- apply(t(df_sub), 1, function(x) sum(x > 0))

# Combine into a data frame
dat <- meta_sub %>%
  mutate(
    Shannon  = shannon_vec[sample_ID],
    richness = richness_vec[sample_ID]
  )

#############################
## 6. Organize grouping info: sampletype2 & diagose1
#############################

# Simple replacement function
replacex <- function(x, org, value) {
  for (ii in seq_along(org)) {
    x[x == org[ii]] <- value[ii]
  }
  x
}

## 6.1 Standardize sample types (4 categories)
dat$sampletype2 <- as.character(dat$sampletype2)
dat$sampletype2 <- replacex(
  dat$sampletype2,
  c("Tongue", "tongue", "Tongue coating",
    "Gjuice", "gastric fluid", "Gastric fluid",
    "Gtissue", "GPtissue", "Gastric tissue",
    "fecal", "Fecal", "Stool"),
  c("Tongue coating", "Tongue coating", "Tongue coating",
    "Gastric fluid", "Gastric fluid", "Gastric fluid",
    "Gastric tissue", "Gastric tissue", "Gastric tissue",
    "Stool", "Stool", "Stool")
)

dat$sampletype2 <- factor(
  dat$sampletype2,
  levels = c("Tongue coating", "Gastric fluid", "Gastric tissue", "Stool")
)

## 6.2 Map group codes to clinical diagnosis names (diagose1)
# group: c("GP","GU","HS","SG","AG","IM","IN","GC")
dat$diagose1 <- as.character(dat$group)

# Map to full names
dat$diagose1 <- replacex(
  dat$diagose1,
  c("GP",               "GU",             "HS","SG",
    "AG",                   "IM",                   "IN",
    "GC"),
  c("Gastric polyps",   "Gastric ulcer",  "Healthy stomach","Superficial gastritis",
    "Atrophic gastritis",                "Intestinal metaplasia","Intraepithelial neoplasia",
    "Gastric cancer")
)

# Set factor order
dat$diagose1 <- factor(
  dat$diagose1,
  levels = c(
    "Gastric polyps",
    "Gastric ulcer",
    "Healthy stomach",
    "Superficial gastritis",
    "Atrophic gastritis",
    "Intestinal metaplasia",
    "Intraepithelial neoplasia",
    "Gastric cancer"
  )
)

#############################
## 7. Compute mean ± SE (Shannon)
#############################

SE_shan <- summarySE(
  dat,
  measurevar = "Shannon",
  groupvars  = c("diagose1", "sampletype2")
)

SE_shan$type  <- "Shannon"
colnames(SE_shan)[which(colnames(SE_shan) == "Shannon")] <- "value"

SE <- SE_shan

#############################
## 8. Handle line connections (optional)
#############################

# e.g., disconnect Gastric polyps and Gastric ulcer from the line
SE_filtered <- SE %>%
  mutate(
    value_line = ifelse(
      diagose1 %in% c("Gastric polyps", "Gastric ulcer"),
      NA,
      value
    )
  )

#############################
## 9. Wilcoxon test and significance stars
#############################

get_significance <- function(dat_all,
                             measure   = "Shannon",
                             ref_group = "Healthy stomach") {
  diagnoses   <- levels(dat_all$diagose1)
  sampletype2s <- levels(dat_all$sampletype2)
  
  res <- data.frame(
    diagose1    = character(),
    sampletype2 = character(),
    label      = character(),
    p_val      = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (st in sampletype2s) {
    ref_values <- dat_all[[measure]][
      dat_all$diagose1 == ref_group & dat_all$sampletype2 == st
    ]
    
    for (dg in diagnoses) {
      if (dg == ref_group) next
      
      test_values <- dat_all[[measure]][
        dat_all$diagose1 == dg & dat_all$sampletype2 == st
      ]
      
      if (length(test_values) > 0 && length(ref_values) > 0) {
        p_val <- wilcox.test(ref_values, test_values)$p.value
        
        star <- ""
        if (p_val < 0.001) {
          star <- "***"
        } else if (p_val < 0.01) {
          star <- "**"
        } else if (p_val < 0.05) {
          star <- "*"
        }
        
        if (star != "") {
          res <- rbind(
            res,
            data.frame(
              diagose1    = dg,
              sampletype2 = st,
              label      = star,
              p_val      = p_val,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }
  res
}

annotation_df <- get_significance(dat)

# Merge stats with SE for y-position of stars
annotation_df <- merge(
  annotation_df,
  SE_filtered,
  by = c("diagose1", "sampletype2"),
  all.x = TRUE
)

# Set star y-position: mean + SE + 0.05
annotation_df$y_pos <- annotation_df$value + annotation_df$se + 0.05

# Export statistical results
write.csv(
  annotation_df,
  file = file.path(out_dir, paste0(out_prefix, ".shannon.sampletype2.wilcox.csv")),
  row.names = FALSE
)

########################################
## Custom sampletype order & colors
########################################
sampletype_order <- c("Tongue coating",
                      "Gastric fluid",
                      "Gastric tissue",
                      "Stool")

sample_colors <- c(
  "Tongue coating" = "#5BA8A2",
  "Gastric fluid"  = "#F6B366",
  "Gastric tissue" = "#6B8FAD",
  "Stool"          = "#F2A099"
)

# Set factor levels
SE_filtered$sampletype2  <- factor(SE_filtered$sampletype2,
                                   levels = sampletype_order)
annotation_df$sampletype2 <- factor(annotation_df$sampletype2,
                                    levels = sampletype_order)

########################################
## Plot: each sampletype gets its own panel (free y-axis)
########################################

f2 <- ggplot(SE_filtered,
             aes(x = diagose1,
                 y = value,
                 group = sampletype2,
                 color = sampletype2)) +
  
  # Each panel draws its own line only
  geom_line(aes(y = value_line), 
            linetype = "dashed") +
  
  geom_point(size = 2) +
  
  geom_errorbar(aes(ymin = value - se,
                    ymax = value + se),
                width = 0.1) +
  
  # Add significance stars
  geom_text(
    data = annotation_df,
    aes(x = diagose1,
        y = y_pos,
        label = label,
        color = sampletype2),
    size = 4,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  # Custom colors
  scale_color_manual(values = sample_colors) +
  
  # Facet by sampletype, free y-axis per panel
  facet_wrap(~ sampletype2,
             ncol = 1,
             scales = "free_y") +
  
  theme_classic() +
  labs(
    title = NULL,
    x     = "",
    y     = "Shannon"
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.text.x      = element_text(angle = 90,
                                    vjust = 0.5,
                                    hjust = 1),
    axis.title.y     = element_text(size = 11, face = "bold")
  )

print(f2)