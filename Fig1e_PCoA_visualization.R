############################################################
## Reproduce PCoA plot from existing PCoA results
############################################################

setwd("~/Desktop/Fig.1c_20260717")

library(readxl)
library(dplyr)
library(ggplot2)

############################
## 1. Read existing results
############################

cmds_dat <- read.csv(
  "bray_distance.pcoa.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pheno <- read_excel(
  "pcoa_pheno.xlsx",
  col_types = "text"
)

pheno <- as.data.frame(pheno)

var.explained <- readRDS("~/Desktop/....../bray_distance_output_0715/var.explained.rds")

############################
## 2. Check column names
############################

colnames(cmds_dat)
colnames(pheno)

# cmds_dat must contain:
# SampleID, PCoA1, PCoA2

# pheno must contain:
# SampleID, sampletype, group

if (!all(c("SampleID", "PCoA1", "PCoA2") %in% colnames(cmds_dat))) {
  stop("bray_distance.pcoa.csv must contain SampleID, PCoA1, PCoA2")
}

if (!all(c("SampleID", "sampletype", "group") %in% colnames(pheno))) {
  stop("pcoa_pheno.xlsx must contain SampleID, sampletype, group")
}

############################
## 3. Merge PCoA coordinates and phenotype info
############################

pcoa_pheno <- cmds_dat %>%
  left_join(pheno, by = "SampleID")

# Check for unmatched samples
cat("sampletype missing count:", sum(is.na(pcoa_pheno$sampletype)), "\n")
cat("group missing count:", sum(is.na(pcoa_pheno$group)), "\n")

table(pcoa_pheno$sampletype, useNA = "ifany")
table(pcoa_pheno$group, useNA = "ifany")

############################
## 4. Set colors and shapes
############################

pcoa_pheno$sampletype <- factor(
  pcoa_pheno$sampletype,
  levels = c(
    "Stool",
    "Gastric fluid",
    "Gastric tissue",
    "Tongue coating"
  )
)

sampletype_colors <- c(
  "Stool"          = "#F2A099",
  "Gastric fluid"  = "#F6B366",
  "Gastric tissue" = "#6B8FAD",
  "Tongue coating" = "#5BA8A2"
)

pcoa_pheno$group <- factor(
  pcoa_pheno$group,
  levels = c("GP", "GU", "HS", "SG", "AG", "IM", "IN", "GC")
)

group_shapes <- c(
  "GP" = 22,
  "GU" = 21,
  "HS" = 24,
  "SG" = 23,
  "AG" = 25,
  "IM" = 3,
  "IN" = 7,
  "GC" = 8
)

############################
## 5. Plot PCoA
############################

p <- ggplot(pcoa_pheno, aes(x = PCoA1, y = PCoA2)) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "grey80"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey80"
  ) +
  geom_point(
    aes(color = sampletype, shape = group),
    size = 2.8,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Sample type",
    values = sampletype_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    name = "Group",
    values = group_shapes,
    drop = FALSE
  ) +
  labs(
    x = paste0("PCoA1 (", round(var.explained[1] * 100, 2), "%)"),
    y = paste0("PCoA2 (", round(var.explained[2] * 100, 2), "%)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(4, "pt"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

print(p)

ggsave(
  "PCoA_reproduce.pdf",
  plot = p,
  width = 8,
  height = 7
)

############################################################
## PCoA plot with species star annotations (requires raw abundance matrix)
############################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggstar)

############################
## 2. Read relative abundance matrix
############################

abun_file <- "~/Desktop/20260715.total_bracken/20260715.total_bracken_profile.relative.txt"

df_raw <- read.table(
  abun_file,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

cat("Raw abundance table dimensions:\n")
print(dim(df_raw))

# First column is sample ID
sample_ids <- as.character(df_raw[[1]])
sample_ids <- trimws(sample_ids)

# Following columns are species names
species_names <- colnames(df_raw)[-1]

# rows = samples, columns = species
mat <- as.matrix(df_raw[, -1, drop = FALSE])
mode(mat) <- "numeric"

rownames(mat) <- sample_ids
colnames(mat) <- species_names

# Transpose to: rows = species, columns = samples
df <- as.data.frame(
  t(mat),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Transposed df dimensions:\n")
print(dim(df))

############################
## 3. Merge PCoA coordinates and phenotype
############################

pcoa_pheno <- cmds_dat %>%
  left_join(pheno, by = "SampleID")

cat("sampletype missing count:", sum(is.na(pcoa_pheno$sampletype)), "\n")
cat("group missing count:", sum(is.na(pcoa_pheno$group)), "\n")

############################
## 4. Set colors and shapes
############################

pcoa_pheno$sampletype <- factor(
  pcoa_pheno$sampletype,
  levels = c("Stool", "Gastric fluid", "Gastric tissue", "Tongue coating")
)

sampletype_colors <- c(
  "Stool"          = "#F2A099",
  "Gastric fluid"  = "#F6B366",
  "Gastric tissue" = "#6B8FAD",
  "Tongue coating" = "#5BA8A2"
)

pcoa_pheno$group <- factor(
  pcoa_pheno$group,
  levels = c("GP", "GU", "HS", "SG", "AG", "IM", "IN", "GC")
)

group_shapes <- c(
  "GP" = 22,
  "GU" = 21,
  "HS" = 24,
  "SG" = 23,
  "AG" = 25,
  "IM" = 3,
  "IN" = 7,
  "GC" = 8
)

############################
## 5. Compute species correlations with PCoA1/PCoA2
############################

# df: rows = species, columns = samples
# species_mat: rows = samples, columns = species
species_mat <- t(df)

# Keep only samples present in PCoA, aligned by pcoa_pheno order
common_samples <- intersect(pcoa_pheno$SampleID, rownames(species_mat))

cat("PCoA sample count:", nrow(pcoa_pheno), "\n")
cat("Abundance table sample count:", nrow(species_mat), "\n")
cat("Common sample count:", length(common_samples), "\n")

pcoa_pheno2 <- pcoa_pheno[
  pcoa_pheno$SampleID %in% common_samples,
  ,
  drop = FALSE
]

species_mat2 <- species_mat[
  match(pcoa_pheno2$SampleID, rownames(species_mat)),
  ,
  drop = FALSE
]

stopifnot(all(rownames(species_mat2) == pcoa_pheno2$SampleID))

# Spearman correlation
pc1_cor <- cor(
  pcoa_pheno2$PCoA1,
  species_mat2,
  method = "spearman",
  use = "pairwise.complete.obs"
)

pc2_cor <- cor(
  pcoa_pheno2$PCoA2,
  species_mat2,
  method = "spearman",
  use = "pairwise.complete.obs"
)

pc.msps <- data.frame(
  Taxon = colnames(species_mat2),
  PCoA1 = as.numeric(pc1_cor),
  PCoA2 = as.numeric(pc2_cor),
  stringsAsFactors = FALSE
)

# If taxon names contain "s__", keep species-level taxa only
if (any(grepl("s__", pc.msps$Taxon))) {
  pc.msps <- pc.msps[grepl("s__", pc.msps$Taxon), ]
}

cat("Species count for annotation selection:", nrow(pc.msps), "\n")

############################
## 6. Select top species correlated with PCoA1/PCoA2
############################

pc1_pos2 <- pc.msps %>%
  arrange(desc(PCoA1)) %>%
  slice(1:4)

pc1_neg2 <- pc.msps %>%
  arrange(PCoA1) %>%
  slice(1:4)

pc2_pos2 <- pc.msps %>%
  arrange(desc(PCoA2)) %>%
  slice(1:9)

pc2_neg2 <- pc.msps %>%
  arrange(PCoA2) %>%
  slice(1:4)

pcoa.anno <- bind_rows(
  pc1_pos2,
  pc1_neg2,
  pc2_pos2,
  pc2_neg2
) %>%
  distinct(Taxon, .keep_all = TRUE)

# Label name
pcoa.anno$label <- sub(".*s__", "s__", pcoa.anno$Taxon)

# If no s__ found, use Taxon directly
pcoa.anno$label[pcoa.anno$label == pcoa.anno$Taxon] <- pcoa.anno$Taxon[pcoa.anno$label == pcoa.anno$Taxon]

write.csv(
  pcoa.anno,
  "pcoa.anno.selected_species.csv",
  row.names = FALSE
)

############################
## 7. Scale species coordinates
############################

scale_x <- diff(range(pcoa_pheno2$PCoA1, na.rm = TRUE)) /
  diff(range(pcoa.anno$PCoA1, na.rm = TRUE)) * 0.75

scale_y <- diff(range(pcoa_pheno2$PCoA2, na.rm = TRUE)) /
  diff(range(pcoa.anno$PCoA2, na.rm = TRUE)) * 0.75

scale_factor <- min(scale_x, scale_y, na.rm = TRUE)

pcoa.anno$PCoA1_scaled <- pcoa.anno$PCoA1 * scale_factor

pcoa.anno$PCoA2_scaled <- pcoa.anno$PCoA2 * scale_factor

# Shift negative PCoA2 values upward proportionally
neg_y_shrink <- 0.75

pcoa.anno$PCoA2_scaled <- ifelse(
  pcoa.anno$PCoA2_scaled < 0,
  pcoa.anno$PCoA2_scaled * neg_y_shrink,
  pcoa.anno$PCoA2_scaled
)

neg_y_shrink <- 1.5

pcoa.anno$PCoA2_scaled <- ifelse(
  pcoa.anno$PCoA2_scaled > 0,
  pcoa.anno$PCoA2_scaled * neg_y_shrink,
  pcoa.anno$PCoA2_scaled
)

############################
## 8. Plot PCoA with species star annotations
############################

p <- ggplot(pcoa_pheno2, aes(x = PCoA1, y = PCoA2)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
  geom_point(
    aes(color = sampletype, shape = group),
    size = 2.8,
    alpha = 0.9
  ) +
  ggstar::geom_star(
    data = pcoa.anno,
    aes(x = PCoA1_scaled, y = PCoA2_scaled),
    inherit.aes = FALSE,
    starshape = 1,
    size = 4,
    color = "black",
    fill = "orange"
  ) +
  geom_text_repel(
    data = pcoa.anno,
    aes(x = PCoA1_scaled, y = PCoA2_scaled, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15,
    max.overlaps = 100,
    force = 10
  ) +
  scale_color_manual(
    name = "Sample type",
    values = sampletype_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    name = "Group",
    values = group_shapes,
    drop = FALSE
  ) +
  labs(
    x = paste0("PCoA1 (", round(var.explained[1] * 100, 2), "%)"),
    y = paste0("PCoA2 (", round(var.explained[2] * 100, 2), "%)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(4, "pt"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

print(p)

ggsave(
  "Fig.PCoA.with_species_stars.20260715.pdf",
  plot = p,
  width = 11,
  height = 10
)