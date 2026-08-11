############################
## 1. Load required R packages
############################
library(data.table)
library(Seurat)
library(vegan)
library(cluster)
library(mclust)
library(readxl)

setwd("~/Desktop")

############################
## 2. Read data
############################
# Species abundance table
df <- fread("~/Desktop/20260715.total_bracken/MAG_COUNT_ge5_filtered_expression/Gastric_tissue_selected_species_relative_abundance.MAG_COUNT_ge5.filtered.csv",
            sep = ",", header = TRUE)

df <- as.data.frame(df)

# Set first column as row names (species/feature names)
rownames(df) <- df[, 1]
df <- df[, -1, drop = FALSE]

grx2 <- read_excel(
  "~/Desktop/20260715.total_bracken/ALL11.17.xlsx",
  col_types = "text"
)
grx2 <- as.data.frame(grx2)
gtissue_ids <- grx2$sample_ID[grx2$sampletype == "Gastric tissue"]

head(gtissue_ids, 20)

############################
## 3. Sample alignment
############################
nm <- intersect(colnames(df), grx2$sample_ID)

if (length(nm) == 0) {
  stop("Column names of df and grx2$sample_ID have no overlap. Please check sample name consistency.")
}

# Align to the same order
df <- df[, match(nm, colnames(df)), drop = FALSE]
grx2 <- grx2[match(nm, grx2$sample_ID), , drop = FALSE]

# Verify full alignment
stopifnot(all(colnames(df) == grx2$sample_ID))

############################
## 4. Split data by sampletype
############################
t_dat <- df[, which(grx2$sampletype == "Gastric tissue"), drop = FALSE]

if (ncol(t_dat) == 0) {
  stop("No samples found with sampletype == 'Gastric tissue'.")
}

############################
## 5. Remove all-zero samples and features
############################
# Remove columns (samples) with sum = 0
t_dat <- t_dat[, apply(t_dat, 2, sum, na.rm = TRUE) > 0, drop = FALSE]

# Remove rows (features) with sum = 0
t_dat <- t_dat[apply(t_dat, 1, sum, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(t_dat) == 0 || ncol(t_dat) == 0) {
  stop("t_dat is empty after filtering. Please check the data.")
}

abun.frac <- t_dat

############################
## 6. Compute Bray-Curtis distance + PCoA
############################
# vegan::vegdist input requires: rows = samples, columns = features
# abun.frac is currently rows = features, columns = samples, so transpose
bray_dist <- vegdist(t(abun.frac), method = "bray")

# PCoA
# k cannot exceed (number of samples - 1)
k_dim <- min(10, ncol(abun.frac) - 1)

if (k_dim < 2) {
  stop("Too few Gastric tissue samples for a valid PCoA.")
}

pcoa_result <- cmdscale(bray_dist, k = k_dim, eig = TRUE)

# Convert to embedding matrix
pcoa_embeddings <- as.matrix(pcoa_result$points)
colnames(pcoa_embeddings) <- paste0("PCoA_", seq_len(ncol(pcoa_embeddings)))

############################
## 7. Fixed k-means clustering
############################
set.seed(123)
fixed_k <- 4

# Use the first few PCoA dimensions for clustering
use_dims <- 1:min(2, ncol(pcoa_embeddings))

km_res <- kmeans(
  pcoa_embeddings[, use_dims, drop = FALSE],
  centers = fixed_k,
  nstart = 100
)

cluster_membership <- factor(km_res$cluster)

names(cluster_membership) <- rownames(pcoa_embeddings)

print(table(cluster_membership))
print(head(cluster_membership))

############################
## 8. Save clustering results
############################
cluster_df <- data.frame(
  SampleID = names(cluster_membership),
  Cluster = as.character(cluster_membership),
  stringsAsFactors = FALSE
)

write.csv(cluster_df,
          file = "Gtissue_cluster_membership.csv",
          row.names = FALSE)

library(pheatmap)

# Extract cluster info
cluster_df <- data.frame(
  SampleID = names(cluster_membership),
  Cluster = as.character(cluster_membership),
  stringsAsFactors = FALSE
)

# Ensure sample order matches the abundance matrix
cluster_df <- cluster_df[match(colnames(abun.frac), cluster_df$SampleID), , drop = FALSE]

# Verify alignment
stopifnot(all(colnames(abun.frac) == cluster_df$SampleID))

# Order samples by cluster
ord <- order(cluster_df$Cluster)
mat <- abun.frac[, ord, drop = FALSE]
anno_col <- data.frame(Cluster = cluster_df$Cluster[ord])
rownames(anno_col) <- cluster_df$SampleID[ord]

# Log-transform abundance to avoid extreme values in display
mat_log <- log10(mat + 1e-10)

# Optionally select the top 50 most variable features for a clearer heatmap
feature_var <- apply(mat_log, 1, var, na.rm = TRUE)
top_n <- min(93, nrow(mat_log))
top_features <- names(sort(feature_var, decreasing = TRUE))[1:top_n]
mat_plot <- mat_log[top_features, , drop = FALSE]

# Draw heatmap
pheatmap::pheatmap(
  mat_plot,
  annotation_col = anno_col,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Heatmap of Top Variable Features by Cluster"
)

bk <- seq(
  -3,
  3,
  length.out = 101
)

pheatmap::pheatmap(
  mat_plot,
  annotation_col = anno_col,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = bk,
  main = "Heatmap of Top Variable Features by Cluster"
)

# Extract cluster result
cluster_df <- data.frame(
  ID = trimws(as.character(names(cluster_membership))),
  Cluster = as.character(cluster_membership),
  stringsAsFactors = FALSE
)

# Clean grx2
grx2$sample_ID <- trimws(as.character(grx2$sample_ID))
grx2$group <- trimws(as.character(grx2$group))

# Merge: cluster_df$ID corresponds to grx2$sample_ID
plot_df <- merge(
  cluster_df,
  grx2[, c("sample_ID", "group")],
  by.x = "ID",
  by.y = "sample_ID"
)

head(plot_df)
table(plot_df$Cluster, plot_df$group)

# View first few rows
head(plot_df)
tab <- table(plot_df$group, plot_df$Cluster)
print(tab)
tab_df <- as.data.frame.matrix(tab)
tab_df <- cbind(Group = rownames(tab_df), tab_df)

write.csv(tab_df, "group_cluster_count.csv", row.names = FALSE)

library(ggplot2)
library(reshape2)

df_plot <- as.data.frame.matrix(tab)
df_plot$Group <- rownames(df_plot)

df_long <- melt(df_plot, id.vars = "Group", variable.name = "Cluster", value.name = "Count")

df_long <- within(df_long, {
  Proportion <- ave(Count, Group, FUN = function(x) x / sum(x) * 100)
})

cluster_labels <- c(
  "1" = "Cluster 1",
  "2" = "Cluster 2",
  "3" = "Cluster 3",
  "4" = "Cluster 4"
)

# Key: levels determine stacking order; the first level is at the bottom
df_long$ClusterName <- factor(
  cluster_labels[df_long$Cluster],
  levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
)

df_long$Group <- factor(df_long$Group, levels = c("GP","GU","HS","SG","AG","IM","IN","GC"))

my_colors <- c(
  "Cluster 1" = "#6B8FAD",
  "Cluster 2" = "#F6B366",
  "Cluster 3" = "#5BA8A2",
  "Cluster 4" = "#F2A099"
)

p <- ggplot(df_long, aes(x = Group, y = Proportion, fill = ClusterName)) +
  geom_bar(stat = "identity", width = 0.95, color = "white", size = 0.2) +
  scale_fill_manual(values = my_colors, name = "Gastric tissue communities type") +
  labs(x = NULL, y = "Proportion (%)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18)
  )

print(p)