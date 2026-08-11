############################################################
## Summarize group1 chisq.test.csv files
## All and H analyzed completely separately
############################################################
library(data.table)
library(dplyr)
library(stringr)
## =========================================================
## 1. Path configuration
## =========================================================
input_dir <- "/Volumes/thinkplus/大队列数据_Nature/Fig3c/output"
out_root <- "~/Desktop/group1_chisq_summary"
tax_file <- "/Volumes/thinkplus/大队列数据_Nature/Fig3c/complete_50.sgb_number_taxonomy.clean.fix_space.tax_sim.rename.tsv"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
qval_cutoff <- 0.05
############################################################
## Load SemiBin species taxonomy annotation
############################################################
tax <- fread(tax_file)
tax[, species_name := sub(".*\\|s__", "s__", GTDB_TAXONOMY)]
tax[
  is.na(species_name) |
    species_name == "" |
    species_name == GTDB_TAXONOMY,
  species_name := GTDB_TAXONOMY
]
tax[, species_name := gsub("_", " ", species_name)]
tax_map <- tax[, .(
  bin_id = RAW_SGB,
  species_name
)]
tax_map <- tax_map[!duplicated(bin_id)]
## =========================================================
## 2. Function to parse GC and Control proportions from string
## =========================================================
parse_gc_con_from_proportions <- function(x) {
  
  if (is.na(x) || x == "") {
    return(c(gc = NA_real_, con = NA_real_))
  }
  
  parts <- unlist(strsplit(x, ";", fixed = TRUE))
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  
  if (length(parts) == 0) {
    return(c(gc = NA_real_, con = NA_real_))
  }
  
  labels <- character(0)
  values <- numeric(0)
  
  for (p in parts) {
    kv <- unlist(strsplit(p, ":", fixed = TRUE))
    
    if (length(kv) < 2) next
    
    lab <- trimws(kv[1])
    val <- suppressWarnings(as.numeric(trimws(kv[2])))
    
    labels <- c(labels, lab)
    values <- c(values, val)
  }
  
  if (length(labels) == 0) {
    return(c(gc = NA_real_, con = NA_real_))
  }
  
  labels_clean <- tolower(labels)
  labels_clean <- gsub("[^a-z0-9]", "", labels_clean)
  
  gc_idx <- which(labels_clean %in% c("gc", "gastriccancer", "cancer"))
  con_idx <- which(labels_clean %in% c(
    "con", "control", "ctrl", "ngc", "nongc",
    "noncancer", "healthy", "healthycontrol"
  ))
  
  gc_val <- if (length(gc_idx) > 0) values[gc_idx[1]] else NA_real_
  con_val <- if (length(con_idx) > 0) values[con_idx[1]] else NA_real_
  
  c(gc = gc_val, con = con_val)
}
## =========================================================
## 3. Function to load single chi-square result file
## =========================================================
read_one_chisq_file <- function(file) {
  
  required_cols <- c(
    "gene",
    "test",
    "pval",
    "effect_size",
    "effect_direction",
    "proportions",
    "proportion_diff",
    "qval",
    "sampletype"
  )
  
  header <- tryCatch(
    names(fread(file, nrows = 0, showProgress = FALSE)),
    error = function(e) {
      warning("Failed to read header: ", file)
      return(NULL)
    }
  )
  
  if (is.null(header)) return(NULL)
  
  missing_cols <- setdiff(required_cols, header)
  
  if (length(missing_cols) > 0) {
    warning(
      "File skipped due to missing required columns: ",
      basename(file),
      " Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  dt <- tryCatch(
    fread(
      file,
      select = required_cols,
      data.table = TRUE,
      showProgress = FALSE
    ),
    error = function(e) {
      warning("File read failed, skipped: ", file)
      return(NULL)
    }
  )
  
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  dt <- dt[tolower(trimws(test)) == "group1"]
  
  if (nrow(dt) == 0) return(NULL)
  
  dt[, source_file := basename(file)]
  dt[, bin_id := sub("\\.gene_presence_absence.*$", "", basename(file))]
  
  dt
}
## =========================================================
## 4. Main function to analyze All or H dataset separately
## =========================================================
run_one_set <- function(set_name, file_pattern, input_dir, out_root, qval_cutoff = 0.1) {
  
  out_dir <- file.path(out_root, set_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n============================================================\n")
  cat("Start analysis:", set_name, "\n")
  cat("File matching pattern:", file_pattern, "\n")
  cat("Output directory:", out_dir, "\n")
  cat("============================================================\n")
  
  all_files <- list.files(
    path = input_dir,
    pattern = file_pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  
  cat("Total matched files found:", length(all_files), "\n")
  
  if (length(all_files) == 0) {
    warning(set_name, " No matched files found, skip this dataset.")
    return(NULL)
  }
  
  group1_list <- lapply(seq_along(all_files), function(i) {
    
    if (i %% 100 == 0) {
      cat(set_name, " Processed files count:", i, "/", length(all_files), "\n")
    }
    
    read_one_chisq_file(all_files[i])
  })
  
  group1_list <- group1_list[!sapply(group1_list, is.null)]
  
  if (length(group1_list) == 0) {
    warning(set_name, " No records with test == group1 in any matched file.")
    return(NULL)
  }
  
  group1_all <- rbindlist(
    group1_list,
    fill = TRUE,
    use.names = TRUE
  )
  
  ## Merge species annotation information
  group1_all <- merge(
    group1_all,
    tax_map,
    by = "bin_id",
    all.x = TRUE
  )
  
  ## Retain bin_id for entries without species annotation
  group1_all[
    is.na(species_name) | species_name == "",
    species_name := bin_id
  ]
  
  
  group1_all[, analysis_set := set_name]
  
  cat(set_name, " Total group1 records:", nrow(group1_all), "\n")
  cat(set_name, " Unique source files involved:", length(unique(group1_all$source_file)), "\n")
  cat(set_name, " Sampletype distribution:\n")
  print(table(group1_all$sampletype))
  cat("\n")
  
  ## Data type conversion
  group1_all[, gene := as.character(gene)]
  group1_all[, sampletype := as.character(sampletype)]
  group1_all[, qval := as.numeric(qval)]
  group1_all[, pval := as.numeric(pval)]
  group1_all[, effect_size := as.numeric(effect_size)]
  group1_all[, proportion_diff := as.numeric(proportion_diff)]
  group1_all[, effect_direction := as.numeric(effect_direction)]
  
  ## Filter statistically significant hits
  sig <- group1_all[
    !is.na(qval) &
      qval < qval_cutoff
  ]
  
  cat(set_name, " Significant records with qval <", qval_cutoff, ":", nrow(sig), "\n")
  
  fwrite(
    group1_all,
    file.path(out_dir, paste0(set_name, "_all_group1_records.csv"))
  )
  
  if (nrow(sig) == 0) {
    warning(set_name, " No significant group1 records pass qval cutoff.")
    return(list(
      all = group1_all,
      sig = sig
    ))
  }
  
  ## Parse group proportions string into separate columns
  prop_mat <- t(
    vapply(
      sig$proportions,
      parse_gc_con_from_proportions,
      numeric(2)
    )
  )
  
  sig[, proportion_gc := prop_mat[, "gc"]]
  sig[, proportion_con := prop_mat[, "con"]]
  sig[, gc_minus_con := proportion_gc - proportion_con]
  
  sig[, enrichment := fifelse(
    is.na(proportion_gc) | is.na(proportion_con),
    NA_character_,
    fifelse(
      proportion_gc > proportion_con,
      "GC_enriched",
      fifelse(
        proportion_gc < proportion_con,
        "Con_enriched",
        "No_difference"
      )
    )
  )]
  
  cat(set_name, " Enrichment direction summary:\n")
  print(table(sig$enrichment, useNA = "ifany"))
  cat("\n")
  
  parse_failed <- sig[
    is.na(proportion_gc) | is.na(proportion_con)
  ]
  
  if (nrow(parse_failed) > 0) {
    cat(set_name, " WARNING: Count of entries with unparseable proportions string:", nrow(parse_failed), "\n")
    print(head(unique(parse_failed$proportions), 10))
  }
  
  ## Export full significant result table
  fwrite(
    sig,
    file.path(out_dir, paste0(set_name, "_group1_significant_qval0.05_with_direction.csv"))
  )
  
  ## Summarize significant hits by sampletype
  summary_by_sampletype <- sig[
    ,
    .(
      significant_record_n = .N,
      significant_gene_n = uniqueN(gene),
      significant_file_n = uniqueN(source_file),
      
      GC_enriched_record_n = sum(enrichment == "GC_enriched", na.rm = TRUE),
      GC_enriched_gene_n = uniqueN(gene[enrichment == "GC_enriched"]),
      
      Con_enriched_record_n = sum(enrichment == "Con_enriched", na.rm = TRUE),
      Con_enriched_gene_n = uniqueN(gene[enrichment == "Con_enriched"]),
      
      No_difference_record_n = sum(enrichment == "No_difference", na.rm = TRUE),
      parse_failed_n = sum(is.na(enrichment))
    ),
    by = sampletype
  ][order(sampletype)]
  
  fwrite(
    summary_by_sampletype,
    file.path(out_dir, paste0(set_name, "_summary_by_sampletype_qval0.05.csv"))
  )
  
  ## Summarize hits grouped by sampletype + enrichment direction
  summary_by_sampletype_direction <- sig[
    !is.na(enrichment),
    .(
      record_n = .N,
      gene_n = uniqueN(gene),
      file_n = uniqueN(source_file),
      min_qval = min(qval, na.rm = TRUE),
      median_qval = median(qval, na.rm = TRUE),
      max_abs_gc_minus_con = max(abs(gc_minus_con), na.rm = TRUE)
    ),
    by = .(sampletype, enrichment)
  ][order(sampletype, enrichment)]
  
  fwrite(
    summary_by_sampletype_direction,
    file.path(out_dir, paste0(set_name, "_summary_by_sampletype_direction_qval0.05.csv"))
  )
  
  ## Extract top most significant enriched genes per sampletype
  sig_non_group <- sig[
    !is.na(enrichment) &
      enrichment %in% c("GC_enriched", "Con_enriched") &
      !grepl("^group", gene, ignore.case = TRUE)
  ]
  
  if (nrow(sig_non_group) == 0) {
    
    warning(set_name, " No remaining significant genes after filtering out genes starting with 'group'.")
    
    top_gene_by_sampletype <- data.table()
    top_gene_overall <- data.table()
    top20_gene_by_sampletype <- data.table()
    
  } else {
    
    sig_non_group[, abs_gc_minus_con := abs(gc_minus_con)]
    sig_non_group[, abs_effect_size := abs(effect_size)]
    
    setorder(
      sig_non_group,
      sampletype,
      enrichment,
      qval,
      -abs_gc_minus_con,
      -abs_effect_size
    )
    
    top_gene_by_sampletype <- sig_non_group[
      ,
      .SD[1],
      by = .(sampletype, enrichment)
    ][
      ,
      .(
        sampletype,
        enrichment,
        top_gene = gene,
        qval,
        pval,
        effect_size,
        effect_direction,
        proportion_gc,
        proportion_con,
        gc_minus_con,
        proportions,
        source_file,
        bin_id,
        species_name
      )
    ][order(sampletype, enrichment)]
    
    setorder(
      sig_non_group,
      enrichment,
      qval,
      -abs_gc_minus_con,
      -abs_effect_size
    )
    
    top_gene_overall <- sig_non_group[
      ,
      .SD[1],
      by = enrichment
    ][
      ,
      .(
        enrichment,
        top_gene = gene,
        sampletype,
        qval,
        pval,
        effect_size,
        effect_direction,
        proportion_gc,
        proportion_con,
        gc_minus_con,
        proportions,
        source_file,
        bin_id,
        species_name
      )
    ][order(enrichment)]
    
    top20_gene_by_sampletype <- sig_non_group[
      ,
      head(.SD, 20),
      by = .(sampletype, enrichment)
    ][
      ,
      .(
        sampletype,
        enrichment,
        gene,
        qval,
        pval,
        effect_size,
        proportion_gc,
        proportion_con,
        gc_minus_con,
        proportions,
        source_file,
        bin_id,
        species_name
      )
    ]
  }
  
  fwrite(
    top_gene_by_sampletype,
    file.path(out_dir, paste0(set_name, "_top_non_group_gene_by_sampletype_enrichment.csv"))
  )
  
  fwrite(
    top_gene_overall,
    file.path(out_dir, paste0(set_name, "_top_non_group_gene_overall_enrichment.csv"))
  )
  
  fwrite(
    top20_gene_by_sampletype,
    file.path(out_dir, paste0(set_name, "_top20_non_group_genes_by_sampletype_enrichment.csv"))
  )
  
  ## Summarize significant hits per bin file + sampletype
  summary_by_file_sampletype <- sig[
    ,
    .(
      significant_record_n = .N,
      significant_gene_n = uniqueN(gene),
      
      GC_enriched_record_n = sum(enrichment == "GC_enriched", na.rm = TRUE),
      GC_enriched_gene_n = uniqueN(gene[enrichment == "GC_enriched"]),
      
      Con_enriched_record_n = sum(enrichment == "Con_enriched", na.rm = TRUE),
      Con_enriched_gene_n = uniqueN(gene[enrichment == "Con_enriched"]),
      
      min_qval = min(qval, na.rm = TRUE)
    ),
    by = .(
      source_file,
      bin_id,
      species_name,
      sampletype
    )
  ][order(sampletype, min_qval)]
  
  fwrite(
    summary_by_file_sampletype,
    file.path(out_dir, paste0(set_name, "_summary_by_file_sampletype_qval0.05.csv"))
  )
  
  cat("\n", set_name, " Analysis finished.\n")
  cat("Output directory:", out_dir, "\n")
  
  invisible(list(
    all = group1_all,
    sig = sig,
    summary_by_sampletype = summary_by_sampletype,
    summary_by_sampletype_direction = summary_by_sampletype_direction,
    top_gene_by_sampletype = top_gene_by_sampletype,
    top_gene_overall = top_gene_overall,
    summary_by_file_sampletype = summary_by_file_sampletype,
    top20_gene_by_sampletype = top20_gene_by_sampletype
  ))
}
## =========================================================
## 5. Execute pipeline for All and H datasets separately
## =========================================================
res_All <- run_one_set(
  set_name = "All",
  file_pattern = "\\.All\\.chisq\\.test\\.csv$",
  input_dir = input_dir,
  out_root = out_root,
  qval_cutoff = qval_cutoff
)
res_H <- run_one_set(
  set_name = "H",
  file_pattern = "\\.H\\.chisq\\.test\\.csv$",
  input_dir = input_dir,
  out_root = out_root,
  qval_cutoff = qval_cutoff
)
cat("\nFull pipeline completed.\n")
cat("All dataset output directory:", file.path(out_root, "All"), "\n")
cat("H dataset output directory:", file.path(out_root, "H"), "\n")
############################################################
## Top species ranked by count of significantly differentially carried genes
##
## Filter criteria:
## 1. qval < 0.05
## 2. abs(GC proportion - Control proportion) > 0.1
## 3. Exclude genes whose name starts with "group"
## 4. Truncate gene name to content before first "~"
## 5. Map SemiBin bin IDs to full species taxonomy names
############################################################
## =========================================================
## 0. Load required packages
## =========================================================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(egg)
## =========================================================
## 1. Parameter configuration
## =========================================================
base_dir <- path.expand("~/Desktop/group1_chisq_summary")
set_name <- "H"      ## Options: "All" or "H"
ftype <- "Gjuice"       ## Options: "fecal", "tongue", "Gjuice", etc.
qval_cutoff_plot <- 0.05
diff_cutoff_plot <- 0
top_n_species <- 10
min_gene_count <- 0
bin_width <- 0.02
tax_file <- "/Volumes/thinkplus/大队列数据_Nature/Fig3c/complete_50.sgb_number_taxonomy.clean.fix_space.tax_sim.rename.tsv"
sig_file <- file.path(
  base_dir,
  set_name,
  paste0(set_name, "_group1_significant_qval0.1_with_direction.csv")
)
out_pdf <- file.path(
  base_dir,
  set_name,
  paste0(
    set_name,
    ".top20_species_gene_direction.",
    ftype,
    ".qval",
    qval_cutoff_plot,
    ".absdiff",
    diff_cutoff_plot,
    ".pdf"
  )
)
cat("Reading significant result file:", sig_file, "\n")
cat("Reading species taxonomy annotation file:", tax_file, "\n")
cat("Output PDF path:", out_pdf, "\n")
## =========================================================
## 2. Load SemiBin species taxonomy table
## =========================================================
tax <- fread(tax_file)
required_tax_cols <- c("RAW_SGB", "GTDB_TAXONOMY")
missing_tax_cols <- setdiff(required_tax_cols, colnames(tax))
if(length(missing_tax_cols) > 0){
  stop(
    "Taxonomy table missing required columns: ",
    paste(missing_tax_cols, collapse = ", ")
  )
}
## Extract species-level taxonomic label
tax[, species_name := sub(".*\\|s__", "s__", GTDB_TAXONOMY)]
## Retain full taxonomy string if no species prefix s__ exists
tax[
  is.na(species_name) |
    species_name == "" |
    species_name == GTDB_TAXONOMY,
  species_name := GTDB_TAXONOMY
]
## Optional: Remove leading s__ from species labels (uncomment line below)
## tax[, species_name := sub("^s__", "", species_name)]
## Replace underscores with spaces for better plot readability
tax[, species_name := gsub("_", " ", species_name)]
## Create mapping table from bin ID to species name
tax_map <- tax[, .(
  bin_id = RAW_SGB,
  species_name
)]
tax_map <- tax_map[!duplicated(bin_id)]
## =========================================================
## 3. Load filtered significant statistical results
## =========================================================
sig <- fread(sig_file)
required_sig_cols <- c(
  "gene",
  "bin_id",
  "sampletype",
  "qval",
  "enrichment",
  "proportion_gc",
  "proportion_con",
  "gc_minus_con"
)
missing_sig_cols <- setdiff(required_sig_cols, colnames(sig))
if(length(missing_sig_cols) > 0){
  stop(
    "Significant result table missing required columns: ",
    paste(missing_sig_cols, collapse = ", ")
  )
}
sig[, gene := as.character(gene)]
sig[, bin_id := as.character(bin_id)]
sig[, sampletype := as.character(sampletype)]
sig[, enrichment := as.character(enrichment)]
sig[, qval := as.numeric(qval)]
sig[, proportion_gc := as.numeric(proportion_gc)]
sig[, proportion_con := as.numeric(proportion_con)]
sig[, gc_minus_con := as.numeric(gc_minus_con)]
## Shorten gene names: remove all text after first tilde symbol
sig[, gene_short := sub("~.*$", "", gene)]
## =========================================================
## 4. Filter differential genes for target sample type
## =========================================================
sig_use <- sig[
  sampletype == ftype &
    enrichment %in% c("GC_enriched", "Con_enriched") &
   # !grepl("^group", gene, ignore.case = TRUE) &
    !is.na(qval) &
    !is.na(gc_minus_con) &
    qval < qval_cutoff_plot &
    abs(gc_minus_con) > diff_cutoff_plot
]
cat("Total filtered differential gene records:", nrow(sig_use), "\n")
if(nrow(sig_use) == 0){
  stop("No differential genes retained after filtering; check sampletype, qval_cutoff_plot or diff_cutoff_plot parameters.")
}
## =========================================================
## 5. Merge species annotation mapping
## =========================================================
sig_use <- merge(
  sig_use,
  tax_map,
  by = "bin_id",
  all.x = TRUE
)
## Fallback to bin ID for entries without matched species label
sig_use[
  is.na(species_name) | species_name == "",
  species_name := bin_id
]
cat("Unique species represented in filtered data:", uniqueN(sig_use$species_name), "\n")
## =========================================================
## 6. Count significant differential genes per species
## =========================================================
top_species <- sig_use %>%
  group_by(species_name) %>%
  summarise(
    total_count = n(),
    GC_enriched_n = sum(enrichment == "GC_enriched", na.rm = TRUE),
    Con_enriched_n = sum(enrichment == "Con_enriched", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_count)) %>%
  filter(total_count > min_gene_count) %>%
  slice_head(n = top_n_species)
cat("Top species retained for plotting:", nrow(top_species), "\n")
if(nrow(top_species) == 0){
  stop("No species meet the minimum gene count threshold; lower min_gene_count value.")
}
print(top_species)
## =========================================================
## 7. Left panel: Stacked bar chart of gene counts
## =========================================================
enrich_cols <- c(
  "GC_enriched" = "#e76254",
  "Con_enriched" = "#376795"
)
top_species_long <- top_species %>%
  pivot_longer(
    cols = c(GC_enriched_n, Con_enriched_n),
    names_to = "enrichment",
    values_to = "gene_count"
  ) %>%
  mutate(
    enrichment = recode(
      enrichment,
      "GC_enriched_n" = "GC_enriched",
      "Con_enriched_n" = "Con_enriched"
    ),
    species_name = factor(species_name, levels = rev(top_species$species_name))
  )
p1 <- ggplot(
  top_species_long,
  aes(x = species_name, y = gene_count, fill = enrichment)
) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = enrich_cols,
    labels = c(
      "Con_enriched" = "Con enriched",
      "GC_enriched" = "GC enriched"
    ),
    name = "Direction"
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Significantly differentially carried genes",
    title = paste0(set_name, " | ", ftype)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7, face = "italic"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
## =========================================================
## 8. Right panel dataset preparation: effect size distribution
## =========================================================
plot_data <- sig_use %>%
  filter(species_name %in% top_species$species_name) %>%
  mutate(
    species_name = factor(species_name, levels = rev(top_species$species_name)),
    direction = ifelse(gc_minus_con > 0, "GC_enriched", "Con_enriched"),
    gene_short = sub("~.*$", "", gene)
  )
## Extract top most significant gene per species & enrichment direction
top_genes_per_species <- plot_data %>%
  filter(!grepl("^group", gene, ignore.case = TRUE)) %>%
  group_by(species_name, direction) %>%
  arrange(qval, desc(abs(gc_minus_con))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    species_name,
    direction,
    gene_short,
    gc_minus_con,
    qval
  )
cat("Top representative gene labels per species & direction:\n")
print(top_genes_per_species)
## =========================================================
## 9. Bin effect size values and calculate within-species proportions
## =========================================================
plot_data_prop <- plot_data %>%
  group_by(species_name) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  mutate(
    bin = round(gc_minus_con / bin_width) * bin_width
  ) %>%
  group_by(species_name, bin, direction, total_count) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = count / total_count
  )
top_genes_label <- top_genes_per_species %>%
  mutate(
    bin = round(gc_minus_con / bin_width) * bin_width,
    is_top_gene = TRUE
  ) %>%
  select(
    species_name,
    bin,
    gene_short,
    is_top_gene
  )
species_levels <- levels(plot_data$species_name)
species_pos_df <- data.frame(
  species_name = factor(species_levels, levels = species_levels),
  y_pos = seq_along(species_levels)
)
plot_data_prop <- plot_data_prop %>%
  mutate(
    species_name = factor(species_name, levels = species_levels)
  ) %>%
  left_join(
    species_pos_df,
    by = "species_name"
  ) %>%
  left_join(
    top_genes_label,
    by = c("species_name", "bin")
  ) %>%
  mutate(
    is_top_gene = ifelse(is.na(is_top_gene), FALSE, TRUE),
    gene_label = ifelse(is_top_gene, gene_short, "")
  )
## =========================================================
## 10. Right panel: Effect size distribution plot
## =========================================================
p2 <- ggplot(
  plot_data_prop,
  aes(
    x = bin,
    y = y_pos,
    color = direction
  )
) +
  geom_segment(
    aes(
      xend = bin,
      y = y_pos - 0.3 * proportion,
      yend = y_pos + 0.3 * proportion,
      size = proportion
    ),
    lineend = "round"
  ) +
  scale_color_manual(
    values = enrich_cols,
    labels = c(
      "Con_enriched" = "Con enriched",
      "GC_enriched" = "GC enriched"
    ),
    name = "Direction"
  ) +
  scale_size_continuous(
    range = c(0, 3),
    guide = "none"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  ggrepel::geom_text_repel(
    data = plot_data_prop %>% filter(gene_label != ""),
    aes(
      x = bin,
      y = y_pos,
      label = gene_label
    ),
    color = "black",
    size = 2.3,
    min.segment.length = 0,
    segment.size = 0.2,
    box.padding = 0.2,
    max.overlaps = Inf,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    breaks = species_pos_df$y_pos,
    labels = species_pos_df$species_name,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  labs(
    x = "GC minus Con proportion",
    y = NULL,
    title = "Effect size distribution"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
## =========================================================
## 11. Combine two panels and export PDF
## =========================================================
combined <- egg::ggarrange(
  p1,
  p2,
  nrow = 1,
  ncol = 2,
  widths = c(2, 4)
)
ggsave(
  out_pdf,
  combined,
  width = 5,
  height = 4
)
cat("\nPlot export finished:\n")
cat(out_pdf, "\n")
