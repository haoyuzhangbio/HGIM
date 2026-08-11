############################################################
## Figure 2c proportion plot for Mac
##
## Input:
##   ~/Desktop/Figure2d.pattern.addprecancer_HSSG.csv
##
## Output:
##   ~/Desktop/figure2c.prop_plot.pdf
##   ~/Desktop/figure2c.prop_plot.png
##   ~/Desktop/figure2c.prop_plot_data.csv
############################################################

## =========================================================
## 1. Basic settings
## =========================================================

input_file <- "~/Desktop/Fig2_20260720/Figure2d.pattern.kraken.csv"
out_dir <- "~/Desktop"

if (!file.exists(input_file)) {
  stop(paste0("Cannot find input file: ", input_file))
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## =========================================================
## 2. Install and load R packages
## =========================================================

packages <- c("ggplot2", "dplyr", "tidyr", "scales")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

## =========================================================
## 3. Read input file
## =========================================================

pat <- read.csv(
  input_file,
  as.is = TRUE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

## =========================================================
## 4. Check required columns
## =========================================================

required_cols <- c(
  "group2",
  "source",
  "gene",
  "fc",
  "glmqval_more",
  "glmpval_more"
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
## 5. Set disease stages and sample types
## =========================================================

diagnosis_order <- c(
  "Atrophic gastritis",
  "Intestinal metaplasia",
  "Intraepithelial Neoplasia"
)

sample_order <- c(
  "Tongue coating",
  "Gastric fluid",
  "Gastric tissue",
  "Stool"
)

available_diagnosis <- diagnosis_order[diagnosis_order %in% unique(pat$group2)]
available_samples <- sample_order[sample_order %in% unique(pat$source)]

if (length(available_diagnosis) == 0) {
  stop("Input file does not contain Atrophic gastritis / Intestinal metaplasia / Intraepithelial Neoplasia in group2.")
}

if (!"GC vs nonGC" %in% unique(pat$group2)) {
  stop("Input file does not contain group2 == 'GC vs nonGC'.")
}

if (length(available_samples) == 0) {
  stop("Input file does not contain Tongue coating / Gastric fluid / Gastric tissue / Stool in source.")
}

## =========================================================
## 6. Generate ph_long
##
## Logic:
## For each precancerous stage x each sample source:
## 1. Take species significant in GC vs nonGC
## 2. Match those species in the precancerous stage's fc and p values
## 3. Determine whether GC direction and precancerous direction are consistent
## =========================================================

ph_long <- data.frame()

for (gg in available_diagnosis) {
  for (ss in available_samples) {
    
    ## Current precancerous stage data
    p_pre <- pat[
      pat$group2 == gg &
        pat$source == ss,
    ]
    
    p_pre <- p_pre[
      !duplicated(p_pre[, c("gene", "source", "group2")]),
    ]
    
    ## GC vs nonGC significant species
    pj <- pat[
      pat$group2 == "GC vs nonGC" &
        pat$glmqval_more < 0.05 &
        pat$source == ss,
    ]
    
    pj <- pj[
      !duplicated(pj[, c("gene", "source", "group2")]),
    ]
    
    ## Skip if no data for this combination
    if (nrow(pj) == 0 | nrow(p_pre) == 0) {
      next
    }
    
    ## Match the same species in the precancerous stage
    idx <- match(pj$gene, p_pre$gene)
    
    dd <- data.frame(
      species = pj$gene,
      GC = as.numeric(pj$fc),
      Precancer = as.numeric(p_pre[idx, "fc"]),
      Precancer_q = as.numeric(p_pre[idx, "glmqval_more"]),
      Precancer_p = as.numeric(p_pre[idx, "glmpval_more"]),
      diagnosis = gg,
      source = ss,
      stringsAsFactors = FALSE
    )
    
    ## Remove species not matched in the precancerous stage
    dd <- dd[!is.na(dd$Precancer), ]
    
    ph_long <- rbind(ph_long, dd)
  }
}

if (nrow(ph_long) == 0) {
  stop("No data generated for plotting. Check whether GC vs nonGC significant species can be matched in the precancerous stages.")
}

## =========================================================
## 7. Classification
##
## Enriched coherent:
##   GC fc > 0 and precancerous stage fc > 0
##
## Depleted coherent:
##   GC fc < 0 and precancerous stage fc < 0
##
## Sig / NonSig:
##   Based on precancerous stage glmpval_more < 0.05
##
## Other:
##   GC and precancerous stage directions disagree, or one fc = 0
## =========================================================

ph_long <- ph_long %>%
  mutate(
    category = case_when(
      GC > 0 & Precancer > 0 & Precancer_p < 0.05 ~
        "Enriched coherent & Sig",
      
      GC > 0 & Precancer > 0 & Precancer_p >= 0.05 ~
        "Enriched coherent & NonSig",
      
      GC < 0 & Precancer < 0 & Precancer_p < 0.05 ~
        "Depleted coherent & Sig",
      
      GC < 0 & Precancer < 0 & Precancer_p >= 0.05 ~
        "Depleted coherent & NonSig",
      
      TRUE ~ "Other"
    )
  )

## =========================================================
## 8. Calculate proportions
## =========================================================

prop_data <- ph_long %>%
  group_by(diagnosis, source, category) %>%
  summarise(
    count = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

## Fill missing categories to ensure complete stacking
category_order <- c(
  "Enriched coherent & Sig",
  "Enriched coherent & NonSig",
  "Depleted coherent & Sig",
  "Depleted coherent & NonSig",
  "Other"
)

prop_data <- prop_data %>%
  complete(
    diagnosis = available_diagnosis,
    source = available_samples,
    category = category_order,
    fill = list(count = 0, proportion = 0)
  )

prop_data$diagnosis <- factor(
  prop_data$diagnosis,
  levels = diagnosis_order
)

prop_data$source <- factor(
  prop_data$source,
  levels = sample_order
)

prop_data$category <- factor(
  prop_data$category,
  levels = category_order
)

print(prop_data)

## =========================================================
## 9. Save plot data
## =========================================================

write.csv(
  prop_data,
  file.path(out_dir, "figure2c.prop_plot_data.csv"),
  row.names = FALSE
)

## =========================================================
## 10. Plot
## =========================================================

category_colors <- c(
  "Enriched coherent & Sig" = "#E64B35",
  "Enriched coherent & NonSig" = "#F39B7F",
  "Depleted coherent & Sig" = "#4DBBD5",
  "Depleted coherent & NonSig" = "#3C5488",
  "Other" = "gray80"
)

prop_plot <- ggplot(
  prop_data,
  aes(x = source, y = proportion, fill = category)
) +
  geom_col(
    position = "stack",
    width = 0.7,
    color = NA
  ) +
  scale_fill_manual(
    values = category_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Sample source",
    y = "Proportion",
    fill = "Category",
    title = "Species composition by precancerous stage and source"
  ) +
  facet_wrap(
    ~ diagnosis,
    ncol = 3
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 11
    ),
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    axis.title = element_text(
      color = "black",
      size = 13
    ),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

print(prop_plot)

## =========================================================
## 11. Output PDF and PNG
## =========================================================

ggsave(
  filename = file.path(out_dir, "figure2c.prop_plot.pdf"),
  plot = prop_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "figure2c.prop_plot.png"),
  plot = prop_plot,
  width = 7,
  height = 5,
  dpi = 300
)

message("Done!")
message("Output files:")
message(file.path(out_dir, "figure2c.prop_plot.pdf"))
message(file.path(out_dir, "figure2c.prop_plot.png"))
message(file.path(out_dir, "figure2c.prop_plot_data.csv"))