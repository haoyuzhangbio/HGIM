############################################################
## Plot intra-person similarity from intra_pd1.csv
##
## intra_pd1.csv:
##   Each row = Bray distance between two sampletypes from the same uid
##   similarity = 1 - bray_distance
##
## Output:
##   Figure3.intra_similarity.boxplot.pdf
##   Figure3.intra_similarity.boxplot.png
##   GLM_vs_HealthySuperficial.csv
############################################################

## ================== 0. Basic settings ==================

setwd("~/Desktop")

input_file <- "~/Desktop/intra_pd1.csv"
out_dir <- "~/Desktop"

if (!file.exists(input_file)) {
  stop(paste0("Cannot find file: ", input_file))
}

## ================== 1. Load packages ==================

packages <- c("ggplot2", "dplyr", "ggsci")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(ggsci)

## ================== 2. Read intra_pd1.csv ==================

intra_pd1 <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(colnames(intra_pd1))
head(intra_pd1)

## ================== 3. Check required columns ==================

required_cols <- c(
  "type_combination",
  "bray_distance",
  "diagose"
)

missing_cols <- setdiff(required_cols, colnames(intra_pd1))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "intra_pd1.csv is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

## ================== 4. Prepare similarity ==================

## Here we use raw similarity: similarity = 1 - bray_distance
## No longer using log(1 - bray_distance + 0.01)

intra_pd1$bray_distance <- as.numeric(intra_pd1$bray_distance)

if (!"similarity" %in% colnames(intra_pd1)) {
  intra_pd1$similarity <- 1 - intra_pd1$bray_distance
} else {
  intra_pd1$similarity <- as.numeric(intra_pd1$similarity)
}

## If you are unsure whether the existing similarity is 1 - bray_distance, force overwrite:
intra_pd1$similarity <- 1 - intra_pd1$bray_distance

cat("similarity range:\n")
print(range(intra_pd1$similarity, na.rm = TRUE))

## ================== 5. Unify covariate column names for GLM ==================

## Your meta may have Age/Gender/Smoking/Drinking
## Here we unify to age/gender/smoke/drink

if (!"age" %in% colnames(intra_pd1) && "Age" %in% colnames(intra_pd1)) {
  intra_pd1$age <- intra_pd1$Age
}

if (!"gender" %in% colnames(intra_pd1) && "Gender" %in% colnames(intra_pd1)) {
  intra_pd1$gender <- intra_pd1$Gender
}

if (!"smoke" %in% colnames(intra_pd1) && "Smoking" %in% colnames(intra_pd1)) {
  intra_pd1$smoke <- intra_pd1$Smoking
}

if (!"drink" %in% colnames(intra_pd1) && "Drinking" %in% colnames(intra_pd1)) {
  intra_pd1$drink <- intra_pd1$Drinking
}

if ("BMI" %in% colnames(intra_pd1)) {
  intra_pd1$BMI <- as.numeric(intra_pd1$BMI)
}

if ("age" %in% colnames(intra_pd1)) {
  intra_pd1$age <- as.numeric(intra_pd1$age)
}

## ================== 6. Diagnosis group order ==================

diagose_levels <- c(
  "Gastric polyps",
  "Gastric ulcer",
  "Healthy stomach",
  "Superficial gastritis",
  "Atrophic gastritis",
  "Intestinal metaplasia",
  "Intraepithelial Neoplasia",
  "Gastric cancer"
)

intra_pd1$diagose <- factor(
  intra_pd1$diagose,
  levels = diagose_levels
)

cat("diagose distribution:\n")
print(table(intra_pd1$diagose, useNA = "always"))

cat("type_combination distribution:\n")
print(table(intra_pd1$type_combination, useNA = "always"))

## ================== 7. Create merged control group variable ==================

## Healthy stomach + Superficial gastritis merged as Healthy_Superficial
## Only used for statistics, does not change x-axis labels on the plot

intra_pd1$diagose2 <- as.character(intra_pd1$diagose)

intra_pd1$diagose2[
  intra_pd1$diagose %in% c("Healthy stomach", "Superficial gastritis")
] <- "Healthy_Superficial"

intra_pd1$diagose2 <- factor(
  intra_pd1$diagose2,
  levels = c(
    "Healthy_Superficial",
    "Gastric polyps",
    "Gastric ulcer",
    "Atrophic gastritis",
    "Intestinal metaplasia",
    "Intraepithelial Neoplasia",
    "Gastric cancer"
  )
)

## ================== 8. Keep only the 6 desired type_combinations ==================

desired_order <- c(
  "Gjuice_tongue",
  "Gtissue_tongue",
  "Gjuice_Gtissue",
  "fecal_tongue",
  "fecal_Gjuice",
  "fecal_Gtissue"
)

plot_data <- intra_pd1 %>%
  filter(type_combination %in% desired_order) %>%
  mutate(
    type_combination = factor(type_combination, levels = desired_order)
  )

if (nrow(plot_data) == 0) {
  stop("No type_combination matched desired_order. Please check the naming of intra_pd1$type_combination.")
}

## ================== 9. GLM: compute p/q values ==================

## For each type_combination:
## Other disease groups vs Healthy stomach + Superficial gastritis
##
## The target variable uses similarity, not bray_distance,
## so that the statistical direction matches the y-axis of the plot.

results_list_glm <- list()

unique_type_combinations <- unique(plot_data$type_combination)

for (type_combo in unique_type_combinations) {
  
  sub <- plot_data %>%
    filter(type_combination == type_combo)
  
  if (!"Healthy_Superficial" %in% sub$diagose2) next
  
  groups_to_compare <- setdiff(levels(intra_pd1$diagose2), "Healthy_Superficial")
  groups_to_compare <- intersect(groups_to_compare, as.character(sub$diagose2))
  
  if (length(groups_to_compare) == 0) next
  
  comparisons <- character()
  diagose_vec <- character()
  p_values <- numeric()
  
  for (g in groups_to_compare) {
    
    tmp <- sub %>%
      filter(diagose2 %in% c("Healthy_Superficial", g))
    
    ## Base model data
    d <- data.frame(
      target = as.numeric(scale(tmp$similarity)),
      group2 = factor(
        tmp$diagose2,
        levels = c("Healthy_Superficial", g)
      ),
      stringsAsFactors = FALSE
    )
    
    ## Add covariates if present
    if ("age" %in% colnames(tmp)) {
      d$age <- as.numeric(tmp$age)
    }
    
    if ("gender" %in% colnames(tmp)) {
      d$gender <- as.factor(tmp$gender)
    }
    
    if ("BMI" %in% colnames(tmp)) {
      d$BMI <- as.numeric(tmp$BMI)
    }
    
    if ("smoke" %in% colnames(tmp)) {
      d$smoke <- as.factor(tmp$smoke)
    }
    
    if ("drink" %in% colnames(tmp)) {
      d$drink <- as.factor(tmp$drink)
    }
    
    d <- na.omit(d)
    
    if (nrow(d) < 5) next
    if (length(unique(d$group2)) < 2) next
    
    ## Automatically build formula based on existing covariates
    covariates <- setdiff(colnames(d), c("target", "group2"))
    
    if (length(covariates) > 0) {
      formula_glm <- as.formula(
        paste("group2 ~ target +", paste(covariates, collapse = " + "))
      )
    } else {
      formula_glm <- group2 ~ target
    }
    
    model <- tryCatch(
      glm(
        formula_glm,
        data = d,
        family = binomial()
      ),
      error = function(e) NULL
    )
    
    if (is.null(model)) next
    
    s <- summary(model)
    
    if (!"target" %in% rownames(coef(s))) next
    
    p <- coef(s)["target", "Pr(>|z|)"]
    
    comparisons <- c(
      comparisons,
      paste("Healthy_Superficial", g, sep = " vs ")
    )
    
    diagose_vec <- c(diagose_vec, as.character(g))
    p_values <- c(p_values, p)
  }
  
  if (length(p_values) == 0) next
  
  q_values <- p.adjust(p_values, method = "BH")
  
  df_glm <- data.frame(
    type_combination = as.character(type_combo),
    diagose = diagose_vec,
    comparison = comparisons,
    p_value = p_values,
    q_value = q_values,
    stringsAsFactors = FALSE
  )
  
  results_list_glm[[as.character(type_combo)]] <- df_glm
}

if (length(results_list_glm) > 0) {
  final_results_glm <- do.call(rbind, results_list_glm)
  rownames(final_results_glm) <- NULL
} else {
  final_results_glm <- data.frame(
    type_combination = character(),
    diagose = character(),
    comparison = character(),
    p_value = numeric(),
    q_value = numeric(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  final_results_glm,
  file = file.path(out_dir, "GLM_vs_HealthySuperficial.csv"),
  row.names = FALSE
)

cat("GLM results:\n")
print(final_results_glm)

## ================== 10. Construct significance stars ==================

sig_results <- final_results_glm %>%
  filter(q_value < 0.05)

anno_df <- data.frame()

if (nrow(sig_results) > 0) {
  
  max_y <- plot_data %>%
    group_by(type_combination, diagose) %>%
    summarise(
      similarity = max(similarity, na.rm = TRUE),
      .groups = "drop"
    )
  
  anno_df <- sig_results %>%
    left_join(
      max_y,
      by = c("type_combination", "diagose")
    ) %>%
    mutate(
      y_position = similarity + 0.04,
      stars = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01  ~ "**",
        q_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      diagose = factor(diagose, levels = diagose_levels),
      type_combination = factor(type_combination, levels = desired_order)
    )
}

## ================== 11. Set colors ==================

box_colors <- c(
  "Gastric polyps"             = "#F49600",
  "Gastric ulcer"              = "#9E9AC8",
  "Healthy stomach"            = "#376795",
  "Superficial gastritis"      = "#4292C6",
  "Atrophic gastritis"         = "#F4BCA4",
  "Intestinal metaplasia"      = "#F09E89",
  "Intraepithelial Neoplasia"  = "#EB806F",
  "Gastric cancer"             = "#E76254"
)

## ================== 12. Draw boxplot ==================

p <- ggplot(
  plot_data,
  aes(
    x = diagose,
    y = similarity
  )
) +
  geom_boxplot(
    fill = NA,
    aes(color = diagose),
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = box_colors,
    drop = FALSE
  ) +
  facet_wrap(
    ~ type_combination,
    ncol = 6,
    scales = "fixed"
  ) +
  ylab("Similarity index (1 - Bray-Curtis distance)") +
  xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = unit(3, "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(
      size = 11,
      color = "black"
    ),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

if (nrow(anno_df) > 0) {
  p <- p +
    geom_text(
      data = anno_df,
      aes(
        x = diagose,
        y = y_position,
        label = stars
      ),
      inherit.aes = FALSE,
      size = 5,
      color = "black"
    )
}

print(p)

## ================== 13. Save figure ==================

ggsave(
  filename = file.path(out_dir, "Figure3.intra_similarity.boxplot.pdf"),
  plot = p,
  width = 13,
  height = 5
)

message("Done!")
message("Output files:")
message(file.path(out_dir, "Figure3.intra_similarity.boxplot.pdf"))
message(file.path(out_dir, "GLM_vs_HealthySuperficial.csv"))

