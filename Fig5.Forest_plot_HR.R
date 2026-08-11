library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(patchwork)

rank_file <- "~/Desktop/Stool_OS_0729/species_survival_ranking_adjusted_Gender_Age.xlsx"

total_n <- 20
target_each <- total_n / 2
fdr_cutoff <- 0.05

risk_sig <- read_excel(rank_file, sheet = "Risk_species") %>%
  filter(!is.na(HR), !is.na(HR_low95), !is.na(HR_high95), !is.na(FDR)) %>%
  filter(FDR < fdr_cutoff) %>%
  mutate(group = "Risk species") %>%
  arrange(FDR, desc(abs_z), desc(c_index))

protect_sig <- read_excel(rank_file, sheet = "Protective_species") %>%
  filter(!is.na(HR), !is.na(HR_low95), !is.na(HR_high95), !is.na(FDR)) %>%
  filter(FDR < fdr_cutoff) %>%
  mutate(group = "Protective species") %>%
  arrange(FDR, desc(abs_z), desc(c_index))

n_risk <- nrow(risk_sig)
n_protect <- nrow(protect_sig)

risk_take <- min(target_each, n_risk)
protect_take <- min(target_each, n_protect)

remaining <- total_n - risk_take - protect_take

if (remaining > 0) {
  if (n_risk - risk_take >= n_protect - protect_take) {
    risk_take <- min(n_risk, risk_take + remaining)
  } else {
    protect_take <- min(n_protect, protect_take + remaining)
  }
}

risk_df <- risk_sig %>%
  slice_head(n = risk_take) %>%
  arrange(desc(HR)) %>%
  mutate(
    species_clean = str_replace_all(species, "_", " "),
    species_clean = factor(species_clean, levels = species_clean),
    species_clean = fct_rev(species_clean),
    HR_label = sprintf("%.2f (%.2f–%.2f)", HR, HR_low95, HR_high95)
  )

protect_df <- protect_sig %>%
  slice_head(n = protect_take) %>%
  arrange(HR) %>%
  mutate(
    species_clean = str_replace_all(species, "_", " "),
    species_clean = factor(species_clean, levels = species_clean),
    species_clean = fct_rev(species_clean),
    HR_label = sprintf("%.2f (%.2f–%.2f)", HR, HR_low95, HR_high95)
  )

if (nrow(risk_df) == 0 & nrow(protect_df) == 0) {
  stop("No species with FDR < 0.05 available for plotting. Please relax the threshold or check input data.")
}

#############################
## Risk species panel parameters
#############################
risk_xlim <- c(0.98, 1.45)
risk_breaks <- c(1.0, 1.2, 1.4)
risk_label_x <- 1.47

#############################
## Protective species panel parameters
#############################
protect_xlim <- c(0.55, 1.02)
protect_breaks <- c(0.6, 0.8, 1.0)
protect_label_x <- 1.04

#############################
## Risk species panel
#############################
p_risk <- ggplot(risk_df, aes(x = HR, y = species_clean)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.6) +
  geom_errorbarh(
    aes(xmin = HR_low95, xmax = HR_high95),
    height = 0,
    linewidth = 0.7,
    color = "#e76254"
  ) +
  geom_point(
    size = 3.2,
    color = "#e76254"
  ) +
  geom_text(
    aes(x = risk_label_x, label = HR_label),
    hjust = 0,
    size = 3.6,
    color = "black"
  ) +
  annotate(
    "text",
    x = risk_label_x,
    y = length(levels(risk_df$species_clean)) + 0.9,
    label = "HR (95% CI)",
    hjust = 0,
    fontface = "bold",
    size = 3.8
  ) +
  scale_x_continuous(
    breaks = risk_breaks,
    labels = c("1.0", "1.2", "1.4"),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = risk_xlim, clip = "off") +
  labs(
    title = "Risk species",
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic", size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#e76254"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 95, 5.5, 5.5)
  )

#############################
## Protective species panel
#############################
p_protect <- ggplot(protect_df, aes(x = HR, y = species_clean)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.6) +
  geom_errorbarh(
    aes(xmin = HR_low95, xmax = HR_high95),
    height = 0,
    linewidth = 0.7,
    color = "#376795"
  ) +
  geom_point(
    size = 3.2,
    color = "#376795"
  ) +
  geom_text(
    aes(x = protect_label_x, label = HR_label),
    hjust = 0,
    size = 3.6,
    color = "black"
  ) +
  annotate(
    "text",
    x = protect_label_x,
    y = length(levels(protect_df$species_clean)) + 0.9,
    label = "HR (95% CI)",
    hjust = 0,
    fontface = "bold",
    size = 3.8
  ) +
  scale_x_continuous(
    breaks = protect_breaks,
    labels = c("0.6", "0.8", "1.0"),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = protect_xlim, clip = "off") +
  labs(
    title = "Protective species",
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic", size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "#376795"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 95, 5.5, 5.5)
  )

#############################
## Side-by-side layout
#############################
p_final <- p_risk + p_protect + plot_layout(ncol = 2)

p_final

ggsave(
  "~/Desktop/Nature_like_HR_forest_top20_FDR005_species_split.pdf",
  p_final,
  width = 12,
  height = 6
)