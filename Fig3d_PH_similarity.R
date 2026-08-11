############################################################
## Fig.3C PH vs similarity
## Spearman rho + P
############################################################

library(ggplot2)
library(dplyr)
library(grid)
library(ggpubr)

############################################################
## 1. Read data
############################################################

intra_pd1 <- read.csv(
  "~/Desktop/正刊可用图_20260716/Fig3a_20260719/intra_pd1.csv",
  stringsAsFactors = FALSE
)

############################################################
## 2. similarity index transformation
############################################################

intra_pd1$similarity_raw <-
  1 - intra_pd1$bray_distance

intra_pd1$similarity_index <-
  log(
    intra_pd1$similarity_raw + 0.01
  )

############################################################
## 3. type_combination order
############################################################

nonfecal_order <- c(
  "Gjuice_tongue",
  "Gtissue_tongue",
  "Gjuice_Gtissue"
)

fecal_order <- c(
  "fecal_tongue",
  "fecal_Gjuice",
  "fecal_Gtissue"
)

all_order <- c(
  nonfecal_order,
  fecal_order
)

############################################################
## 4. Data cleaning
############################################################

plot_data_all <- intra_pd1 %>%
  filter(
    type_combination %in% all_order
  ) %>%
  filter(
    !is.na(PH),
    !is.na(similarity_index)
  ) %>%
  mutate(
    type_combination =
      factor(
        type_combination,
        levels = all_order
      )
  )

if(nrow(plot_data_all)==0){
  stop(
    "No matching type_combination found"
  )
}

############################################################
## 5. Spearman rho + P
############################################################

get_spearman <- function(data){
  data %>%
    group_by(type_combination) %>%
    summarise(
      test = list(
        cor.test(
          PH,
          similarity_index,
          method="spearman",
          exact=FALSE
        )
      ),
      .groups="drop"
    ) %>%
    mutate(
      rho =
        sapply(
          test,
          function(x)
            as.numeric(x$estimate)
        ),
      p =
        sapply(
          test,
          function(x)
            x$p.value
        ),
      label =
        paste0(
          "\u03C1=",
          round(rho,2),
          "\nP=",
          format.pval(
            p,
            digits=2,
            eps=0.001
          )
        )
    )
}

cor_df <- get_spearman(plot_data_all)

print(cor_df)

############################################################
## 6. Disease colors
############################################################

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

############################################################
## 7. Plotting function
############################################################

plot_PH_similarity <- function(data, order_use, title_text){
  data <- data %>%
    filter(
      type_combination %in% order_use
    ) %>%
    mutate(
      type_combination =
        factor(
          type_combination,
          levels=order_use
        )
    )

  cor_use <- cor_df %>%
    filter(
      type_combination %in% order_use
    )

  ggplot(
    data,
    aes(
      x=PH,
      y=similarity_index
    )
  )+
    geom_point(
      aes(
        color=diagose
      ),
      size=1.8,
      alpha=0.9
    )+
    geom_smooth(
      method="glm",
      color="black",
      fill="grey70",
      se=TRUE,
      linewidth=0.6
    )+
    geom_text(
      data=cor_use,
      aes(
        x=-Inf,
        y=Inf,
        label=label
      ),
      inherit.aes=FALSE,
      hjust=-0.1,
      vjust=1.3,
      size=3.2
    )+
    facet_wrap(
      ~type_combination,
      ncol=3,
      scales="fixed"
    )+
    scale_color_manual(
      values=box_colors,
      name=NULL
    )+
    labs(
      x="PH",
      y="log(1 - Bray distance + 0.01)",
      title=title_text
    )+
    theme_classic(base_size=12)+
    theme(
      strip.background=
        element_rect(
          fill="white",
          color="black"
        ),
      strip.text=
        element_text(
          size=12
        ),
      axis.text=
        element_text(
          color="black"
        ),
      axis.title=
        element_text(
          size=12,
          face="bold"
        ),
      legend.position="right"
    )
}

############################################################
## 8. Two plots
############################################################

p_nonfecal <- plot_PH_similarity(
  plot_data_all,
  nonfecal_order,
  "Non-fecal samples"
)

p_fecal <- plot_PH_similarity(
  plot_data_all,
  fecal_order,
  "Fecal-related samples"
)

############################################################
## 9. Output
############################################################

print(p_nonfecal)
print(p_fecal)

ggsave(
  "~/Desktop/Fig3C_PH_similarity_nonfecal.pdf",
  p_nonfecal,
  width=10,
  height=4
)

ggsave(
  "~/Desktop/Fig3C_PH_similarity_fecal.pdf",
  p_fecal,
  width=10,
  height=4
)

############################################################
## Export Spearman correlation results
############################################################

cor_output <- cor_df %>%
  select(
    type_combination,
    rho,
    p,
    label
  )

write.csv(
  cor_output,
  "~/Desktop/Fig3C_PH_similarity_spearman_results.csv",
  row.names = FALSE
)

print(cor_output)