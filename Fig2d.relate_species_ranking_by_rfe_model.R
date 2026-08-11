ROC_with_label <- function(df_auc, highlight_points = NULL, opt_size = NULL, startfrom = 0.5) {
  colnames(df_auc) <- c("Variables", "AUC", "ROCSD", "Sens", "Spec")
  

  df_auc$label <- if (!is.null(highlight_points)) {
    ifelse(df_auc$Variables %in% highlight_points,
           paste0(sprintf("%.3f", df_auc$AUC), "(", round(df_auc$Variables, 0), ")"),
           NA)
  } else {
    NA
  }
  

  

  if (is.null(opt_size)) {
    opt_size <- df_auc$Variables[which.max(df_auc$AUC)]
  }
  print(df_auc)
  print(opt_size)
  stats_row <- df_auc[df_auc$Variables == opt_size, ]
  
  auc_text <- paste0("AUC = ", sprintf("%.4f", stats_row$AUC))
  sens_text <- paste0("Sensitivity: ", sprintf("%.1f%%", stats_row$Sens * 100))
  spec_text <- paste0("Specificity: ", sprintf("%.1f%%", stats_row$Spec * 100))

  

  p <- ggplot(df_auc, aes(x = Variables, y = AUC)) +

    geom_ribbon(aes(ymin = pmax(AUC - ROCSD, startfrom), ymax = pmin(AUC + ROCSD, 1)),
                fill = "#2c7bb6", alpha = 0.2) +
    

    geom_line(linewidth = 1.2, color = "#2c7bb6") +
    

    geom_point(data = subset(df_auc, !is.na(label)), shape = 1, size = 2) +
    geom_text_repel(
      data = subset(df_auc, !is.na(label)),
      aes(label = label),
      size = 3,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.4,
      segment.color = "grey50",
      segment.size = 0.4,
      segment.curvature = -0.1
    ) +
    
 
    annotate("text", x = max(df_auc$Variables)*0.75, y = startfrom + 0.03, 
             label = paste(auc_text, sens_text, spec_text,  sep = "\n"),
             hjust = 0, vjust = 0, size = 4.2, color = "black") +
    
    labs(
      x = "Number of Features",
      y = "ROC AUC Score",
      title = "Recursive Feature Elimination Performance",
      subtitle = "Model performance vs. number of selected features"
    ) +
    
    scale_y_continuous(limits = c(startfrom, 1), breaks = pretty_breaks(6)) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.title = element_text(face = "bold"),
      # plot.margin = margin(1, 1, 1, 1, "cm")
      plot.margin = unit(c(1, 1, 1, 1), "cm")  
    )
  
  return(p)
}




dir_path <- "tissue_select_filter_200/"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
setwd(dir_path)

prefix="tissue_select_filter_200"



total_df = as.data.frame(readRDS("../tissue_mag.total_df.rds"))
rownames(total_df) = total_df$samples
library(openxlsx)

#### read the tissue-associated microbial species
species_list = read.xlsx("Tissue_93species.xlsx")
selected_species =  species_list
length(intersect(colnames(total_df), species_list$Species))

### keep tissue-associated microbial species
total_df = as.data.frame(total_df[,c("samples", "group", intersect(colnames(total_df), species_list$Species))])

#___________________________________________
con_rfe = as.data.frame(total_df$group)
rownames(con_rfe) = rownames(total_df)
colnames(con_rfe) = c("group")
total_df_tmp = total_df
total_df_tmp$group <- NULL
total_df_tmp$samples <- NULL
rownames(total_df_tmp) = total_df$samples
rownames(con_rfe) = total_df$samples
# 
library(dplyr)
# total_df_tmp <- total_df_tmp[, -grep("chm13|UniVec", colnames(total_df_tmp))]
total_df_tmp[] <- lapply(total_df_tmp, function(x) as.numeric(trimws(x)))
# total_df_tmp <- total_df_tmp %>% select(where(~ mean(. > 0, na.rm = TRUE) > 0.05))








fiveStats <- function(...) c(twoClassSummary(...),
                             defaultSummary(...))

library(caret)
predVars <- colnames(total_df_tmp)
varSeq <- c(c(1,2,3,4,5,6,7,8,9),seq(10, 90, by = 5), length(predVars)-1)
# varSeq <- c(1,2,3,4,5)
varSeq

## The built-in random forest functions are in rfFuncs.
str(rfFuncs)
newRF <- rfFuncs
newRF$summary <- fiveStats
set.seed(1120)
index <- createMultiFolds(factor(con_rfe$group,levels=c("cancer","control")), times = 5)

## The control function is similar to trainControl():
ctrl <- rfeControl(method = "repeatedcv",
                   repeats = 5,
                   verbose = TRUE,
                   functions = newRF,
                   index = index,
                   saveDetails = TRUE
                   # ,
                   # allowParallel = TRUE 
)

set.seed(1120)
rfRFE <- rfe(x = total_df_tmp,
             y = factor(con_rfe$group,levels=c("cancer","control")),
             sizes = varSeq,
             metric = "ROC",
             rfeControl = ctrl,
             sampsize = c(min(table(con_rfe$group)), min(table(con_rfe$group))),
             strata = con_rfe$group,
             ntree = 1000)


# stopCluster(cl)
# registerDoSEQ()


library(openxlsx)
library(dplyr)
rfRFE$results[,1:2]
write.xlsx(rfRFE$results, paste0(format(Sys.Date(), "%Y%m%d"), ".", prefix, ".rfe_result.xlsx"))



species_importance <- rfRFE$variables %>%
  group_by(var) %>%
  summarise(
    Mean_Importance = mean(Overall),
    SD_Importance = sd(Overall),
    N_Folds = n()
  ) %>%
  arrange(desc(Mean_Importance)) 
head(species_importance, 10)
write.xlsx(species_importance, paste0(format(Sys.Date(), "%Y%m%d"), ".", prefix, ".species_importance.xlsx"))





library(scales)
library("ggrepel")

highlight_point <- c(1,2,3,4,5,6,7,8,9,10,35,60,100,350,1000,5000)
df_auc_input <- rfRFE$results[, c("Variables", "ROC", "ROCSD", "Sens", "Spec")]

roc_with_label <- ROC_with_label(df_auc_input, highlight_point, opt_size = 92, startfrom = 0.7)
roc_with_label
ggsave(plot = roc_with_label,  filename = paste0(format(Sys.Date(), "%Y%m%d"), ".", prefix, ".lineplot_with_label.withband.pdf"), width = 8, height = 5, dpi=600)


highlight_point <- varSeq
df_auc_input <- rfRFE$results[, c("Variables", "ROC", "ROCSD", "Sens", "Spec")]
roc_with_label <- ROC_with_label(df_auc_input, highlight_point, opt_size = 92, startfrom = 0.7)
roc_with_label
ggsave(plot = roc_with_label,  filename = paste0(format(Sys.Date(), "%Y%m%d"), ".", prefix, ".total.lineplot_with_label.withband.pdf"), width = 8, height = 5, dpi=600)

save.image(paste0(format(Sys.Date(), "%Y%m%d"), ".", prefix, ".Rdata"),)

