## =========================
## 1. Load packages
## =========================

library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

## =========================
## 2. Read data
## =========================

file_path <- "~/Desktop/Stool_OS_0729/Stool_OS_0729.xlsx"

dat <- read_excel(
  file_path,
  sheet = 1
)

print(colnames(dat))

## =========================
## 3. Data cleaning
## =========================

colnames(dat) <- trimws(colnames(dat))

required_cols <- c(
  "Group",
  "PFS_time",
  "PFS_status",
  "OS_time",
  "OS_status"
)

missing_cols <- setdiff(
  required_cols,
  colnames(dat)
)

if(length(missing_cols)>0){
  stop(
    paste(
      "Missing:",
      paste(missing_cols,collapse=", ")
    )
  )
}

dat <- dat %>%
  mutate(
    Group = trimws(as.character(Group)),
    PFS_time = as.numeric(PFS_time),
    PFS_status = as.numeric(PFS_status),
    OS_time = as.numeric(OS_time),
    OS_status = as.numeric(OS_status)
  )

## Three groups

dat <- dat %>%
  filter(
    Group %in% c(
      "Low",
      "Middle",
      "High"
    )
  )

## Low as reference

dat$Group <- factor(
  dat$Group,
  levels=c(
    "Low",
    "Middle",
    "High"
  )
)

cat("===== Group =====\n")
print(table(dat$Group))

## =========================
## 4. PFS / OS
## =========================

dat_pfs <- dat %>%
  filter(
    !is.na(Group),
    !is.na(PFS_time),
    !is.na(PFS_status),
    PFS_status %in% c(0,1),
    PFS_time>=0
  )

dat_os <- dat %>%
  filter(
    !is.na(Group),
    !is.na(OS_time),
    !is.na(OS_status),
    OS_status %in% c(0,1),
    OS_time>=0
  )

cat("===== PFS =====\n")
print(table(dat_pfs$Group))

cat("===== OS =====\n")
print(table(dat_os$Group))

## =========================
## 5. Kaplan-Meier
## =========================

fit_pfs <- survfit(
  Surv(PFS_time,PFS_status)~Group,
  data=dat_pfs
)

fit_os <- survfit(
  Surv(OS_time,OS_status)~Group,
  data=dat_os
)

## =========================
## 6. Cox HR
## =========================

cox_pfs <- coxph(
  Surv(PFS_time,PFS_status)~Group,
  data=dat_pfs
)

cox_os <- coxph(
  Surv(OS_time,OS_status)~Group,
  data=dat_os
)

get_cox_result <- function(model){
  s <- summary(model)
  data.frame(
    comparison = rownames(s$coefficients),
    HR = s$coefficients[,"exp(coef)"],
    low95 = s$conf.int[,"lower .95"],
    high95 = s$conf.int[,"upper .95"],
    p = s$coefficients[,"Pr(>|z|)"]
  )
}

cox_pfs_res <- get_cox_result(cox_pfs)
cox_os_res <- get_cox_result(cox_os)

print(cox_pfs_res)
print(cox_os_res)

make_hr_label <- function(df){
  paste0(
    gsub("^Group","",df$comparison),
    " vs Low: HR=",
    sprintf("%.2f",df$HR),
    " (",
    sprintf("%.2f",df$low95),
    "-",
    sprintf("%.2f",df$high95),
    "), p=",
    sprintf("%.3f",df$p),
    collapse="\n"
  )
}

label_pfs_hr <- make_hr_label(cox_pfs_res)
label_os_hr <- make_hr_label(cox_os_res)

## =========================
## 7. Pairwise log-rank
## =========================

pairwise_logrank <- function(data,time,status){
  formula <- as.formula(
    paste0(
      "Surv(",
      time,
      ",",
      status,
      ")~Group"
    )
  )

  res <- pairwise_survdiff(
    formula,
    data=data,
    p.adjust.method="BH"
  )

  return(res$p.value)
}

pair_pfs <- pairwise_logrank(
  dat_pfs,
  "PFS_time",
  "PFS_status"
)

pair_os <- pairwise_logrank(
  dat_os,
  "OS_time",
  "OS_status"
)

print(pair_pfs)
print(pair_os)

format_pairwise <- function(mat){
  txt <- c()
  for(i in seq_len(nrow(mat))){
    for(j in seq_len(ncol(mat))){
      if(!is.na(mat[i,j])){
        txt <- c(
          txt,
          paste0(
            rownames(mat)[i],
            " vs ",
            colnames(mat)[j],
            ": ",
            sprintf("%.3f",mat[i,j])
          )
        )
      }
    }
  }
  paste(
    txt,
    collapse="\n"
  )
}

label_pfs_pair <- format_pairwise(pair_pfs)
label_os_pair <- format_pairwise(pair_os)

label_pfs <- paste0(
  label_pfs_hr,
  "\n\nPairwise log-rank:\n",
  label_pfs_pair
)

label_os <- paste0(
  label_os_hr,
  "\n\nPairwise log-rank:\n",
  label_os_pair
)

## =========================
## 8. Colors
## =========================

group_levels <- c(
  "Low",
  "Middle",
  "High"
)

group_cols <- c(
  "#5BA8A2",
  "#6B8FAD",
  "#F2A099"
)

## =========================
## 9. PFS curve
## =========================

p1 <- ggsurvplot(
  fit_pfs,
  data = dat_pfs,
  pval = TRUE,              # Overall log-rank p-value for all three groups
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab = "Time (months)",
  ylab = "Progression-free survival probability",
  title = "PFS Kaplan-Meier Curve by Group",
  legend.title = "Group",
  legend.labs = group_levels,
  palette = group_cols,
  xlim = c(0, 38.5),
  break.time.by = 6
)

p1$plot <- p1$plot +
  scale_x_continuous(
    limits = c(0, 38.5),
    breaks = c(0, 6, 12, 18, 24, 30, 36),
    expand = c(0,0)
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.25,
    label = paste0(
      label_pfs_hr,
      "\n\n",
      "Pairwise log-rank:\n",
      label_pfs_pair
    ),
    hjust = 0,
    size = 3.6
  ) +
  theme(
    axis.text = element_text(
      color="black"
    ),
    axis.title = element_text(
      color="black"
    ),
    plot.title = element_text(
      hjust=0.5,
      face="bold"
    )
  )

p1$table <- p1$table +
  scale_x_continuous(
    limits = c(0,38.5),
    breaks = c(0,6,12,18,24,30,36),
    expand=c(0,0)
  )

print(p1)

## =========================
## 10. OS curve
## =========================

p2 <- ggsurvplot(
  fit_os,
  data = dat_os,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  title = "OS Kaplan-Meier Curve by Group",
  legend.title = "Group",
  legend.labs = group_levels,
  palette = group_cols,
  xlim = c(0,38.5),
  break.time.by = 6
)

p2$plot <- p2$plot +
  scale_x_continuous(
    limits=c(0,38.5),
    breaks=c(0,6,12,18,24,30,36),
    expand=c(0,0)
  ) +
  annotate(
    "text",
    x=5,
    y=0.25,
    label=paste0(
      label_os_hr,
      "\n\n",
      "Pairwise log-rank:\n",
      label_os_pair
    ),
    hjust=0,
    size=3.6
  ) +
  theme(
    axis.text = element_text(
      color="black"
    ),
    axis.title = element_text(
      color="black"
    ),
    plot.title = element_text(
      hjust=0.5,
      face="bold"
    )
  )

p2$table <- p2$table +
  scale_x_continuous(
    limits=c(0,38.5),
    breaks=c(0,6,12,18,24,30,36),
    expand=c(0,0)
  )

print(p2)