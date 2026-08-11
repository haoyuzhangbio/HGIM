############################################################
## Figure1f Three-site Venn diagram
## Oral / Stomach / Gut MAG overlap
############################################################

library(readr)
library(dplyr)
library(VennDiagram)
library(grid)

############################################################
## 1. Input file
############################################################

input_file <- "~/Desktop/Figure1f.data.csv"

out_pdf <- "~/Desktop/Figure1f_three_site_Venn.pdf"

############################################################
## 2. Read data
############################################################

df <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(df)[is.na(colnames(df)) | colnames(df)==""] <- "MAG_ID"
head(df)
colnames(df)

############################################################
## 3. Define three body sites
############################################################

# Stomach:
# Gjuice or Gtissue any non-zero

df <- df %>%
  mutate(
    Stomach = ifelse(
      Gjuice > 0 | Gtissue > 0,
      TRUE,
      FALSE
    ),
    Gut = ifelse(
      fecal > 0,
      TRUE,
      FALSE
    ),
    Oral = ifelse(
      tongue > 0,
      TRUE,
      FALSE
    )
  )

############################################################
## 4. Count each set size
############################################################

cat("Stomach MAG:", sum(df$Stomach), "\n")
cat("Gut MAG:", sum(df$Gut), "\n")
cat("Oral MAG:", sum(df$Oral), "\n")

############################################################
## 5. Build Venn input
############################################################

venn_list <- list(
  Stomach = which(df$Stomach),
  Gut = which(df$Gut),
  Oral = which(df$Oral)
)

############################################################
## 6. Output Venn diagram
############################################################

pdf(
  out_pdf,
  width = 6,
  height = 6
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c(
    "#376795",   # Stomach blue
    "#E76254",   # Gut red
    "#299D8F"    # Oral green
  ),
  alpha = 0.45,
  lwd = 1.5,
  cex = 1.4,
  cat.cex = 1.3,
  cat.fontface = "bold",
  fontface = "bold",
  category.names = c(
    "Stomach",
    "Gut",
    "Oral"
  ),
  margin = 0.1
)

grid.draw(venn.plot)

dev.off()

cat("Done:", out_pdf, "\n")