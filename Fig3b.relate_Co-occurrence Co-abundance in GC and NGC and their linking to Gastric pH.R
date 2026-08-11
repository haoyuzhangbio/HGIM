#####------------------------------------------------------------------------------------------

	##  Co-occurrence,Co-abundance
	## Differences in species co-abundance patterns between GC and NGC groups
	#### t_dat2 and j_dat2 were already matched within the same individuals. 
	### species has been filtered by 10% occurence in the corresponding body site.  
	##t_dat2 is one body site profile, such as tissue
	## j_dat2 is another body site profile, such as juice.

#####------------------------------------------------------------------------------------------

	getcor = function(t_dat2,tissueid, j_dat2,juiceid){


		##### species in GC samples
		tmp = grx2[match(colnames(j_dat2) ,grx2$ID),]

		a1 = t_dat2[,which(tmp$diagose=="Gastric cancer")]
		b1= j_dat2[,which(tmp$diagose=="Gastric cancer")]
		cancer = rownames(a1)[apply(a1,1,function(x){mean(x>0)})>0.01 & apply(b1,1,function(x){mean(x>0)})>0.01]

		##### species in NGC samples
		a1 = t_dat2[,which(tmp$diagose!="Gastric cancer")]
		b1= j_dat2[,which(tmp$diagose!="Gastric cancer")]
		noncancer = rownames(a1)[apply(a1,1,function(x){mean(x>0)})>0.01 & apply(b1,1,function(x){mean(x>0)})>0.01]


		##### species in both GC and NGC samples.
		kep = intersect(cancer,noncancer)
		kep =kep[kep%in%rownames(occ)]

		### co- abundance analysis
		c2 = cor(t(t_dat2[kep,which(tmp$diagose=="Gastric cancer")]),t(j_dat2[kep,which(tmp$diagose=="Gastric cancer")]),method="spearman")
		c2[is.na(c2)]=0


		c3 = cor(t(t_dat2[kep,which(tmp$diagose!="Gastric cancer")]),t(j_dat2[kep,which(tmp$diagose!="Gastric cancer")]),method="spearman")
		c3[is.na(c3)]=0



### fisher test the co-abundance difference between GC and NGC samples.
n1 <- sum(tmp$diagose == "Gastric cancer")
n2 <- sum(tmp$diagose != "Gastric cancer")

r_cancer <- diag(c2)
r_non_cancer <- diag(c3)

z_cancer <- 0.5 * log((1 + r_cancer) / (1 - r_cancer))
z_non_cancer <- 0.5 * log((1 + r_non_cancer) / (1 - r_non_cancer))

se_cancer <- 1 / sqrt(n1 - 3)
se_non_cancer <- 1 / sqrt(n2 - 3)

z_scores <- (z_cancer - z_non_cancer) / sqrt(se_cancer^2 + se_non_cancer^2)
p_values <- 2 * pnorm(-abs(z_scores))
p_adjusted <- p.adjust(p_values, method = "BH")

# Combine with your original d_plot_rank_location
d_plot_rank_location <- data.frame(
    Stool = diag(c2),      # Cancer group correlations
    Device = diag(c3),     # Non-cancer group correlations
    Symbol = kep,          # Species names
    species = kep,          # Species names
    diff = diag(c2) - diag(c3),  # Raw difference
    z_score = z_scores,    # Fisher's z-test score
    p_value = p_values,    # Unadjusted p-value
    p_adj = p_adjusted    # Adjusted p-value (BH correction)
) %>%
  mutate(
	abs_diff = abs(diff),
	  rank = rank(-abs_diff),
    # Label based on correlation difference and magnitude
    label = ifelse(
      rank <= 20 & abs(diff) >= 0.2 & (Stool >= 0.2 | Device > 0.2) & p_adj<0.05,
	# p_adj<1e-10,
      TRUE,
      FALSE
    ),
    # Color based on direction of difference
    color = ifelse(
      diff > 0.1 & p_adj<0.05, '#e76254',  # Stronger in cancer (Stool)
      ifelse(diff <= -0.1 & p_adj<0.05, '#376795', NA)  # Stronger in non-cancer (Device)
    ),
    # Significance flag (optional: e.g., p_adj < 0.05)
    is_significant = ifelse(p_adj < 0.05, TRUE, FALSE)
  )

	d_plot_rank_location_sp =d_plot_rank_location

		### jaccard index in NGC
		
		library(vegan)
		a1 = t_dat2[rownames(t_dat2)%in%rownames(occ),which(tmp$diagose!="Gastric cancer")]
		b1= j_dat2[rownames(j_dat2)%in%rownames(occ),which(tmp$diagose!="Gastric cancer")]


		rownames(a1)= paste0("a|",rownames(a1))
		rownames(b1)= paste0("b|",rownames(b1))
		colnames(a1)=colnames(b1)
		c1 = rbind(a1,b1)
		c1[c1>0]=1
		

		jc2 = vegdist(c1, method="jaccard") 
		jc2= as.matrix(jc2)
		jc2_1 = 1- jc2[match(rownames(a1),rownames(jc2)), match(gsub("^a","b",rownames(a1)),colnames(jc2))]
		jc2_1[1:4,1:4]
	

		### jaccard index in GC

		a1 = t_dat2[rownames(t_dat2)%in%rownames(occ),which(tmp$diagose=="Gastric cancer")]
		b1= j_dat2[rownames(j_dat2)%in%rownames(occ),which(tmp$diagose=="Gastric cancer")]

		a1 = df_r[match(rownames(a1),rownames(df_r)), match(colnames(a1),colnames(df_r))]
		b1 = df_r[match(rownames(b1),rownames(df_r)), match(colnames(b1),colnames(df_r))]


		rownames(a1)= paste0("a|",rownames(a1))
		rownames(b1)= paste0("b|",rownames(b1))
		colnames(a1)=colnames(b1)
		c1 = rbind(a1,b1)
		c1[c1>0]=1
		jc3 = vegdist(c1, method="jaccard") 
		jc3= as.matrix(jc3)
		jc3_1 = 1- jc3[match(rownames(a1),rownames(jc3)), match(gsub("^a","b",rownames(a1)),colnames(jc3))]
		jc3_1[1:4,1:4]
		

		dim(jc2_1)
		dim(jc3_1)


##### glm test jaccard similarity difference between GC and NGC
library(parallel)
library(dplyr)

process_species <- function(species) {
  tryCatch({
    # Get variables from parent environment
    t_dat2 <- get("t_dat2", envir = parent.frame())
    j_dat2 <- get("j_dat2", envir = parent.frame())
    tmp <- get("tmp", envir = parent.frame())
    
    clean_species <- gsub("^a\\|", "", species)
    
    stool_presence <- t_dat2[clean_species, ] > 0
    device_presence <- j_dat2[clean_species, ] > 0
    
    glm_data <- data.frame(
      y = as.numeric(stool_presence & device_presence),
      group = ifelse(tmp$diagose == "Gastric cancer", 1, 0)
    )
    
    if (length(unique(glm_data$y)) <= 1) {
      return(data.frame(
        Symbol = clean_species,
        beta_group = NA,
        p_value = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    glm_fit <- glm(y ~ group, data = glm_data, family = binomial)
    
    data.frame(
      Symbol = clean_species,
      beta_group = coef(glm_fit)["group"],
      p_value = coef(summary(glm_fit))["group", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Symbol = gsub("^a\\|", "", species),
      beta_group = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    )
  })
}

# Set up parallel cluster
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export required variables to all workers
clusterExport(cl, c("t_dat2", "j_dat2", "tmp"), envir = environment())

# Run parallel processing
species_list <- gsub("^a\\|", "", rownames(jc3_1))
results_list <- parLapply(cl, species_list, function(x) {
  process_species(x)
})

# Clean up
stopCluster(cl)

# Combine results
glm_results <- bind_rows(results_list)
glm_results$p_adj <- p.adjust(glm_results$p_value, method = "BH")
glm_results[grepl("nucl",glm_results[,1]),]
library(dplyr)
library(vegan)


	   rhocut1 =0.1
		rhocut=0.25
# Merge GLM results with your original d_plot_rank_location
d_plot_rank_location <- data.frame(
  Stool = diag(jc3_1),        # Jaccard similarity (cancer)
  Device = diag(jc2_1),       # Jaccard similarity (non-cancer)
  Symbol = gsub("^a\\|", "", rownames(jc3_1)),
  species = gsub("^a\\|", "", rownames(jc3_1))
) %>%
  left_join(glm_results, by = "Symbol") %>%
  mutate(
	abs_diff = abs(Stool - Device),
	  rank = rank(-abs_diff),
    diff = Stool - Device,
    # Label based on difference, magnitude, AND significance
    label = ifelse(
        rank <= 20 & abs(diff) > rhocut & (Stool >= 0.2 | Device > 0.2) & p_adj < 0.05,
      TRUE,
      FALSE
    ),
    # Color by direction of difference (only if significant)
    color = ifelse(
      diff > rhocut1 & p_adj < 0.05 , '#e76254',  # Stronger in cancer
      ifelse(diff <= -rhocut1  & p_adj < 0.05, '#376795', NA)  # Stronger in control
    )
  )




d_plot_rank_location_jaccard = d_plot_rank_location




	write.csv(d_plot_rank_location_jaccard,paste0(outpath,paste0(tissueid," and ", juiceid),".jaccard.similarity.csv"))

	write.csv(d_plot_rank_location_sp,paste0(outpath,paste0(tissueid," and ", juiceid),".spearman.csv"))

	}


#####------------------------------------------------------------------------------------------

	##  gastric pH distributions were compared between samples with and without species sharing
	## Differences in species co-abundance patterns between GC and NGC groups
	#### t_dat2 and j_dat2 were already matched within the same individuals. 
	### species has been filtered by 10% occurence in the corresponding body site.  

#####------------------------------------------------------------------------------------------


PH_jaccard = function(t_dat2,tissueid, j_dat2,juiceid){

		tmp = grx2[match(colnames(j_dat2) ,grx2$ID),]

		j_dat2[j_dat2>0]=1
		t_dat2[t_dat2>0]=1



results <- apply(j_dat2, 1, function(species_oral) {
  # Get the corresponding gut data for the same species (rows are aligned)
  row_idx <- which(apply(j_dat2, 1, identical, species_oral))
  species_gut <- t_dat2[row_idx, ]
  
  # Define groups
  group1 <- species_oral == 1 & species_gut == 1  # Both 1 (A ∩ B)
  group2 <- (species_oral == 1 | species_gut == 1) & !group1  # Either 1 (A ∪ B \ A ∩ B)
  
  # Extract pH values (exclude NAs)
  pH_group1 <- na.omit(tmp$PH[group1])
  pH_group2 <- na.omit(tmp$PH[group2])
  
  # Skip if insufficient data
  if (length(pH_group1) < 20 || length(pH_group2) < 20) {
    return(NULL)  # Skip instead of returning partial results
  }
  
  # Wilcoxon test
  test <- wilcox.test(pH_group1, pH_group2, exact = FALSE)
  
  # Calculate metrics
  median1 <- median(pH_group1)
  median2 <- median(pH_group2)
  fold_change <- median1 / median2
  
  # Effect size (Cliff’s Delta)
  library(effsize)
  cliff <- cliff.delta(pH_group1, pH_group2)
  
  # Return as a NAMED LIST (critical for correct DF conversion)
  list(
    p.value = test$p.value,
    W = as.numeric(test$statistic),  # Ensure numeric (not "htest" object)
    median_group1 = median1,
    median_group2 = median2,
    fold_change = fold_change,
    cliff_delta = cliff$estimate,
    n_group1 = length(pH_group1),
    n_group2 = length(pH_group2)
  )
})

# Remove NULL entries (species with insufficient data)
results <- results[!sapply(results, is.null)]

# Convert to data frame (with proper column names)
results_df <- do.call(rbind, lapply(results, as.data.frame))
results_df$species <- names(results)  # Add species names from the original list

# Multiple testing correction
results_df$p.adj <- p.adjust(results_df$p.value, method = "BH")
results_df$testgroup = paste0(tissueid,"-",juiceid)

results_df

}




	