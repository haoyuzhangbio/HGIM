


library("tidyverse")
library(dplyr)
library(ggsci)

library(dplyr)
library(ggplot2)
library(ggplot2)
library(ggsci)

library(ggbeeswarm)
library(ggessentials)

colors <- pal_npg("nrc")(10)
# Define the function
library("tidyverse")
library("coin")
library("pROC")
library("yaml")
log.n0 =1e-5


	### read phenotype

	grx2 = read.csv("D:/Work/GAC2022/Population/Geneset/final.config.2025.11.17.csv",head=T)

	dim(grx2)

	table(grx2$diagose)
	grx2$diagose  = factor(grx2$diagose,levels=c("Gastric polyps","Gastric ulcer",
	"Healthy stomach","Superficial gastritis","Atrophic gastritis","Intestinal metaplasia",	
	"Intraepithelial Neoplasia","Gastric cancer"))

		table(grx2$group2)
		grx2$group2 = as.character(grx2$diagose)
		grx2$group2[which(grx2$TNM.type=="I")]="Gastric cancer I-II"
		grx2$group2[which(grx2$TNM.type=="II")]="Gastric cancer I-II"
		grx2$group2[which(grx2$TNM.type=="III")]="Gastric cancer III-IV"
		grx2$group2[which(grx2$TNM.type=="IV")]="Gastric cancer III-IV"
		
		table(grx2[,c("group1","diagose")],useNA ="always")
		 
		disease =read.csv("C:\\Users\\jiezhuye\\Downloads\\data\\config_disease_0919.csv",as.is=T,head=T,row.names=1)
	      table(disease[,2],disease[,3],useNA ="always")
		grx2=data.frame(grx2,disease[match(grx2$ID,disease$ID),-1])

		apply(grx2[,c("age","gender","BMI","smoke","drink","disease_HTN","disease_T2D","disease_CVD")],2,function(x){ mean(!is.na(x)) })
		
		


		grx2_org = grx2  ### for gc vs non gc test
		grx2$group2 = as.character(grx2$group2)

			table(grx2[,c("group","group2")])
		
		table(grx2$group2,grx2$sampletype)
		dim(grx2)
		dim(grx2_org)

		### fill missing with median
# Define variables
vars <- c("age","gender","BMI","smoke","drink","disease_HTN","disease_T2D","disease_CVD")

# Impute missing values
for(var in vars) {
  if(is.numeric(grx2[[var]])) {
    # For numeric variables (age, BMI, disease indicators)
    grx2[is.na(grx2[,var]), var] <- median(grx2[,var], na.rm = TRUE)
  } else if(is.character(grx2[[var]])) {
    # For categorical variables (gender, smoke, drink)
    # Use mode (most frequent value) instead of median
    mode_val <- names(sort(table(grx2[,var]), decreasing = TRUE))[1]
    grx2[is.na(grx2[,var]), var] <- mode_val
  }
}


readx = function(file){
		library(data.table)
	 #d1 = read.table(file,sep="\t",row.names=1,head=T)
	d1<-fread(file,sep="\t",head=T)
	d1 = as.data.frame(d1)
	rownames(d1)=d1[,1]
	d1=d1[,-1]


	
	colnames(d1) = gsub(".mp4","",gsub("^X","",colnames(d1)))
	# keep species
	d1  = d1[grepl("s__",rownames(d1)) & (!grepl("t__",rownames(d1))),]
	rownames(d1) = unlist(lapply(strsplit(rownames(d1),"|",fixed=T),function(x){x[length(x)]}))
	#d1=d1[,apply(d1,2,sum)>0]
	print(dim(d1))	 #4789
	print(d1[1:4,1:4]	)
	d1
	
	}



######## read data


df = data.table::fread("C:\\Users\\jiezhuye\\Downloads\\data\\braken_data\\20260715.total_bracken_profile.relative.txt",sep="\t")

	df = as.data.frame(df)
	rownames(df)=df[,1]
	df=df[,-1]
	df=as.data.frame(t(df))
	df_org= df ## for gc vs non gc test

	nm = intersect(colnames(df),grx2$ID)
	length(nm)	
	dim(grx2)
	grx2 = grx2[match(nm,grx2$ID),]
	df =df[,match(nm,colnames(df))]
	dim(df)
	dim(df_org)

	rownames(df)[grepl("pylor",rownames(df))]

	

	
###################### read count
df_r = data.table::fread("C:\\Users\\jiezhuye\\Downloads\\data\\braken_data\\20260715.total_bracken_profile.reads.txt",sep="\t")

	df_r = as.data.frame(df_r)
	rownames(df_r)=df_r[,1]
	df_r=df_r[,-1]
	df_r =as.data.frame(t(df_r))

	library(openxlsx)
	reads = read.xlsx("C:\\Users\\jiezhuye\\Downloads\\data\\braken_data\\MP4.mapping.summary.xlsx",sheet="Sheet 1")

	rc_org = apply(df_r,2,sum)
	sum(is.na(match(grx2_org$ID,reads$id)))
	match(grx2_org$ID[is.na(match(grx2_org$ID,reads$id))],colnames(df_r))
	grx2_org$rc =  reads[match(grx2_org$ID,reads$id),"reads_processed" ]
	grx2$rc =  reads[match(grx2$ID,reads$id),"reads_processed" ]
	grx2$rc[which(grx2$rc==0)]=1
	grx2_org$rc[which(grx2_org$rc==0)]=1
	range(grx2_org$rc,na.rm=T)
	


#library(phyloseq)
#library(ANCOMBC)


################## filtre species with less than 10% occurence
	table(grx2$sampletype)
	df_r =df_r[,match(grx2$ID,colnames(df_r))]
	oral=df_r[,which(grx2$sampletype=="tongue")]
	gut=df_r[,which(grx2$sampletype=="fecal")]
	tdat=df_r[,which(grx2$sampletype=="Gtissue")]
	jdat=df_r[,which(grx2$sampletype=="Gjuice")]

	occ=data.frame(fecal=apply(gut,1,function(x){mean(x>200)}),
		Gjuice=apply(jdat,1,function(x){mean(x>200)}),
		Gtissue = apply(tdat,1,function(x){mean(x>200)}),
		tongue=apply(oral,1,function(x){mean(x>200)}) )
			occ_org = occ
	occ = occ[apply(occ,1,function(x){ any(x>0.1) }),]

occ_gtdb = taxo[match(rownames(occ),taxo$sp),"GTDB_TAXONOMY"]

library(dplyr)

taxo1 <- taxo %>%
  group_by(sp) %>%
  summarise(max_MAG_COUNT = max(MAG_COUNT, na.rm = TRUE))
taxo1 =as.data.frame(taxo1)
	occ$MAG_COUNT = taxo1[match(rownames(occ),taxo1[,1]),"max_MAG_COUNT"]
	occ =occ[occ$MAG_COUNT>=5,] #2509
	table( occ[,"fecal"]>0.1 )
	table( occ[,"tongue"]>0.1 )
	table( occ[,"Gjuice"]>0.1 )




### sgb taxo
	taxo= data.table::fread("C:\\Users\\jiezhuye\\Desktop\\GC\\complete_50.sgb_number_taxonomy.clean.fix_space.tax_sim.rename.tsv",sep="\t")
	taxo=as.data.frame(taxo)
	taxo[grepl("SGB4567",taxo[,4]),]
	taxo$sp =  gsub("s__","",unlist(lapply(strsplit(taxo[,4],"|",fixed=T),function(x){x[length(x)]})))
	head(taxo)



	

##-----------------------------------------------
## glm and wilcox test-----  stage vs healthy
##-------------------------------------------------
	
		#############


 stat1 = function(x1,conf,g1,g2){ 
    x <- as.numeric( x1[which(conf$group2==g2)] )
    y <- as.numeric( x1[which(conf$group2==g1)] )

    
    # Wilcoxon
    p.val <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # AUC
    aucs.all <- c(roc(controls=y, cases=x, 
                                direction='<', ci=TRUE, auc=TRUE)$ci)
    aucs.mat<- c(roc(controls=y, cases=x, 
                           direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    
    # FC

	q.p <- quantile(x, probs=seq(.1, .9, .05))
	q.n <- quantile(y, probs=seq(.1, .9, .05))

    fc<- sum(q.p - q.n)/length(q.p)
	newd = c(x,y)
	newg = c(rep(1,length(x)), rep(0,length(y)))
	tmp = by(rank(newd), newg, mean,na.rm=T) 
	 
	# log rank fc
	fc_rank = log(tmp[2]/tmp[1])
 	abun = mean(c(x,y),na.rm=T)


	formula = "group2 ~ age+gender+target"
	d= data.frame(target=as.numeric(scale(x1)) ,
			conf[,c("group2","age","gender")] );
	d$group2 = factor(d$group2,levels=c(g1,g2))

	u = summary(glm(as.formula(formula), data=na.omit(d), family=binomial() ))
	bx = coef(u)[grepl("target",rownames(coef(u))),]

	formula = "group2 ~ age+gender+BMI+smoke+drink+rc+disease_HTN+disease_T2D+disease_CVD+target"
	d= data.frame(target=as.numeric(scale(x1)) ,
			conf[,c("group2","age","gender","BMI","smoke","drink","rc","disease_HTN","disease_T2D","disease_CVD")] );
	d$BMI= log(d$BMI)
	d$rc = log(d$rc)
	d=data.frame(na.omit(d))
	d$group2 = factor(d$group2,levels=c(g1,g2))
	# wilcox.test(d$target~d$group2)
	u = summary(glm(as.formula(formula), data=na.omit(d), family=binomial() ))
	bx1 = coef(u)[grepl("target",rownames(coef(u))),]

	c( p.val, fc , bx[1],bx[4],bx1[1],bx1[4],nrow(d) )

  }


library(compositions)

group_stat = function(df, grx2,sampletype="fecal"){
	#### extract body site sample
	nm =intersect( grx2_org[which(grx2_org$sampletype==sampletype),"ID"],colnames(df))
	print(length(nm))
	df_subset =df_org[,match(nm,colnames(df_org))]
	grx1 = grx2[match(nm,grx2$ID),]

	mean1 = apply(df_subset,1,base::mean)
	occ1 = apply(df_subset,1,function(x){base::mean(x>0)})

	### keep >10 samples species
	df_subset = df_subset[apply(df_subset,1,function(x){sum(x>0)})>0,]
	dim(df_subset)

	### CLR transformation to account composition effect
	df_subset <- t(clr(t(df_subset) + min(df_subset[df_subset > 0]) * 0.5 )) 
	df_subset = as.data.frame(df_subset) 


	

	### keep species consist with the RFE model, which filter species with less 10% occurence withn balanced 1:1 using propensity score matching (PSM) train set
		if(sampletype=="Gtissue"){	ff = "C:/Users/jiezhuye/Downloads/data/braken_data/20260719.rfe/20260719.rfe/20260718.tissue_select_filter_200.species_importance.xlsx"   }
		if(sampletype=="Gjuice"){	ff = "C:/Users/jiezhuye/Downloads/data/braken_data/20260719.rfe/20260719.rfe/20260718.juice_kraken_filter_200.species_importance.xlsx"  }
		if(sampletype=="fecal"){	ff = "C:/Users/jiezhuye/Downloads/data/braken_data/20260719.rfe/20260719.rfe/20260719.stool_kraken_filter_200.species_importance.xlsx" }
		if(sampletype=="tongue"){	ff = "C:/Users/jiezhuye/Downloads/data/braken_data/20260719.rfe/20260719.rfe/20260719.oral_kraken_filter_200.species_importance.xlsx"   }



		library(openxlsx)
		imp = read.xlsx(ff,sheet="Sheet 1")
		imp[,1]=gsub("___kraken","",imp[,1])

	### filter species with less 10% occurence in all samples in the body site
		imp=imp[imp[,1]%in%rownames(occ)[occ[,sampletype]>0.1],]
	### for tissue, just keep the species that pass multiple classfication.
		if(sampletype=="Gtissue"){
			tissue_sp = read.csv("C:\\Users\\jiezhuye\\Downloads\\data\\braken_data\\braken.tissue.species.list.csv",as.is=T)
			tissue_sp =tissue_sp[,1]
			tissue_sp=tissue_sp[tissue_sp%in%rownames(occ)[occ[,"Gtissue"]>0.1]]
			tissue_sp=tissue_sp[!tissue_sp%in%c("Pseudomonas_aeruginosa")]	
			imp=imp[imp[,1]%in%tissue_sp,]				
		}


	df_subset =df_subset[rownames(df_subset)%in%imp$var,]

		tmp2 =c()
		for(g2 in c("Superficial gastritis",
				"Atrophic gastritis","Intestinal metaplasia",
				"Intraepithelial Neoplasia","Gastric cancer I-II","Gastric cancer III-IV")){
				df_subset1 = df_subset[,which(grx1$group2%in%c("Healthy stomach", g2))]
				conf = grx1[which(grx1$group2%in%c("Healthy stomach", g2)),]
				conf = conf[match(colnames(df_subset1),conf$ID),]
				tmp1 = t(apply(df_subset1,1,function(x){ stat1(x ,conf, "Healthy stomach",g2 ) }))
				colnames(tmp1) = c("pvalue","fc","estimate","glmpval","estimate_more","glmpval_more","glm_n")
				tmp1=data.frame(tmp1)
				tmp1$group2=g2
				tmp1$gene = rownames(df_subset1)
				tmp2=rbind(tmp2,tmp1)
		}

				df_org_subset <- df_org[,colnames(df_org)%in%grx2_org[which(grx2_org$sampletype==sampletype),"ID"]]
				df_org_subset=df_org_subset[apply(df_org_subset,1,function(x){sum(x>0)})>0,]

				df_org_subset <- t(clr(t(df_org_subset) + min(df_org_subset[df_org_subset > 0]) * 0.5 )) 
				df_org_subset = as.data.frame(df_org_subset) 
				df_org_subset[1:4,1:4]

				nm1 =intersect(colnames(df_org_subset), grx2_org[which(grx2_org$sampletype==sampletype),"ID"])
				if(sampletype=="Gjuice"){ ## juice keep sample >10w reads
					nm1 = nm1[!nm1%in%grx2_org[which(grx2_org[,"rc"] <= 100000),"ID"]]
				}
				
					df_subset1 = df_org_subset[rownames(df_org_subset)%in%imp$var,match(nm1,colnames(df_org_subset))]

				conf = grx2_org[match(nm1,grx2_org$ID),]
				conf$group2 = conf$group1
				 table(conf$group2,conf$cite, useNA="always")
				dim(df_subset1)
				tmp1 = t(apply(df_subset1,1,function(x){ stat1(x ,conf, "con","gc" ) }))
				colnames(tmp1) = c("pvalue","fc","estimate","glmpval","estimate_more","glmpval_more","glm_n")
				tmp1=data.frame(tmp1)
				tmp1$gene = rownames(df_subset1)
				tmp1$group2="GC vs nonGC"



				tmp2=rbind(tmp2,tmp1)
				
	
	tmp2$source =sampletype
	tmp2$qvalue = p.adjust(tmp2$pvalue,method="BH")
	tmp2$glmqval = p.adjust(tmp2$glmpval,method="BH")
	tmp2$glmqval_more = p.adjust(tmp2$glmpval_more,method="BH")
	tmp2$pat =(abs(tmp2$qvalue)<0.05) * sign(tmp2$fc) 
	tmp2$glmpat =(abs(tmp2$glmqval)<0.05) * sign(tmp2$estimate) 
	tmp2$glmpat_more =(abs(tmp2$glmqval_more)<0.05) * sign(tmp2$estimate_more) 
		tmp2$mean1 = mean1[pmatch(tmp2$gene,names(mean1),duplicates.ok=T)]
		tmp2$occ1 = occ1[pmatch(tmp2$gene,names(occ1),duplicates.ok=T)]

	tmp2

}

library("tidyverse")
library("coin")
library("pROC")
library("yaml")
log.n0 =1e-5


		checklist =c("s__Fusobacterium_nucleatum","s__Streptococcus_anginosus",
				"s__Lancefieldella_parvula","s__Gemella_morbillorum",
				"s__Dialister_pneumosintes","s__Pseudomonas_aeruginosa",
				"s__Helicobacter_pylori","s__Escherichia_coli",
				"s__Neisseria_subflava","s__Porphyromonas_pasteri",
				"s__Streptococcus_mitis","s__Treponema_vincentii",
				"s__Prevotella_buccae")

set.seed(0)
pat = data.frame( rbind(group_stat(df, grx2,sampletype="fecal"),
			group_stat(df, grx2,sampletype="Gjuice"),
			group_stat(df, grx2,sampletype="Gtissue"),
			group_stat(df, grx2,sampletype="tongue")))



          pat$label = ifelse(pat$glmqval_more > 0.001, "*", NA)
           pat$label = ifelse(pat$glmqval_more < 0.001, "**", pat$label)
           pat$label = ifelse(pat$glmqval_more > 0.05, NA, pat$label)

pat$label1 <- NA

# Apply directional labels based on significance and fold change
pat$label1[pat$glmqval_more <= 0.01 & pat$fc > 0] <- "++"  # Highly significant positive
pat$label1[pat$glmqval_more <= 0.01 & pat$fc < 0] <- "--"  # Highly significant negative
pat$label1[pat$glmqval_more > 0.01 & pat$glmqval_more < 0.05 & pat$fc > 0] <- "+"  # Significant positive
pat$label1[pat$glmqval_more > 0.01 & pat$glmqval_more < 0.05 & pat$fc < 0] <- "-"  # Significant negative




replacex =function(x,org, value){
	for(ii in 1:length(org)){
		x[which(x==org[ii])]=value[ii]
	}
	x
}
table(pat$source)
pat$source = replacex(pat$source,c("fecal","Gjuice" , "Gtissue", "tongue"),
	c("Stool","Gastric fluid","Gastric tissue","Tongue coating") )
pat$source = factor(pat$source,levels=c("Tongue coating","Gastric fluid","Gastric tissue","Stool"))

table(pat$group2)
	pat$group2 =factor(pat$group2, levels = c("Superficial gastritis",
				"Atrophic gastritis","Intestinal metaplasia",
				"Intraepithelial Neoplasia","Precancer vs HS_SG" ,"Gastric cancer I-II","Gastric cancer III-IV","GC vs nonGC"))

write.csv(pat,"C:/Users/jiezhuye/Downloads/data/braken_data/Figure2d.pattern.kraken.csv")

		
