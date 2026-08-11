result=read.csv("~/Desktop/...../Isolate_20260720.csv",as.is=T,head=T,row.names=1)

mpa_tree <- ape::read.tree("~/Desktop/...../20260716.isolate_PlusPa.unrooted.tree")

library(ggtree)
library(treeio)
library(ggsci)
library("scales")
library(ggplot2)
library(ggtreeExtra)
library(dplyr)

factor_data <- data.frame(ID=mpa_tree$tip.label,
                          label = mpa_tree$tip.label, 
                          label1 = ifelse(result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"tissue_sp"]==1, 
                                          gsub("_"," ",gsub("s__","",result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"profid"])),""),
                          phylum  = result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"phylum"],
                          Novel = grepl("g__", result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"rename"]) &
                            is.na( result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"profid"]),
                          name= result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"newtax1"],
                          num = result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"culture_num"],
                          mp4 = result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"significant_rows"],
                          mags = result[match(mpa_tree$tip.label,gsub("isolate_new---","",result[,1]) ),"gastric.tissue"]>0 )

keep_labels <- factor_data %>%
  group_by(label1) %>%
  mutate(keep = ifelse(label1 == "", FALSE, num == max(num))) %>%
  ungroup()

# Keep only label1 entries with the largest num, others set to ""
factor_data$label1 <- ifelse(keep_labels$keep, factor_data$label1, "")

phylum_colors <- c("p__Actinomycetota" = "#F2A099",
                   "p__Bacillota" = "#F6B366",
                   "p__Bacteroidota" ="#6B8FAD", 
                   "p__Fusobacteriota" = "#5BA8A2",
                   "p__Pseudomonadota" = "#70A5D9")

mp4_colors <- c("FALSE" = "white",
                "TRUE" = "black")

mp4_vector <- c(
  "FALSE" = 0,   # square
  "TRUE" = 19   # circle
)

nhx_reduced <- full_join(mpa_tree, factor_data, by = 'label')

clade1 <- MRCA(mpa_tree, mpa_tree$tip.label[mpa_tree$tip.label%in%factor_data[factor_data$phylum=="p__Actinomycetota","label"]])
clade2 <- MRCA(mpa_tree, mpa_tree$tip.label[mpa_tree$tip.label%in%factor_data[factor_data$phylum=="p__Bacillota","label"]])
clade3 <- MRCA(mpa_tree, mpa_tree$tip.label[mpa_tree$tip.label%in%factor_data[factor_data$phylum=="p__Bacteroidota","label"]])
clade4 <- MRCA(mpa_tree, mpa_tree$tip.label[mpa_tree$tip.label%in%factor_data[factor_data$phylum=="p__Fusobacteriota","label"]])
clade5 <- MRCA(mpa_tree, mpa_tree$tip.label[mpa_tree$tip.label%in%factor_data[factor_data$phylum=="p__Pseudomonadota","label"]])

p1= ggtree(nhx_reduced, layout="fan", size=0.01, open.angle=5, branch.length="none")+
  ggtitle("culture")+
  geom_tree(size=0.01)  +
  geom_hilight(node = clade1,  fill = "#F2A099", extend = 0.0017) +
  geom_hilight(node = clade2,  fill = "#F6B366", extend = 0.0017) +
  geom_hilight(node = clade3,  fill = "#6B8FAD", extend = 0.0017) +
  geom_hilight(node = clade4,  fill = "#5BA8A2", extend = 0.0017) +
  geom_hilight(node = clade5,  fill = "#70A5D9", extend = 0.0017) +
  geom_tippoint(aes(shape = mp4, color = mp4), size=1) +
  geom_tippoint(aes(fill = mags), shape = 21, size=1, color="white", stroke=0.5, 
                position = position_nudge(x = 2.5)) +
  geom_tippoint(aes(fill = Novel), shape = 21, size=1, color="white", stroke=0.5, 
                position = position_nudge(x = 4.5)) +
  scale_color_manual(values= mp4_colors )+
  scale_fill_manual(values= mp4_colors, na.value = "white" )+
  scale_shape_manual(values = mp4_vector )+
  geom_fruit(
    geom = geom_bar,
    mapping=aes(y=label, x=num),
    orientation = "y",
    stat="identity",
    offset = 0.15,
    pwidth = 0.5
  )  +
  geom_tiplab(aes(label = label1), linesize=0, offset = 12.5, size=2) 

ggsave("~/Desktop/20260108.gtdb_tree.unrooted.tree.pdf",p1,width=7,height=7)

factor_data2 <- factor_data %>%
  select(label, label1, phylum, Novel, name, num, mp4, mags)

p0 <- ggtree(mpa_tree, layout="fan", open.angle=5, branch.length="none", linewidth=0.2) %<+% factor_data2

phylum_fill <- c(
  "p__Actinomycetota"  = "#F2A099",
  "p__Bacillota"       = "#F6B366",
  "p__Bacteroidota"    = "#6B8FAD",
  "p__Fusobacteriota"  = "#5BA8A2",
  "p__Pseudomonadota"  = "#70A5D9"
)

p1 <- p0 +
  ggtitle("culture") +
  
  # Highlight phyla with transparency
  geom_hilight(node=clade1, fill=phylum_fill["p__Actinomycetota"], alpha=0.22, extend=0.0012) +
  geom_hilight(node=clade2, fill=phylum_fill["p__Bacillota"],     alpha=0.22, extend=0.0012) +
  geom_hilight(node=clade3, fill=phylum_fill["p__Bacteroidota"],  alpha=0.22, extend=0.0012) +
  geom_hilight(node=clade4, fill=phylum_fill["p__Fusobacteriota"],alpha=0.22, extend=0.0012) +
  geom_hilight(node=clade5, fill=phylum_fill["p__Pseudomonadota"],alpha=0.22, extend=0.0012) +
  
  # Phylum legend (invisible points, only for legend display)
  geom_tippoint(aes(fill = phylum), shape = 22, size = 3, alpha = 0, show.legend = TRUE) +
  scale_fill_manual(values = phylum_fill, name = "Phylum") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 22, size = 4, color = NA))) +
  
  # Important: reset fill scale to avoid conflict with mags/Novel
  ggnewscale::new_scale_fill() +
  
  # Ring 1: mp4 (shape + grey/black)
  geom_tippoint(aes(shape=mp4, color=mp4), size=1.2) +
  scale_shape_manual(values=c("FALSE"=1, "TRUE"=16)) +
  scale_color_manual(values=c("FALSE"="grey80", "TRUE"="black"), name="mp4") +
  
  # Ring 2: mags (independent fill, blue)
  ggnewscale::new_scale_fill() +
  geom_tippoint(aes(fill=mags), shape=21, size=1.2, color="white", stroke=0.4,
                position=position_nudge(x=2.2)) +
  scale_fill_manual(values=c("FALSE"="white", "TRUE"="#2C7FB8"), na.value="white", name="mags") +
  
  # Phylum color ring (close to tree)
  ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = phylum),
    offset = 0.01,
    pwidth = 0.02
  ) +
  scale_fill_manual(values = phylum_fill, name = "Phylum")+
  
  # Ring 3: Novel (independent fill, red)
  ggnewscale::new_scale_fill() +
  geom_tippoint(aes(fill=Novel), shape=21, size=1.2, color="white", stroke=0.4,
                position=position_nudge(x=3.2)) +
  scale_fill_manual(values=c("FALSE"="white", "TRUE"="#D7301F"), na.value="white", name="Novel") +
  
  # Outer ring bar chart (subdued)
  geom_fruit(
    geom=geom_bar,
    mapping=aes(x=num),
    orientation="y",
    stat="identity",
    offset=0.10,
    pwidth=0.38,
    alpha=0.75,
    fill="grey30"
  ) +
  
  # Labels (cleaner alignment)
  geom_tiplab2(aes(label=label1), align=TRUE, linesize=0,
               offset=6.5, size=3, fontface="italic") +
  
  # Theme and legend
  theme_void() +
  theme(
    plot.title = element_text(hjust=0.5, size=11, face="bold"),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text  = element_text(size=7)
  )

ggsave("~/Desktop/20260715.gtdb_tree.unrooted.tree.beautified.pdf",
       p1, width=8, height=8, device="pdf")