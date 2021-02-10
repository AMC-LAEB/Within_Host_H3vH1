### R script using ggtree to create phylogenetic tree figures  

# load libraries 
require(ggtree)
library(cowplot)
library(ggplot2)
library(ggsci)
library("scales")
library(RColorBrewer)
library(scales)
library("ggthemes")
library("LaCroixColoR")
library(dplyr)
library(sommer)
library(stringr)

## inputs
# set working directory accordingly 
setwd("/path/to/Within_Host_H3vH1/H3N2/reference/")
# set output filename
outfname <- "H3N2_SFig_concatenated_wgs.pdf"

# newick tree fname (must be in working directory)
tree<-read.tree("ggtree_concatenated_wgs.nwk")

# read meta data file 
metadf<-read.csv("ggtree_concatenated_wgs.meta.csv", stringsAsFactors=FALSE)
metadf[is.na(metadf)] = ""
row.names(metadf)=metadf$index 
metadf$show_name = as.numeric(metadf$show_name) # show_name is a binary column whether to show print_name
metadf$show_tip = as.numeric(metadf$show_tip) # show_tip is a binary column whether to show tip shape
metadf$subject[metadf$subject_id==""] <- NA

# generate distinct colours 
n <- 53
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
all_colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = sample(all_colours, n)

p1 <- ggtree(tree, size=.7, ladderize=T, right=F)%<+%metadf 
p1 <- p1 + theme_tree() + xlim(0, 0.011) + # higher xlim max = more compact 
  aes(color=subject_id)  + 
  geom_tiplab(aes(color=subject_id, label=print_name, subset=(show_name==1)), size=1.5, offset=0., linesize = 0.2, align=T) + # color="black"
  scale_fill_manual(na.value="black", values = c(col_vector)) + 
  geom_tippoint(aes(color=subject_id, subset=(show_tip==1)), size=2, pch=I(19)) + # colour tips; change fill=<var> to coloumn name to be coloured
  geom_treescale(x=0, y=40, fontsize=3, linesize=0.5, offset=1)
    
p1
ggsave(outfname, width = 11.7, height = 16.5, units = "in", limitsize = FALSE)
