#This is a first try at GOplot an R package that can be used to visualize genes and enriched categories.
#First try only is only using two files EC$david.csv and EC$genelist.csv.  The vignette calls for more files but it doesn't appear to use them
#genelist is missing some columns

library(GOplot)
library(tidyverse)
library(ggplot2)

#david output with based on the top 1000 genes in genelist with all unmappable genes removed
david3 <- read.csv("david3.csv", header = T, stringsAsFactors = T)
dim(david3)
colnames(david3)
View(david3)
length(which(david3$adj_pval>0.05))
length(david3$adj_pval)

#filter david output for adj_value < 0.05


#gene list of top 1000 genes with unmappable genes removed
genelist <- read.csv("genelist2.csv", header = T, stringsAsFactors = T)
dim(genelist)
colnames(genelist)
View(genelist)

#plotting object
circ2 <- circle_dat(david3, genelist)
dim(circ2)
View(circ2)
colnames(circ2)
colnames(circ)
summary(david3)
#david3$adj_pval <- david3$adj_pval + 0.000001
View(circ)
circ2$ID <- circ2$term
View(circ2)

#Bar plot
GOBar(circ2) 
GOBar(subset(circ2, category == 'UP_KEYWORDS'))

drop_these <- which(circ2$category=="COG_ONTOLOGY")
circ3 <- circ2[-drop_these,]
keep_these <- which(circ3$category == "GOTERM_MF_DIRECT" |
        circ3$category == "GOTERM_BP_DIRECT" |
      circ3$category == "KEGG_PATHWAY")
circ4 <- circ3[keep_these,]


unique(circ4$category)

cat_list <- unique(circ4$category)
for(category_item in cat_list){
  print(category_item)
  print(GOBar(subset(circ3, category == category_item)))
}

GOBar(circ4, labels=3)

GOCircle(circ4)

View(circ4)

#Bubble plot
GOBubble(circ4) 
#unique(circ2$)



circ <- circle_dat(EC$david, EC$genelist)
View(circ)
summary(circ)
summary(circ2)
GOBar(subset(circ, category == 'BP'))

circ2$adj_pval <- circ2$adj_pval + 0.00000001


data(EC)
EC$genes
EC$genelist
View(EC$genelist)
