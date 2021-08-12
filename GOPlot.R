#This is a first try at GOplot an R package that can be used to visualize genes and enriched categories.
#First try only is only using two files EC$david.csv and EC$genelist.csv.  The vignette calls for more files but it doesn't appear to use them
#genelist is missing some columns


## setup root directory path
# tom epa windows
if(Sys.info()[4]=="DZ2626UTPURUCKE"){
  goplot_root <- file.path("c:", "git", "goplot")
}else
  # tom windows 
  if(Sys.info()[4]=="LZ2626UTPURUCKE"){
    goplot_root <- file.path("c:","git","goplot")
  } else
    # julia epa windows
    if(Sys.info()[4]=="LZ26JSTELMAN"){
      goplot_root <- file.path("c:","Users","jstelman","Git","goplot")
    } else {
      # default path: assumes your wd is set to this file's location
      goplot_root <- ""
    }

## packages to load
library(tidyverse)
library(GOplot)

## objects to load: david3, genelist

## load david3
# david output based on the top 1000 genes in genelist with all unmappable genes removed
david3 <- read.csv(file = file.path(goplot_root, "data_in", "david3.csv"), 
                   header = T, stringsAsFactors = T)

## load genelist
# gene list of top 1000 genes with unmappable genes removed
genelist <- read.csv(file = file.path(goplot_root, "data_in", "genelist2.csv"), 
                     header = T, stringsAsFactors = T)


# dim(david3)
# colnames(david3)
# View(david3)
# length(which(david3$adj_pval>0.05))
# length(david3$adj_pval)

#filter david output for adj_value < 0.05


# dim(genelist)
# colnames(genelist)
# View(genelist)


# what is with those three last columns? Get rid of them!
genelist[,4:6]<- NULL

#' #' the david data has different kinds of terms. some have ID conventions others don't
#' #' lets make it easier to make use of those that do by separating them out
#' #' anything in the ..._DIRECT category has a `GO:`ID`~`Term convention
#' david3direct <- david3 %>% filter(str_ends(Category, "_DIRECT")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = "~")
#' #' the _PATHWAY categories use the ID`:`Term convention
#' david3pathway <- david3 %>% filter(str_ends(Category, "_PATHWAY")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = ":")
#' #' INTERPRO - ID`:`Term
#' david3interpro <- david3 %>% filter(str_detect(Category, "^INTERPRO$")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = ":", extra = "merge")
#' #' BIOCARTA - ID`Pathway:`Term
#' david3biocarta <- david3 %>% filter(str_detect(Category, "^BIOCARTA$")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = "Pathway:")
#' #' OMIN_DISEASE - ID`~`Term
#' david3disease <- david3 %>% filter(str_detect(Category, "^OMIN_DISEASE$")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = "~")
#' #' SMART - ID`:`Term
#' david3smart <- david3 %>% filter(str_detect(Category, "^SMART$")) %>% 
#'   separate(col = Term, into = c("ID", "Term"), sep = ":")
#' # some of them have no id

david3new <- (function(){
  david3_w_idx <- rownames_to_column(david3)
  return(rbind(david3_w_idx %>% filter(str_ends(Category, "_DIRECT") |
                                         str_detect(Category, "^OMIN_DISEASE$")) %>% 
                 separate(col = Term, into = c("ID", "Term"), sep = "~"),
               david3_w_idx %>% filter(str_ends(Category, "_PATHWAY") |
                                         str_detect(Category, "^INTERPRO$") |
                                         str_detect(Category, "^BIOCARTA$") |
                                         str_detect(Category, "^SMART$")) %>% 
                 separate(col = Term, into = c("ID", "Term"), sep = "(Pathway)?:", extra = "merge")) %>% 
           right_join(y = david3_w_idx, 
                      by = c("rowname", "Category", "Genes", "adj_pval"), 
                      suffix = c("","_")) %>% 
           mutate(Term = ifelse(is.na(Term),Term_,Term)) %>% 
           select(-rowname, -Term_))})()

# Generate the plotting object
circ2 <- circle_dat(david3new, genelist)

# circ4
circ4 <- circ2 %>% filter(category %in% c("GOTERM_MF_DIRECT", "GOTERM_BP_DIRECT", "KEGG_PATHWAY"))


# Generate the plotting object
# circ2 <- circle_dat(david3, genelist)
# dim(circ2)
# View(circ2)
# colnames(circ2)
# # colnames(circ)
# summary(david3)
# #david3$adj_pval <- david3$adj_pval + 0.000001
# # View(circ)
# circ2$ID <- circ2$term
# View(circ2)

#Bar plot
# GOBar(circ2) 
# GOBar(subset(circ2, category == 'UP_KEYWORDS'))
# 
# drop_these <- which(circ2$category=="COG_ONTOLOGY")
# circ3 <- circ2[-drop_these,]
# keep_these <- which(circ2$category == "GOTERM_MF_DIRECT" |
#         circ3$category == "GOTERM_BP_DIRECT" |
#       circ3$category == "KEGG_PATHWAY")
# circ4 <- circ2[keep_these,]


# unique(circ4$category)

cat_list <- unique(circ4$category)
for(category_item in cat_list){
  print(category_item)
  print(GOBar(subset(circ4, category == category_item), title = category_item))
}

# GOBar(circ4, labels=3)

GOCircle(circ4)

View(circ4)

#Bubble plot
GOBubble(circ4) 
#unique(circ2$)



# circ <- circle_dat(EC$david, EC$genelist)
# View(circ)
# summary(circ)
# summary(circ2)
# GOBar(subset(circ, category == 'BP'))
# 
# circ2$adj_pval <- circ2$adj_pval + 0.00000001
# 
# 
# data(EC)
# EC$genes
# EC$genelist
# View(EC$genelist)
