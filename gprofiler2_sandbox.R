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

# packages to load
library(gprofiler2)
library(tidyverse)
library(cowplot)

## load genelist
# gene list of top 1000 genes with unmappable genes removed
genelist <- read.csv(file = file.path(goplot_root, "data_in", "genelist2.csv"), 
                     header = T, stringsAsFactors = T)

## load david3
# david output based on the top 1000 genes in genelist with all unmappable genes removed
david3 <- read.csv(file = file.path(goplot_root, "data_in", "david3.csv"),
                   header = T, stringsAsFactors = F)

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
           mutate(Term = ifelse(!is.na(Term),Term,Term_)) %>%
           select(-rowname, -Term_)
         )
  })()

# save the detlist filepaths
detlists <-list.files(path = "C:/Users/jstelman/Git/ww2dw/DET_lists",full.names = T)


# this works
plots <- list()
# this mostly works :/
for (x in detlists){
  detname <- str_sub(x,43,-5)
  dlist <- read.table(file=x, sep = "\t", quote = "\"", header = T)
  gostreslist <- gost(query = unique(drop_na(dlist)$hs_gene),
                      organism = "hsapiens", 
                      ordered_query = FALSE, 
                      multi_query = FALSE, 
                      significant = TRUE, 
                      exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, 
                      evcodes = FALSE, 
                      user_threshold = 0.01, 
                      correction_method = "g_SCS", 
                      #domain_scope = "annotated", 
                      custom_bg = NULL, 
                      numeric_ns = "", 
                      sources = c("GO:BP","GO:MF","KEGG","REAC"), # just the 4
                      as_short_link = FALSE)
  plots[detname] <- gostplot(gostreslist, capped = FALSE, interactive = FALSE) + ggtitle(detname)
}

# this doesn't do what it should
# use cowplot's plot_grid to arrange all 5 plots and the legend
gprofiler3x5plot <- plot_grid(plotlist = plots, 
                         nrow = 3, labels = LETTERS[1:15],
                         label_y = 0.85, label_x = 0.85) 

# this is a sad result
# export it to a file
png(filename = file.path(goplot_root, "plots", "3x5gprofiler.png"),
    width = 1666, height = 1080, units = "px", bg = "white")
print(gprofiler3x5plot)
dev.off()


links <- 
  lapply(list(1:5,6:10,11:15), function(x){
    detnames <- sapply(detlists[x], function(a) str_sub(a, 43, -5))
    dlists <- lapply(detlists[x], function(b){
      dl <- read.table(file=b, sep = "\t", quote = "\"", header = T)
      return(unique(drop_na(dl)$hs_gene))
      })
    #dlist <- read.table(file=x, sep = "\t", quote = "\"", header = T)
    gostreslist <- gost(query = dlists,
                        organism = "hsapiens", 
                        ordered_query = FALSE, 
                        multi_query = TRUE, 
                        significant = TRUE, 
                        exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, 
                        evcodes = FALSE, 
                        user_threshold = 0.01, 
                        correction_method = "g_SCS", 
                        #domain_scope = "annotated", 
                        custom_bg = NULL, 
                        numeric_ns = "", 
                        sources = c("GO:BP","GO:MF","KEGG","REAC"), # just the 4
                        as_short_link = TRUE)
  return(gostreslist)
  })
