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
library(pheatmap)

# save the detlist filepaths
detlists <-list.files(path = "C:/Users/jstelman/Git/ww2dw/DET_lists",full.names = T)

# change the names of th

for (x in detlists[1:1]){
  detname <- str_sub(x, 43, -5)
  dlist <- read.table(file=x, sep = "\t", quote = "\"", header = T) %>% 
    dplyr::filter(sign(des.log2FoldChange) == sign(edg.logFC)) %>% 
    transmute(transcriptID, avg.rank, direction = sign(edg.logFC), gene, hs_gene, name) %>%
    drop_na() %>% # do the uniquing before moving on, (take the row w/ least avg.rank)
    arrange(avg.rank)
  uplist <- dlist %>% dplyr::filter(direction == 1)
  downlist <- dlist %>% dplyr::filter(direction == -1)
  multi_gp2 <- gost(query = list(upregulated = unique(uplist$hs_gene), # now u don't need unique
                                 downregulated = unique(downlist$hs_gene)),
                    ordered_query = TRUE, 
                    multi_query = FALSE, 
                    user_threshold = 0.05, 
                    sources = c("KEGG"))
  if (!is.null(multi_gp2)){
    # make heatmap table
    fullhmtab <- multi_gp2$result %>% select(term_id, term_name, query, p_value) %>% 
      mutate(hmval. = -log10(p_value)*(ifelse(query=="upregulated",1,-1))) %>% 
      select(-query, -p_value) 
  }
}

for (x in detlists[-1]){
  detname <- str_sub(x, 43, -5)
  dlist <- read.table(file=x, sep = "\t", quote = "\"", header = T) %>% 
    dplyr::filter(sign(des.log2FoldChange) == sign(edg.logFC)) %>% 
    transmute(transcriptID, avg.rank, direction = sign(edg.logFC), gene, hs_gene, name) %>%
    drop_na() %>%
    arrange(avg.rank)
 uplist <- dlist %>% dplyr::filter(direction == 1)
 downlist <- dlist %>% dplyr::filter(direction == -1)
 multi_gp2 <- gost(query = list(upregulated = unique(uplist$hs_gene),
                                downregulated = unique(downlist$hs_gene)),
                   ordered_query = TRUE, 
                   multi_query = FALSE, 
                   user_threshold = 0.05, 
                   sources = c("KEGG"), # just the 1
 )
 if (!is.null(multi_gp2)){
   # make heatmap table
   hmtabx <- multi_gp2$result %>% select(term_id, term_name, query, p_value) %>% 
     mutate(hmval. = -log10(p_value)*(ifelse(query=="upregulated",1,-1))) %>% 
     select(-query, -p_value)
   fullhmtab <- full_join(x = fullhmtab, y = hmtabx, by = c("term_id", "term_name"),
                          suffix = c("",detname))
 }
}


# fix that first column's name
colnames(fullhmtab)[3] <- "hmval.ctrlvefp_field"


# replace nas with 0s
fullhmtab <- fullhmtab %>% mutate_all(replace_na,0)

# make the table with the kegg ids as rownames and make heatmap
fullhmtab2 <- fullhmtab %>% group_by(term_id, term_name) %>% summarize_all(mean) # this can go
fullhmtab2 <- column_to_rownames(fullhmtab2, var = "term_id")
png(filename = file.path(goplot_root, "plots", "keggheatmapIDs.png"),
    width = 670, height = 1940, units = "px", bg = "white")
    # width = 560, height = 1615, units = "px", bg = "white")
pheatmap(fullhmtab2[,-1],
         # cluster_rows=F,
         # cluster_cols = F,
         main = "Kegg pathways",
         cex = 1.3)
dev.off()

# make the table with term names as rownames and make heatmap
fullhmtab3 <- fullhmtab %>% group_by(term_name, term_id) %>% summarize_all(mean)
fullhmtab3 <- column_to_rownames(fullhmtab3, var = "term_name")
png(filename = file.path(goplot_root, "plots", "keggheatmapDETAILs.png"),
    width = 800, height = 1940, units = "px", bg = "white")
    # width = 560, height = 1615, units = "px", bg = "white")
pheatmap(fullhmtab3[,-1],
         # cluster_rows=F,
         # cluster_cols = F,
         main = "Kegg pathways details",
         cex = 1.3)
dev.off()




# there are a lot of genes that map to both up- and down- regulated transcripts
# there are a lot of genes that map to more than one transcript
# option 1: ignore this and push forward
# option 2: redo the workflow, at the gene level from the start
# option 3: redo the workflow, at the hs_gene level from the start

# ----------------------------
# assuming option 1:
#' now we have to decide when and how to aggregate
#' how: we can take 
#'   (i) avg-rank 
#'   (ii) min-rank
#'   (iii) unique-gene, and don't indicate rank ordering in the lists
#' when: we can do this 
#'   (i) before splitting into up/down which will eliminate double agent genes
#'   (ii) within up/down each separately which will not eliminate double agents
