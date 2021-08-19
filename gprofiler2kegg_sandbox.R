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
#! this will be fixed once this is integrated into the ww2dw repo
detlists <-list.files(path = "C:/Users/jstelman/Git/ww2dw/DET_lists",full.names = T)


for (x in detlists[1:1]){
  # extract the name substring from the file path string
  detname <- str_sub(x, 43, -5)
  #! ^ use regex to change names to what Tom suggested
  
  # read in the table from the file and change stuff in it using dplyr
  #! name it something better than xdata
  xdata <- read.table(file=x, sep = "\t", quote = "\"", header = T) %>% 
    # delete all cases when edgeR and DESeq didn't agree on direction of regulation
    dplyr::filter(sign(des.log2FoldChange) == sign(edg.logFC)) %>% 
    # keep only 6 columns, 3 of which are ever used again, 3 of which I'm just hoarding
    #   direction tells direction of regulation
    transmute(transcriptID, avg.rank, direction = sign(edg.logFC), gene, hs_gene, name) %>%
    # drops rows that have NA in the hs_gene column
    #! do the uniquing before moving on, (take the row w/ least avg.rank)
    drop_na() %>%
    # order the rows by most to least significant
    arrange(avg.rank)
  # make list of hs_genes, separate into upregulated and downregulated, keeping them ordered
  # in each one, keep only the first occurance of each gene (using unique())
  uplist <- xdata %>% dplyr::filter(direction == 1) %>% .$hs_gene %>% unique()
  downlist <- xdata %>% dplyr::filter(direction == -1) %>% .$hs_gene %>% unique()
  # get the multiple gprofiler2 result
  multi_gp2 <- gost(query = list(upregulated = uplist, 
                                 downregulated = downlist),
                    ordered_query = TRUE,   # assign meaning to the ordering of genes
                    multi_query = FALSE,   # kind of like saying "use rbind, not cbind"
                    user_threshold = 0.05, 
                    sources = c("KEGG"))
  if (!is.null(multi_gp2)){
    # make heatmap table
    # The 1st of up to 15 tables to be full_join()ed, so convention differs slightly from below
    fullhmtab <- multi_gp2$result %>% select(term_id, term_name, query, p_value) %>% 
      # assign heatmap value column and calculate its elements
      mutate(hmval. = -log10(p_value)*(ifelse(query=="upregulated",1,-1))) %>% 
      # take out the trash
      select(-query, -p_value) 
  }
}

for (x in detlists[-1]){
  # extract the name substring from the file path string
  detname <- str_sub(x, 43, -5)
  #! ^ use regex to change names to what Tom suggested
  
  # read in the table from the file and change stuff in it using dplyr
  #! name it something better than xdata
  xdata <- read.table(file=x, sep = "\t", quote = "\"", header = T) %>% 
    # delete all cases when edgeR and DESeq didn't agree on direction of regulation
    dplyr::filter(sign(des.log2FoldChange) == sign(edg.logFC)) %>% 
    # keep only 6 columns, 3 of which are ever used again, 3 of which I'm just hoarding
    #   direction tells direction of regulation
    transmute(transcriptID, avg.rank, direction = sign(edg.logFC), gene, hs_gene, name) %>%
    # drops rows that have NA in the hs_gene column
    #! do the uniquing before moving on, (take the row w/ least avg.rank)
    drop_na() %>%
    # order the rows by most to least significant
    arrange(avg.rank)
  # make list of hs_genes, separate into upregulated and downregulated, keeping them ordered
  # in each one, keep only the first occurance of each gene (using unique())
  uplist <- xdata %>% dplyr::filter(direction == 1) %>% .$hs_gene %>% unique()
  downlist <- xdata %>% dplyr::filter(direction == -1) %>% .$hs_gene %>% unique()
  # get the multiple gprofiler2 result
  multi_gp2 <- gost(query = list(upregulated = uplist, 
                                 downregulated = downlist),
                    ordered_query = TRUE,   # assign meaning to the ordering of genes
                    multi_query = FALSE,   # kind of like saying "use rbind, not cbind"
                   user_threshold = 0.05, 
                   sources = c("KEGG"), # just the 1
 )
 if (!is.null(multi_gp2)){
   # make this rotation's heatmap table
   hmtabx <- multi_gp2$result %>% select(term_id, term_name, query, p_value) %>% 
     # assign heatmap value column and calculate its elements (signed -log10 adj pval)
     mutate(hmval. = -log10(p_value)*(ifelse(query=="upregulated",1,-1))) %>% 
     # take out the trash
     select(-query, -p_value)
   # join it to the full heatmap data table
   fullhmtab <- full_join(x = fullhmtab, y = hmtabx, by = c("term_id", "term_name"),
                          suffix = c("",detname))
 }
}

# currently, all numeric columns except 1st one are called hmval.[ctrl]v[treatment][cohort tag]
# so let's manually go back and fix that first numeric column's name so it matches the rest
colnames(fullhmtab)[3] <- "hmval.ctrlvefp_field"

# the mathematical integrity is questionable, but the only way clustering can work is if we
# replace nas with 0s
fullhmtab <- fullhmtab %>% mutate_all(replace_na,0)

# --------------------
# time to make 2 versions of the table: 1 with shorter rownames, 1 with more informative rownames

# make the table with the kegg IDs as rownames and make heatmap
fullhmtab2 <- fullhmtab %>% 
  # some things were listed as both up- and down- regulated
  # that means term aren't strictly unique, some are found twice
  # that won't fly - only one row per term! 
  # therefore, we must aggregate all double-sightings (using mean for now)
  group_by(term_id, term_name) %>% summarize_all(mean) %>%   #! this can go if aggregating moves
  # now that there are no repeats in the term columns, we can assign them to be the rownames
  column_to_rownames(var = "term_id")
# designate a png file for the heatmap
png(filename = file.path(goplot_root, "plots", "keggheatmapIDs.png"),
    width = 670, height = 1940, units = "px", bg = "white")
    # width = 560, height = 1615, units = "px", bg = "white")
# make the heatmap
pheatmap(fullhmtab2[,-1], # first column is term names, exclude it
         # cluster_rows = F, cluster_cols = F,  # if no coalescing of NAs, uncomment this
         main = "Kegg pathways",
         cex = 1.3)
dev.off()

# make the table with term NAMEs as rownames and make heatmap
fullhmtab3 <- fullhmtab %>% group_by(term_name, term_id) %>% summarize_all(mean) %>% 
  column_to_rownames(var = "term_name")
png(filename = file.path(goplot_root, "plots", "keggheatmapDETAILs.png"),
    width = 800, height = 1940, units = "px", bg = "white")
    # width = 560, height = 1615, units = "px", bg = "white")
pheatmap(fullhmtab3[,-1],  # first column is term IDs, exclude it
         # cluster_rows = F, cluster_cols = F,  # if no coalescing of NAs, uncomment this
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
