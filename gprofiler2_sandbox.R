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

# prepare pdf file
pdf(file = "plots/gprofilerPlots.pdf")
# save names so we know which ones didn't fail
plotnames <- list()
# make them and print them in pdf, also save them to global env
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
  assign(x = detname,
         value = gostplot(gostreslist, capped = TRUE, interactive = FALSE) + ggtitle(detname) +
           theme(plot.title = element_text(size = 30, hjust = 0.5, margin = margin(10,0,-10,0))),
         envir = globalenv())
  plotnames[detname] <- detname
  print(get(detname, envir = globalenv()))
}
#close pdf connection
dev.off()

# make the legend
gp2legend <- as.data.frame(cbind(
  values = NA, colors = c("GO:MF", "GO:BP", "KEGG", "REAC"))) %>% 
  mutate(colors = relevel(as.factor(colors), "GO:MF")) %>% 
  ggplot(aes(values, fill = colors)) + geom_bar() + 
  scale_fill_manual(values = c("#DC3912","#FF9900","#DD4477","#3366CC")) + 
  theme(axis.title = element_blank(), legend.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"), panel.ontop = T,
        legend.justification = "center", legend.position = c(0.5,0.5),
        legend.text = element_text(size = 20, margin = margin(6,6,6,6)),
        legend.key.size = unit(1, 'cm'))

# use cowplot's plot_grid to arrange all 14 successful plots and the legend in place of the 15th
gprofiler5x3plot <- plot_grid(ctrlvefp1, ctrlvefp2, ctrlvefp3, ctrlvefp4, ctrlvefp_field,
                              ctrlvup1, ctrlvup2, ctrlvup3, ctrlvup4, ctrlvup_field,
                              upvefp1, upvefp2, upvefp3, 
                              gp2legend,
                              # ggplot()+geom_blank()+theme_nothing(), 
                              upvefp_field,
                         nrow = 5, labels = LETTERS[1:15], byrow = FALSE,
                         label_y = 0.95, label_x = 0.95) +
  theme(plot.margin = margin(0,0,0,40)) +
  draw_label(x = -0.01, angle = 90, label = "-log10(p-adj)", size = 30,)


# export it to a file
png(filename = file.path(goplot_root, "plots", "5x3gprofiler.png"),
    width = 2240, height = 2200, units = "px", bg = "white")
print(gprofiler5x3plot)
dev.off()


links <- 
  lapply(list(1:5,6:10,11:15), function(x){
    detnames <- sapply(detlists[x], function(a) str_sub(a, 43, -5))
    dlists <- lapply(detlists[x], function(b){
      dl <- read.table(file=b, sep = "\t", quote = "\"", header = T)
      return(unique(drop_na(dl)$hs_gene))
      })
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
