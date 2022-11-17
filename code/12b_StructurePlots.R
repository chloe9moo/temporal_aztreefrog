##TITLE: AZ treefrog: STRUCTURE output management + plot creation
##AUTHOR: C. E. Moore
##Updated on 5 Oct 2021

library(tidyverse); library(ggh4x) #; library(adegenet); library(paletteer); library(cowplot)

## home path
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"

## load data
# source(file = paste0(PATH, "/code/01c_DataLoad.R"))
pop.key <- read.csv(paste0(PATH, "/microsat_data/year_pop_number_key4structure.csv"))

aztf <- read.csv(paste0(PATH, "/microsat_data/final_aztf_loci_nosibs.csv"))
aztf <- aztf %>% mutate(year = ifelse(year == 2018, 2019, year)) %>%
  mutate(year = ifelse(year == 2015, 2014, year))

## delta K function
#modified from STRUCTURE HARVESTER python script (https://taylor0.biology.ucla.edu/structureHarvester/faq.html)
# 1/ average the L(K) over the x (say 20) replicates
# 2/ estimate from these averages L''(K) as abs( L(K+1) - 2L(K) + L(K-1) )
# 3/ divide by the standard deviation of L(K) (sd of the different replicates for the same K)
delta_K <- function (x) {
  if(!all(c("K", "Rep", "MeanLnP_K") %in% colnames(x))) {
    stop("Dataframe needs 'K', 'Rep', and 'MeanLnP_K' columns.")
  }
  
  dKc <- x %>% 
    group_by(K) %>%
    rename(LnP_K = MeanLnP_K) %>%
    mutate(MeanLnP_K = mean(LnP_K, na.rm=T),
           stdLnP_K = sd(LnP_K, na.rm=T),
           LnDP_K = as.numeric(NA), #Ln'(K)
           LnDDP_K = as.numeric(NA), #|Ln"(K)|
           deltaK = as.numeric(NA)) %>% #delta K
    select(-LnP_K, -Rep) %>%
    distinct() %>%
    arrange(K) %>% #make sure its arranged numerically by K
    ungroup() 
  
  #calc Ln'(K)
  for (i in 2:nrow(dKc)) {
    thisK <- dKc[i, ]$MeanLnP_K
    prevK <- dKc[i-1, ]$MeanLnP_K
    dKc[i, "LnDP_K"]$LnDP_K <- thisK - prevK
  }
  #calc Ln"(K)
  for (i in 2:nrow(dKc)) {
    dKc[i, ]$LnDDP_K <- abs(dKc[i+1, ]$MeanLnP_K - 
                                  2*dKc[i, ]$MeanLnP_K +
                                  dKc[i-1, ]$MeanLnP_K)
    dKc[i, ]$deltaK <- dKc[i, ]$LnDDP_K / dKc[i, ]$stdLnP_K
  }
  
  cat("Most likely K =", dKc %>% filter(!is.na(deltaK)) %>% slice_max(deltaK, n=1, with_ties = F) %>% pull(K),"\n")
  return(dKc)
}

set_cluster_df <- function(x, subset = 'all') {
  # set year to grab the correct keys for relabeling the inferred ancestry
  # if nothing is specified, it will load the full set (all years + pops)
  if (subset == 'all') {
    ind.key <- read.csv(paste0(PATH, "/microsat_data/ind_key4structure.csv"))
    cat("Using all years, all pops for key\n")
  } else {
    ind.key <- read.csv(paste0(PATH, "/microsat_data/ind_key4structure_", subset, ".csv"))
    cat("Using", subset, "for key\n")
  }
  
  #reformat dataframe from structure output and 12c_str_text2csv.sh
  cluster_df <- x[,-c(1, ncol(x))] %>%
    rename_with(.cols = everything(), ~paste0("clust", seq_along(.))) %>%
    mutate(Label = row_number()) %>%
    left_join(., ind.key) %>%
    left_join(., aztf %>% select(ID, year, pop)) %>%
    relocate(ID, year, pop)
  
  return(cluster_df)
}

structure_plot_it <- function(cluster.assignments, color.assignments, plot.title, n.yrs = FALSE, cluster.relevel = NULL) {

#cluster.assignments should be the inferred ancestry csv already gone through set_cluster_df
  if(!all(c("ID", "year", "pop", "clust1", "Label") %in% colnames(cluster.assignments))) {
    stop("Dataframe needs 'ID', 'year', 'pop', 'Label', and at least one 'clust' column.")
  }
  
#set up cluster assignment dataframe for ggplot
  pop.fact <- unique(cluster.assignments[, "pop"])
  pop.fact <- as.character(pop.fact[order(pop.fact)])
  
  yr.fact <- unique(cluster.assignments[, "year"])
  yr.fact <- as.character(yr.fact[order(yr.fact)])
  
  long.clust <- pivot_longer(cluster.assignments, cols = starts_with("clust"), names_to = "cluster", values_to = "prob")
  long.clust$pondfact <- fct_relevel(as.factor(long.clust$pop), pop.fact)
  long.clust$yearfact <- fct_relevel(as.factor(long.clust$year), yr.fact)
  if(is.null(cluster.relevel) == FALSE) {
    long.clust <- long.clust %>%
      mutate(cluster = fct_relevel(as.factor(cluster), cluster.relevel))
  } else { 
    long.clust <- long.clust %>% mutate(cluster = factor(cluster))
    }
  
#group.colors should be character setting cluster color assignments, must be ggplot readable colors
  if(length(color.assignments) != length(unique(long.clust$cluster))) {
    stop("Not enough or too many colors assigned. Need ", length(unique(long.clust$cluster)), " colors.")
  }
#order of colors set matters (will change which cluster it's assigned to)
  names(color.assignments) <- unique(long.clust$cluster)
  
#prep graphics
  facetstrips <- strip_nested(
    text_x = elem_list_text(size = c(12, 6)),
    by_layer_x = TRUE, clip = "off"
  )
  
  if(n.yrs == TRUE) { 
    graphic.list <- list(
      geom_col(size = 0.01, width = 1),
      facet_nested(~ yearfact+pondfact,
                   #~ pondfact+yearfact, #switch which is on top or bottom in the fig
                   switch = "x",
                   nest_line = element_line(size = 1, lineend = "round"),
                   scales = "free", space = "free", strip = facetstrips),
      #theme_minimal(),
      labs(x = "", y = "membership probability"),
      scale_y_continuous(expand = c(0, 0)),
      scale_x_discrete(expand = expansion(add = 0.5)),
      scale_fill_manual(values = group.colors),
      theme(
        axis.ticks      = element_blank(),
        panel.background  = element_rect(fill = "white", colour = NA),
        #panel.border      = element_blank(),
        strip.background  = element_blank(),
        #plot.background   = element_blank(),
        panel.spacing.x = unit(0.18, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
      ))
  } else {
    graphic.list <- list(
      geom_col(size = 0.01, width = 1),
      facet_nested(~ pondfact,
                   switch = "x",
                   nest_line = element_line(size = 1, lineend = "round"),
                   scales = "free", space = "free", strip = facetstrips),
      #theme_minimal(),
      labs(x = "", y = "membership probability"),
      scale_y_continuous(expand = c(0, 0)),
      scale_x_discrete(expand = expansion(add = 0.5)),
      scale_fill_manual(values = group.colors),
      theme(
        axis.ticks      = element_blank(),
        panel.background  = element_rect(fill = "white", colour = NA),
        #panel.border      = element_blank(),
        strip.background  = element_blank(),
        #plot.background   = element_blank(),
        panel.spacing.x = unit(0.18, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
      ))
  }
  
  #plot
  p <- ggplot(long.clust, aes(factor(ID), prob, fill = cluster)) +
    graphic.list +
    ggtitle( paste0( plot.title, "; K = ", length(unique(long.clust$cluster)) ) )
  
  return(p)
}

#All pops + all years ----
##+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/K24/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK_K24.csv"), row.names = F)

##+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK_K24.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize

#convert the output final from structure to readable here using 12c_str_text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/K24/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust)
write.csv(clust, file=paste0(PATH, "/structure_results/AllInd_dK2_clusterassignments.csv"), row.names = F)

##+ make ggplot ----
group.colors <- c("#12263A", "#F4D1AE")
p.all <- structure_plot_it(clust, group.colors, "All years, All pops", n.yrs = T)
ggsave(paste0(PATH, "/structure_results/figures/pop_on_year_strct.png"), plot = p.all, width = 10, height = 6)

# ggsave(paste0(PATH, "/structure_results/figures/year_on_pop_strct.png"), plot = last_plot(), width = 10, height = 6)

#2014 ----
##+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2015/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2015_maxK9.csv"), row.names = F)

##+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2015_maxK9.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize
# str_out %>% filter(K == 3) %>% slice_max(MeanLnP_K, n=1)

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2015/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2014')
write.csv(clust, file=paste0(PATH, "/structure_results/r2014_dK2_clusterassignments.csv"), row.names = F)

# clust <- read.delim(paste0(PATH, "/structure_results/r2015/inferred_ancestry_K3.csv"), header = F, sep = ",")
# clust <- set_cluster_df(clust, '2014')
# write.csv(clust, file=paste0(PATH, "/structure_results/r2014_dK3_clusterassignments.csv"), row.names = F)

##+ make ggplot ----
group.colors <- c("#12263A", "#F4D1AE")
structure_plot_it(clust, group.colors, "2014")
ggsave(paste0(PATH, "/structure_results/figures/strct_2014_K2.png"), plot = last_plot(), width = 10, height = 6)

clust <- read.csv(paste0(PATH, "/structure_results/r2014_dK3_clusterassignments.csv"))
group.colors <- c("#00C6B8", "#12263A", "#F4D1AE")
structure_plot_it(clust, group.colors, "2014")
ggsave(paste0(PATH, "/structure_results/figures/strct_2014_K3.png"), plot = last_plot(), width = 10, height = 6)

#// pops 1, 3, 4 ----
##//+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2015/r134/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2015_pop134.csv"), row.names = F)

##//+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2015_pop134.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize
# str_out %>% filter(K == 3) %>% slice_max(MeanLnP_K, n=1)

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2015/r134/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2014_134')
write.csv(clust, file=paste0(PATH, "/structure_results/r2014_134_dK2_clusterassignments.csv"), row.names = F)

# clust <- read.delim(paste0(PATH, "/structure_results/r2015/r134/inferred_ancestry_K3.csv"), header = F, sep = ",")
# clust <- set_cluster_df(clust, '2014_134')
# write.csv(clust, file=paste0(PATH, "/structure_results/r2014_dK3_clusterassignments.csv"), row.names = F)

##//+ make ggplot ----
clust <- read.csv(paste0(PATH, "/structure_results/r2014_134_dK2_clusterassignments.csv"))
group.colors <- c("#2B847C", "#DF9F1F")
structure_plot_it(clust, group.colors, "2014, Pops 1, 3, & 4")
ggsave(paste0(PATH, "/structure_results/figures/strct_2014_134_K2.png"), plot = last_plot(), width = 10, height = 6)

# clust <- read.csv(paste0(PATH, "/structure_results/r2014_134_dK3_clusterassignments.csv"))
# group.colors <- c("#2B847C", "#DF9F1F", "#B2C4DF")
# structure_plot_it(clust, group.colors, "2014, Pops 1, 3, & 4")
# ggsave(paste0(PATH, "/structure_results/figures/strct_2014_134_K3.png"), plot = last_plot(), width = 10, height = 6)

#//pops 6, 7, 8, 9, 10 ----
##//+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2015/r678910/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2015_pop678910.csv"), row.names = F)

##//+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2015_pop678910.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize
# str_out %>% filter(K == 2) %>% slice_max(MeanLnP_K, n=1)

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2015/r678910/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2014_678910')
write.csv(clust, file=paste0(PATH, "/structure_results/r2014_678910_dK3_clusterassignments.csv"), row.names = F)

# clust <- read.delim(paste0(PATH, "/structure_results/r2015/r678910/inferred_ancestry_K2.csv"), header = F, sep = ",")
# clust <- set_cluster_df(clust, '2014_678910')
# write.csv(clust, file=paste0(PATH, "/structure_results/r2014_678910_dK2_clusterassignments.csv"), row.names = F)

##//+ make ggplot ----
clust <- read.csv(paste0(PATH, "/structure_results/r2014_678910_dK3_clusterassignments.csv"))
group.colors <- c("#2B847C", "#DF9F1F", "#B2C4DF")
structure_plot_it(clust, group.colors, "2014, Pops 6, 7, 8, 9, & 10")
ggsave(paste0(PATH, "/structure_results/figures/strct_2014_678910_K3.png"), plot = last_plot(), width = 10, height = 6)

# clust <- read.csv(paste0(PATH, "/structure_results/r2014_678910_dK2_clusterassignments.csv"))
# group.colors <- c("#2B847C", "#DF9F1F")
# structure_plot_it(clust, group.colors, "2014, Pops 6, 7, 8, 9, & 10")
# ggsave(paste0(PATH, "/structure_results/figures/strct_2014_678910_K2.png"), plot = last_plot(), width = 10, height = 6)

#2019 ----
##+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2019/iter_structure_results.csv"), header = T, sep = "\t")
str_out <- str_out[c(143:212),] %>% mutate(across(.cols = everything(), as.numeric))
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2019_maxK7.csv"), row.names = F)

##+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2019_maxK7.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize
str_out %>% filter(K == 2) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize for close dK

#convert the output final from structure to readable here using 12c_str_text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2019/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2019')
write.csv(clust, file=paste0(PATH, "/structure_results/r2019_dK3_clusterassignments.csv"), row.names = F)

#for delta K=2, which was close
# clust <- read.delim(paste0(PATH, "/structure_results/r2019/inferred_ancestry_K2.csv"), header = F, sep = ",")
# clust <- set_cluster_df(clust, '2019')
# write.csv(clust, file=paste0(PATH, "/structure_results/r2019_dK2_clusterassignments.csv"), row.names = F)

##+ make ggplot ----
# K = 3
group.colors <- c("#12263A", "#F4D1AE", "#00C6B8")
clust <- read.csv(paste0(PATH, "/structure_results/r2019_dK3_clusterassignments.csv"))
structure_plot_it(clust, group.colors, "2019")

# K = 2
clust <- read.csv(paste0(PATH, "/structure_results/r2019_dK2_clusterassignments.csv"))
group.colors <- c("#12263A", "#F4D1AE")
structure_plot_it(clust, group.colors, "2019", cluster.relevel = c("clust2", "clust1"))
ggsave(paste0(PATH, "/structure_results/figures/strct_2019_K3.png"), plot = last_plot(), width = 10, height = 6)

#//pops 6, 7, 8, 9 ----
##//+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2019/r6789/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2019_pop6789.csv"), row.names = F)

##//+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2019_pop6789.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2019/r6789/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2019_6789')
write.csv(clust, file=paste0(PATH, "/structure_results/r2019_6789_dK2_clusterassignments.csv"), row.names = F)

##//+ make ggplot ----
clust <- read.csv(paste0(PATH, "/structure_results/r2019_6789_dK2_clusterassignments.csv"))
group.colors <- c("#2B847C", "#DF9F1F")
structure_plot_it(clust, group.colors, "2019, Pops 6, 7, 8, & 9")
ggsave(paste0(PATH, "/structure_results/figures/strct_2019_6789_K2.png"), plot = last_plot(), width = 10, height = 6)

#2021 ----
##+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2021/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2021_maxK10.csv"), row.names = F)

##+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2021_maxK10.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize
str_out %>% filter(K == 3) %>% slice_max(MeanLnP_K, n=1)

#convert the output final from structure to readable here using 12c_str_text2csv.sh

#read in from bash
#K = 2
clust <- read.delim(paste0(PATH, "/structure_results/r2021/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2021')
write.csv(clust, file=paste0(PATH, "/structure_results/r2021_dK2_clusterassignments.csv"), row.names = F)

#K = 3 
clust <- read.delim(paste0(PATH, "/structure_results/r2021/inferred_ancestry_K3.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2021')
write.csv(clust, file=paste0(PATH, "/structure_results/r2021_dK3_clusterassignments.csv"), row.names = F)

##+ make ggplot ----
#K = 2
clust <- read.csv(paste0(PATH, "/structure_results/r2021_dK2_clusterassignments.csv"))
group.colors <- c("#12263A", "#F4D1AE")
structure_plot_it(clust, group.colors, "2021")
ggsave(paste0(PATH, "/structure_results/figures/strct_2021_K2.png"), plot = last_plot(), width = 10, height = 6)

#K = 3
clust <- read.csv(paste0(PATH, "/structure_results/r2021_dK3_clusterassignments.csv"))
group.colors <- c("#F4D1AE", "#12263A", "#00C6B8")
structure_plot_it(clust, group.colors, "2021")
ggsave(paste0(PATH, "/structure_results/figures/strct_2021_K3.png"), plot = last_plot(), width = 10, height = 6)

#//pops 1, 3, 4, 10, 16 ----
##//+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2021/r1341016/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2021_pop1341016.csv"), row.names = F)

##//+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2021_pop1341016.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2021/r1341016/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2021_1341016')
write.csv(clust, file=paste0(PATH, "/structure_results/r2021_1341016_dK2_clusterassignments.csv"), row.names = F)

##//+ make ggplot ----
clust <- read.csv(paste0(PATH, "/structure_results/r2021_1341016_dK2_clusterassignments.csv"))
group.colors <- c("#2B847C", "#DF9F1F")
structure_plot_it(clust, group.colors, "2021, Pops 1, 3, 4, 10, 16")
ggsave(paste0(PATH, "/structure_results/figures/strct_2021_1341016_K2.png"), plot = last_plot(), width = 10, height = 6)

#//pops 6, 7, 8, 9 ----
##//+ delta K calc----
str_out <- read.delim(paste0(PATH, "/structure_results/r2021/r6789/iter_structure_results.csv"), header = T, sep = "\t")
K_out <- delta_K(str_out)
write.csv(K_out, paste0(PATH, "/results_tables/structure_deltaK2021_pop6789.csv"), row.names = F)

##//+ read STRUCTURE output ----
#based on highest delta K
K_out <- read.csv(paste0(PATH, "/results_tables/structure_deltaK2021_pop1341016.csv"))
maxK <- K_out %>% slice_max(deltaK) %>% pull(K)
str_out %>% filter(K == maxK) %>% slice_max(MeanLnP_K, n=1) #selecting which rep to visualize

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/r2021/r6789/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust, '2021_6789')
write.csv(clust, file=paste0(PATH, "/structure_results/r2021_6789_dK2_clusterassignments.csv"), row.names = F)

##//+ make ggplot ----
clust <- read.csv(paste0(PATH, "/structure_results/r2021_6789_dK2_clusterassignments.csv"))
group.colors <- c("#2B847C", "#DF9F1F")
structure_plot_it(clust, group.colors, "2021, Pops 6, 7, 8, 9")
ggsave(paste0(PATH, "/structure_results/figures/strct_2021_6789_K2.png"), plot = last_plot(), width = 10, height = 6)
