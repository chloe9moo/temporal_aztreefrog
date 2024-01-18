##TITLE: AZ treefrog: STRUCTURE output management + plot creation

# Take STRUCTURE output from 12a_iter_structure.sh and find most likely K based on delta K method
# Find rep within K run with highest likelihood
# >> 12c_str_text2csv.sh : Take STRUCTURE output from K x rep run and get inferred ancestry, save for use in R
# Take inferred ancestry file and plot with ggplot

library(tidyverse); library(ggh4x) #; library(adegenet); library(paletteer); library(cowplot)

## home path
PATH <- getwd()

## load data
source(file = paste0(PATH, "/code/01c_DataLoad.R"))
# pop.key <- read.csv(paste0(PATH, "/microsat_data/structure_files/year_pop_number_key4structure.csv"))

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
#cluster.assignments: should be the inferred ancestry csv already gone through set_cluster_df
#color.assignments: set which colors to have in structure plot, must equal number of clusters 
#plot.title: set plot title
#n.yrs: if TRUE, all of the years are included (adds the year under or over pops in the plot). if false, does not show a year label
#cluster.relevel: can specify how to order clusters, in case numerical order isn't wanted
  
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
    text_x = elem_list_text(size = c(12, 10)),
    by_layer_x = TRUE, clip = "off"
  )
  
  if(n.yrs == TRUE) { 
    graphic.list <- list(
      geom_col(size = 0.01, width = 1),
      facet_nested(~ yearfact+pondfact,
                   #~ pondfact+yearfact, #switch which is on top or bottom in the fig
                   switch = "x",
                   nest_line = element_line(linewidth = 1, lineend = "round"),
                   scales = "free", space = "free", strip = facetstrips),
      #theme_minimal(),
      labs(x = "", y = ""),
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
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "none"
      ))
  } else {
    graphic.list <- list(
      geom_col(size = 0.01, width = 1),
      facet_nested(~ pondfact,
                   switch = "x",
                   nest_line = element_line(linewidth = 1, lineend = "round"),
                   scales = "free", space = "free", strip = facetstrips),
      #theme_minimal(),
      labs(x = "", y = ""),
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
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "none"
      ))
  }
  
  #plot
  p <- ggplot(long.clust, aes(factor(ID), prob, fill = cluster)) +
    graphic.list #+
    #ggtitle( paste0( plot.title, "; K = ", length(unique(long.clust$cluster)) ) )
  # p
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
# clust <- read.csv(paste0(PATH, "/structure_results/AllInd_dK2_clusterassignments.csv"))

##+ make ggplot ----
group.colors <- c("#12263A", "#F4D1AE")
p.all <- structure_plot_it(clust, group.colors, "All years, All pops", n.yrs = T)
ggsave(paste0(PATH, "/structure_results/figures/pop_on_year_strct.png"), plot = p.all, width = 10, height = 6)

p14 <- structure_plot_it(clust %>% filter(year == 2014), group.colors, "", n.yrs = T)
p19 <- structure_plot_it(clust %>% filter(year == 2019), group.colors, "", n.yrs = T)
p21 <- structure_plot_it(clust %>% filter(year == 2021), group.colors, "", n.yrs = T)

cowplot::plot_grid(p14, p19, p21, align = "v", axis = "l", ncol=1) #stack them
ggsave(paste0(PATH, "/structure_results/figures/pop_on_year_strct_stacked.png"), width = 6, height = 10)
# ggsave(paste0(PATH, "/structure_results/figures/year_on_pop_strct.png"), plot = last_plot(), width = 10, height = 6)
