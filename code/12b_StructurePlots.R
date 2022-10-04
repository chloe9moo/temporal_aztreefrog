##TITLE: AZ treefrog: STRUCTURE output management + plot creation
##AUTHOR: C. E. Moore
##Updated on 5 Oct 2021

library(tidyverse); library(ggh4x) #; library(adegenet); library(paletteer); library(cowplot)

## home path
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"

## load data
# source(file = paste0(PATH, "/code/01c_DataLoad.R"))
pop.key <- read.csv(paste0(PATH, "/microsat_data/year_pop_number_key4structure.csv"))
ind.key <- read.csv(paste0(PATH, "/microsat_data/ind_key4structure.csv"))

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

set_cluster_df <- function(x) {
  cluster_df <- x[,-c(1, ncol(x))] %>%
    rename_with(.cols = everything(), ~paste0("clust", seq_along(.))) %>%
    mutate(Label = row_number()) %>%
    left_join(., ind.key) %>%
    left_join(., aztf %>% select(ID, year, pop)) %>%
    relocate(ID, year, pop)
  return(cluster_df)
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

#convert the output final from structure to readable here using text2csv.sh

#read in from bash
clust <- read.delim(paste0(PATH, "/structure_results/K24/inferred_ancestry.csv"), header = F, sep = ",")
clust <- set_cluster_df(clust)
write.csv(clust, file=paste0(PATH, "/structure_results/AllInd_dK2_clusterassignments.csv"), row.names = F)

##+ make ggplot ----
long.clust <- pivot_longer(clust, cols = starts_with("clust"), names_to = "cluster", values_to = "prob")
long.clust$pondfact <- fct_relevel(as.factor(long.clust$pop), "1", "3", "4", "6", "7", "8", "9", "10", "16")
long.clust$yearfact <- fct_relevel(as.factor(long.clust$year), "2014", "2019", "2021")

#pick colors for clusters
group.colors <- c('clust1' = "#12263A", 'clust2' = "#F4D1AE")
#prep graphics
facetstrips <- strip_nested(
  text_x = elem_list_text(size = c(12, 4)),
  by_layer_x = TRUE, clip = "off"
)
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

#plot
ggplot(long.clust, aes(factor(ID), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("All years, all pops, delta K = 2")
ggsave(paste0(PATH, "/structure_results/figures/year_on_pop_strct.png"), plot = last_plot(), width = 10, height = 6)
ggsave(paste0(PATH, "/structure_results/figures/pop_on_year_strct.png"), plot = last_plot(), width = 10, height = 6)


