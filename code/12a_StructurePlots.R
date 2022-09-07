# STRUCTURE like plot creation
library(tidyverse); library(adegenet); library(ggh4x); library(paletteer); library(cowplot)

## home path
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"

## load data
source(file = paste0(PATH, "/code/01c_DataLoad.R"))

# set graphics -----
group.colors <- c('1' = aztf.pal[1], '2' = aztf.pal[2], '3' = aztf.pal[3], '4' = aztf.pal[4], '5' = aztf.pal[5])
facetstrips <- strip_nested(
  text_x = elem_list_text(size = c(12, 4)),
  by_layer_x = TRUE, clip = "off"
)

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~ pondfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = group.colors),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ))

# 2014 ----
# find clusters
grp <- find.clusters(ind.list[[3]], max.n.clust = length(unique(ind.list[[3]]$pop)), pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_2014.png"))
plot(grp$Kstat, type="b", col="blue", main="2014 - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind.list[[3]], n.clust = 3, pca.select = "percVar", perc.pca = 95)
grp2 <- find.clusters(ind.list[[3]], n.clust = 4, pca.select = "percVar", perc.pca = 95)
grp3 <- find.clusters(ind.list[[3]], n.clust = 5, pca.select = "percVar", perc.pca = 95)

# dapc
#from tutorial: The method displays the same graph of cumulated variance as in find.cluster. However,
# unlike k-means, DAPC can benefit from not using too many PCs. Indeed, retaining too many
# components with respect to the number of individuals can lead to over-fitting and unstability
# in the membership probabilities returned by the method

#k=3
az.dapc <- dapc(ind.list[[3]], grp$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind.list[[3]], grp$grp, n.pca = 7, n.da = 10)
scatter(az.dapc)

#k=4
az.dapc2 <- dapc(ind.list[[3]], grp2$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc2) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc2 <- dapc(ind.list[[3]], grp2$grp, n.pca = 7, n.da = 10)
scatter(az.dapc2)

# k=5
az.dapc3 <- dapc(ind.list[[3]], grp3$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc3) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc3 <- dapc(ind.list[[3]], grp3$grp, n.pca = 8, n.da = 10)
scatter(az.dapc3)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2015'))$pop)

postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2015'))$pop)

postprobs3 <- as.data.frame(round(az.dapc3$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters3 <- tibble::rownames_to_column(postprobs3, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2015'))$pop)


# longer
az.long <- azclusters %>% pivot_longer(2:4, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "1", "3", "4", "6", "7", "8", "9", "10")
az.long2 <- azclusters2 %>% pivot_longer(2:5, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "1", "3", "4", "6", "7", "8", "9", "10")
az.long3 <- azclusters3 %>% pivot_longer(2:6, names_to="cluster", values_to="prob")
az.long3$pondfact <- fct_relevel(as.factor(az.long3$pond), "1", "3", "4", "6", "7", "8", "9", "10")

p1 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=3")
p2 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=4")
p3 <- ggplot(az.long3, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")
title <- ggdraw() + draw_label("2014 Samples Only", fontface = 'bold')
p <- plot_grid(title, p1, p2, p3, nrow = 4, rel_heights = c(0.1, 1, 1, 1))
ggsave(paste0(PATH,"/figures/stacked_structure_2014.png"), plot = p)

# 2019 ----
# find clusters
grp <- find.clusters(ind.list[[4]], max.n.clust = length(unique(ind.list[[4]]$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_2019.png"))
plot(grp$Kstat, type="b", col="blue", main="2019 - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind.list[[4]], n.clust = 2, pca.select = "percVar", perc.pca = 99)

#k=2
az.dapc <- dapc(ind.list[[4]], grp$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind.list[[4]], grp$grp, n.pca = 2, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2019'))$pop)
# longer
az.long <- azclusters %>% pivot_longer(2:3, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "4", "6", "7", "8", "9")

p4 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("2019 Samples Only; K=2")
ggsave(paste0(PATH,"/figures/structure_2019.png"), plot = p4, height = 7, width = 12)

# 2021 ----
# find clusters
grp <- find.clusters(ind.list[[5]], max.n.clust = length(unique(ind.list[[5]]$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_2021.png"))
plot(grp$Kstat, type="b", col="blue", main="2021 - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind.list[[5]], n.clust = 2, pca.select = "percVar", perc.pca = 99)
grp2 <- find.clusters(ind.list[[5]], n.clust = 3, pca.select = "percVar", perc.pca = 99)
grp3 <- find.clusters(ind.list[[5]], n.clust = 4, pca.select = "percVar", perc.pca = 99)

#k=2
az.dapc <- dapc(ind.list[[5]], grp$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind.list[[5]], grp$grp, n.pca = 2, n.da = 10)
scatter(az.dapc)

#k=3
az.dapc2 <- dapc(ind.list[[5]], grp2$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc2) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc2 <- dapc(ind.list[[5]], grp2$grp, n.pca = 6, n.da = 10)
scatter(az.dapc2)

# k=4
az.dapc3 <- dapc(ind.list[[5]], grp3$grp, pca.select = "percVar", perc.pca = 95, n.da = 10)
temp <- optim.a.score(az.dapc3) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc3 <- dapc(ind.list[[5]], grp3$grp, n.pca = 9, n.da = 10)
scatter(az.dapc3)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2021'))$pop)

postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2021'))$pop)

postprobs3 <- as.data.frame(round(az.dapc3$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters3 <- tibble::rownames_to_column(postprobs3, var = "ind") %>% 
  mutate(pond = (aztf %>% filter(year == '2021'))$pop)


# longer
az.long <- azclusters %>% pivot_longer(2:3, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "1", "4", "6", "7", "9", "10", "16")
az.long2 <- azclusters2 %>% pivot_longer(2:4, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "1", "4", "6", "7", "9", "10", "16")
az.long3 <- azclusters3 %>% pivot_longer(2:5, names_to="cluster", values_to="prob")
az.long3$pondfact <- fct_relevel(as.factor(az.long3$pond), "1", "4", "6", "7", "9", "10", "16")

p5 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=2")
p6 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=3")
p7 <- ggplot(az.long3, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=4")
title <- ggdraw() + draw_label("2021 Samples Only", fontface = 'bold')
pb <- plot_grid(title, p5, p6, p7, nrow = 4, rel_heights = c(0.1, 1, 1, 1))
ggsave(paste0(PATH,"/figures/stacked_structure_2021.png"), plot = pb, width = 12, height = 7)

# ALL years ----
ind_aztf_yr <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(aztf$year_pop)

# find clusters
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 5, pca.select = "percVar", perc.pca = 99)
grp2 <- find.clusters(ind_aztf_yr, n.clust = 6, pca.select = "percVar", perc.pca = 99)
grp3 <- find.clusters(ind_aztf_yr, n.clust = 7, pca.select = "percVar", perc.pca = 99)
grp4 <- find.clusters(ind_aztf_yr, n.clust = 8, pca.select = "percVar", perc.pca = 99)

#k=5
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 14, n.da = 10)
scatter(az.dapc)

#k=6
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.da = 10)
temp <- optim.a.score(az.dapc2) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.pca = 13, n.da = 10)
scatter(az.dapc2)

# k=7
az.dapc3 <- dapc(ind_aztf_yr, grp3$grp, n.da = 10)
temp <- optim.a.score(az.dapc3) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc3 <- dapc(ind_aztf_yr, grp3$grp, n.pca = 14, n.da = 10)
scatter(az.dapc3)

# k=8
az.dapc4 <- dapc(ind_aztf_yr, grp4$grp, n.da = 10)
temp <- optim.a.score(az.dapc4) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc4 <- dapc(ind_aztf_yr, grp4$grp, n.pca = 14, n.da = 10)
scatter(az.dapc4)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = aztf$pop, pond_yr = aztf$year_pop, year= aztf$year)

postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = aztf$pop, pond_yr = aztf$year_pop, year= aztf$year)

postprobs3 <- as.data.frame(round(az.dapc3$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters3 <- tibble::rownames_to_column(postprobs3, var = "ind") %>% 
  mutate(pond = aztf$pop, pond_yr = aztf$year_pop, year= aztf$year)

postprobs4 <- as.data.frame(round(az.dapc4$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters4 <- tibble::rownames_to_column(postprobs4, var = "ind") %>% 
  mutate(pond = aztf$pop, pond_yr = aztf$year_pop, year= aztf$year)

# longer
az.long <- azclusters %>% pivot_longer(2:6, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "1", "3", "4", "6", "7", "8", "9", "10", "16")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

az.long2 <- azclusters2 %>% pivot_longer(2:7, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "1", "3", "4", "6", "7", "8", "9", "10", "16")
az.long2$yearfact <- fct_relevel(as.factor(az.long2$year), "2014", "2019", "2021")

az.long3 <- azclusters3 %>% pivot_longer(2:8, names_to="cluster", values_to="prob")
az.long3$pondfact <- fct_relevel(as.factor(az.long3$pond), "1", "3", "4", "6", "7", "8", "9", "10", "16")
az.long3$yearfact <- fct_relevel(as.factor(az.long3$year), "2014", "2019", "2021")

az.long4 <- azclusters4 %>% pivot_longer(2:9, names_to="cluster", values_to="prob")
az.long4$pondfact <- fct_relevel(as.factor(az.long4$pond), "1", "3", "4", "6", "7", "8", "9", "10", "16")
az.long4$yearfact <- fct_relevel(as.factor(az.long4$year), "2014", "2019", "2021")

#set graphics for this set
group.colors <- c('1'=aztf.pal[1], '2'=aztf.pal[9], '3'=aztf.pal[2], '4'=aztf.pal[8], 
                  '5'=aztf.pal[3], '6'=aztf.pal[7], '7'=aztf.pal[4], '8'=aztf.pal[6])
facetstrips <- strip_nested(
  text_x = elem_list_text(size = c(10, 10)),
  by_layer_x = TRUE, clip = "off"
)
graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  # scale_fill_manual(values = group.colors),
  scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p8 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")
p9 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=6")
p10 <- ggplot(az.long3, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=7")
p11 <- ggplot(az.long4, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=8")

title <- ggdraw() + draw_label("All pops x years", fontface = 'bold')
pc <- plot_grid(title, p8, p9, p10, p11, nrow = 5, rel_heights = c(0.1, 0.9, 0.9, 0.9, 0.9))
ggsave(paste0(PATH,"/figures/stacked_structure_all.png"), plot = pc, width = 15, height = 8)
ggsave(paste0(PATH, "/figures/structure_all.png"), plot = p10, width = 16, height = 8)

ggplot(az.long %>% filter(pond == '7'), aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")

# hierarchical structure ----
## // pop 9 ----
sub.az <- aztf %>% filter(pop == "9")
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop9.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 2, pca.select = "percVar", perc.pca = 99)

# longer
az.long <- azclusters %>% pivot_longer(2:3, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "9")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

#k=2
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 2, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[1], aztf.pal[7])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop9 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=2")
p.pop9
ggsave(paste0(PATH, "/figures/structure_all_pop9.png"), plot = p.pop9, width = 16, height = 8)



## // pop 10 + 16 ----
sub.az <- aztf %>% filter(pop %in% c("10", "16"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop10_16.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 3, pca.select = "percVar", perc.pca = 99)

#k=3
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 5, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:4, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "10", "16")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[1], aztf.pal[7], aztf.pal[9])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop1016 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=3")
p.pop1016
ggsave(paste0(PATH, "/figures/structure_all_pop1016.png"), plot = p.pop1016, width = 16, height = 8)

## // pop 3 ----
sub.az <- aztf %>% filter(pop %in% c("3"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop3.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 2, pca.select = "percVar", perc.pca = 99)

#k=3
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 5, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:3, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "3")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[9], aztf.pal[7])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop3 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=2; most likely K is 1")
p.pop3
ggsave(paste0(PATH, "/figures/structure_all_pop3.png"), plot = p.pop3, width = 16, height = 8)

## // pop 1 + 4 ----
sub.az <- aztf %>% filter(pop %in% c("1", "4"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop1_4.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 2, pca.select = "percVar", perc.pca = 99)

#k=2
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 2, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:3, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "1", "4")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[9], aztf.pal[7], aztf.pal[4])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop14 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=2")
p.pop14
ggsave(paste0(PATH, "/figures/structure_all_pop1-4.png"), plot = p.pop14, width = 16, height = 8)

## // pop 6 + 7 + 8 ----
sub.az <- aztf %>% filter(pop %in% c("6", "7", "8"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop678.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 3, pca.select = "percVar", perc.pca = 99)

#k=2
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 16, n.da = 10)
scatter(az.dapc)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:4, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "6", "7", "8")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[9], aztf.pal[7], aztf.pal[8])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop678 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=3")
p.pop678
ggsave(paste0(PATH, "/figures/structure_all_pop678.png"), plot = p.pop678, width = 16, height = 8)

## // pop 1 + 4 + 6 + 7 + 8 ----
sub.az <- aztf %>% filter(pop %in% c("1", "4", "6", "7", "8"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop14678.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 4, pca.select = "percVar", perc.pca = 99)
grp2 <- find.clusters(ind_aztf_yr, n.clust = 5, pca.select = "percVar", perc.pca = 99)

#k=4
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 13, n.da = 10)
scatter(az.dapc)

#k=5
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.da = 10)
temp <- optim.a.score(az.dapc2) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.pca = 13, n.da = 10)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# membership probabilities
postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:5, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "1", "4", "6", "7", "8")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

az.long2 <- azclusters2 %>% pivot_longer(2:6, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "1", "4", "6", "7", "8")
az.long2$yearfact <- fct_relevel(as.factor(az.long2$year), "2014", "2019", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[9], aztf.pal[7], aztf.pal[8], aztf.pal[1], aztf.pal[2])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop14678 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=4")
p.pop14678
p.pop14678.k5 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")
p.pop14678.k5
ggsave(paste0(PATH, "/figures/structure_all_pop678.png"), plot = p.pop678, width = 16, height = 8)

## // pop 6 + 7 + 8 + 9 ----
sub.az <- aztf %>% filter(pop %in% c("6", "7", "8", "9"))
ind_aztf_yr <- df2genind(sub.az[,c(6:22)], sep = "/", ind.names = sub.az$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(sub.az$year_pop)
grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))#, pca.select = "percVar", perc.pca = 95)

png(paste0(PATH, "/figures/BIC_plot_All_pop6789.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 4, pca.select = "percVar", perc.pca = 99)
grp2 <- find.clusters(ind_aztf_yr, n.clust = 5, pca.select = "percVar", perc.pca = 99)

#k=4
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 13, n.da = 10)

#k=5
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.da = 10)
temp <- optim.a.score(az.dapc2)
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.pca = 14, n.da = 10)


# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = sub.az$pop, pond_yr = sub.az$year_pop, year= sub.az$year)

# longer
az.long <- azclusters %>% pivot_longer(2:5, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "6", "7", "8", "9")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2014", "2019", "2021")

az.long2 <- azclusters2 %>% pivot_longer(2:6, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "6", "7", "8", "9")
az.long2$yearfact <- fct_relevel(as.factor(az.long2$year), "2014", "2019", "2021")

graphic.list <- list(
  geom_col(color="black", size = 0.01, width = 1),
  facet_nested(~pondfact+yearfact,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips),
  theme_minimal(),
  labs(x = "Individuals", y = "membership probability"),
  scale_y_continuous(expand = c(0, 0)),
  scale_x_discrete(expand = expansion(add = 0.5)),
  scale_fill_manual(values = c(aztf.pal[9], aztf.pal[4], aztf.pal[8], aztf.pal[3], aztf.pal[1])),
  # scale_fill_viridis_d(option = "magma"),
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=7),
    axis.title.x = element_blank()
  ))

p.pop6789.k4 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=4")
p.pop6789.k4
p.pop6789.k5 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")
p.pop6789.k5
ggsave(paste0(PATH, "/figures/structure_all_pop6789_k4.png"), plot = p.pop6789.k4, width = 16, height = 8)
ggsave(paste0(PATH, "/figures/structure_all_pop6789_k5.png"), plot = p.pop6789.k5, width = 16, height = 8)

# only ponds with all 3 years ----
# find clusters
aztf.a <- aztf %>% filter(pop %in% c("4", "6", "7", "9"))
ind_aztf_yr <- df2genind(aztf.a[,c(6:22)], sep = "/", ind.names = aztf.a$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(aztf.a$year_pop)

grp <- find.clusters(ind_aztf_yr, max.n.clust = length(unique(ind_aztf_yr$pop)))

png(paste0(PATH, "/figures/BIC_plot_acrossyrs.png"))
plot(grp$Kstat, type="b", col="blue", main="All yrs - Value of BIC vs K Clusters", xlab = "Number of clusters", ylab="BIC")
dev.off()

grp <- find.clusters(ind_aztf_yr, n.clust = 4, pca.select = "percVar", perc.pca = 99)
grp2 <- find.clusters(ind_aztf_yr, n.clust = 5, pca.select = "percVar", perc.pca = 99)
grp3 <- find.clusters(ind_aztf_yr, n.clust = 6, pca.select = "percVar", perc.pca = 99)

#k=4
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.da = 10)
temp <- optim.a.score(az.dapc) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc <- dapc(ind_aztf_yr, grp$grp, n.pca = 13, n.da = 10)
scatter(az.dapc)

#k=5
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.da = 10)
temp <- optim.a.score(az.dapc2) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc2 <- dapc(ind_aztf_yr, grp2$grp, n.pca = 13, n.da = 10)
scatter(az.dapc2)

# k=6
az.dapc3 <- dapc(ind_aztf_yr, grp3$grp, n.da = 10)
temp <- optim.a.score(az.dapc3) 
#^displays plot that recommends a number of PCs to retain based on proportion of successful reassignment of the analysis 
#and values obtained using random groups/random discrimination
az.dapc3 <- dapc(ind_aztf_yr, grp3$grp, n.pca = 15, n.da = 10)
scatter(az.dapc3)

# membership probabilities
postprobs <- as.data.frame(round(az.dapc$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(pond = aztf.a$pop, pond_yr = aztf.a$year_pop, year= aztf.a$year)

postprobs2 <- as.data.frame(round(az.dapc2$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters2 <- tibble::rownames_to_column(postprobs2, var = "ind") %>% 
  mutate(pond = aztf.a$pop, pond_yr = aztf.a$year_pop, year= aztf.a$year)

postprobs3 <- as.data.frame(round(az.dapc3$posterior, 4))
# put probs in tibble w/ IDs and labels for sites
azclusters3 <- tibble::rownames_to_column(postprobs3, var = "ind") %>% 
  mutate(pond = aztf.a$pop, pond_yr = aztf.a$year_pop, year= aztf.a$year)


# longer
az.long <- azclusters %>% pivot_longer(2:5, names_to="cluster", values_to="prob")
az.long$pondfact <- fct_relevel(as.factor(az.long$pond), "4", "6", "7", "9")
az.long$yearfact <- fct_relevel(as.factor(az.long$year), "2015", "2019", "2021")

az.long2 <- azclusters2 %>% pivot_longer(2:6, names_to="cluster", values_to="prob")
az.long2$pondfact <- fct_relevel(as.factor(az.long2$pond), "4", "6", "7", "9")
az.long2$yearfact <- fct_relevel(as.factor(az.long2$year), "2015", "2019", "2021")

az.long3 <- azclusters3 %>% pivot_longer(2:7, names_to="cluster", values_to="prob")
az.long3$pondfact <- fct_relevel(as.factor(az.long3$pond), "4", "6", "7", "9")
az.long3$yearfact <- fct_relevel(as.factor(az.long3$year), "2015", "2019", "2021")

p12 <- ggplot(az.long, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=4")
p13 <- ggplot(az.long2, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=5")
p14 <- ggplot(az.long3, aes(factor(ind), prob, fill = factor(cluster))) +
  graphic.list +
  ggtitle("K=6")
title <- ggdraw() + draw_label("Multiyear Samples Only", fontface = 'bold')
pd <- plot_grid(title, p12, p13, p14, nrow = 4, rel_heights = c(0.1, 1, 1, 1))
ggsave(paste0(PATH,"/figures/stacked_structure_multiyear.png"), plot = pd, width = 16, height = 8)
