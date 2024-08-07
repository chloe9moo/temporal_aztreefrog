##TITLE: DAPC Analyses
##AUTHOR: C. E. Moore
##Updated on 15 SEP 2022

library(tidyverse); library(adegenet); library(ggpubr)
# library(ggh4x); library(paletteer); library(cowplot)

#set up ----
##home path
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"

##load data
source(file = paste0(PATH, "/code/01c_DataLoad.R"))

#redo the yr grouping genind
yr.coord <- coord %>% dplyr::select(year_pop, lon, lat) %>% unique() %>% remove_rownames() %>% column_to_rownames("year_pop")

ind_aztf_yr <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(aztf$year_pop)
ind_aztf_yr$other$xy <- subset(yr.coord, rownames(yr.coord) %in% aztf$year_pop) #assign coordinates

ind.list[[6]] <- ind_aztf_yr
pop.list[[6]] <- genind2genpop(ind_aztf_yr, quiet = T)
ind.names[[6]] <- "Year_Pop"
pop.names[[6]] <- "Year_Pop"

rm(yr.coord, ind_aztf_yr)

#find clusters ----
#will need to do this one by one
p <- 6 #select set to run

#no real reason to choose low number of pc's
grp <- find.clusters(ind.list[[p]], max.n.clust = length(levels(ind.list[[p]]@pop))+1, n.pca = 200)
grp$Kstat

table.value(table(pop(ind.list[[p]]), grp$grp), col.lab=paste("inf", 1:24),
            row.lab=paste("ori", 1:24))

##+ K BIC plot ----
#grab bic values for each value of K, make a dataframe and turn it into a ggplot
bic.m <- matrix(nrow = length(levels(ind.list[[p]]@pop))+1, ncol = length(levels(ind.list[[p]]@pop))+1)
for(i in 1:length(levels(ind.list[[p]]@pop))+1){
  grp <- find.clusters(ind.list[[p]], n.pca = 200, choose.n.clust = FALSE,  max.n.clust = length(levels(ind.list[[p]]@pop))+1)
  bic.m[i,] <- grp$Kstat
}

bic.m <- bic.m %>%
  as.data.frame() %>%
  rownames_to_column(., var = "Group") %>%
  pivot_longer(cols=starts_with("V"), names_to = "K", values_to = "BIC") %>%
  mutate(K = as.factor(as.numeric(gsub("V", "", K)))) %>%
  arrange(K)
p1 <- ggplot(bic.m, aes(x = K, y = BIC)) + geom_boxplot() + theme_bw() + xlab("Number of groups (K)")
p1
ggsave(paste0(PATH, "/figures/DAPC_Year_Pop_BIC.png"), plot = p1, width = 7, height = 5)
#retaining too many components with respect to the number of individuals can lead to over-fitting and unstability
##in the membership probabilities returned by the method
grp <- find.clusters(ind.list[[p]], n.clust = 7, n.pca = 200)
dapc1 <- dapc(ind.list[[p]], grp$grp)
#look at plots
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:24, cstar=0, scree.pca=TRUE,
        posi.pca="bottomleft")
assignplot(dapc1, subset=1:50)

dapc1 <- dapc(ind.list[[p]], grp$grp, n.da = 100, n.pca = 80)

##+ optimal pcs to retain based on a score ----
dapc100 <- dapc(ind.list[[p]], grp$grp, n.da = 100, n.pca = 200)
temp <- optim.a.score(dapc100)
dapc1 <- dapc(ind.list[[p]], grp$grp, n.da=100, n.pca = temp$best)

scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:24, cstar=0, scree.pca=TRUE,
        posi.pca="bottomleft")

#plot k=7 change over time ----
temp <- as.data.frame(dapc1$ind.coord)
temp$Group <- dapc1$grp
dapc2plot <- temp %>% 
  mutate(ID = row.names(.)) %>% 
  left_join(aztf, by = "ID") %>% 
  select(ID, year, pop, Group, starts_with("LD")) %>%
  mutate(pop = as.factor(as.numeric(pop)),
         year = fct_relevel(factor(year), c("2014", "2019", "2021")),
         year_pop = as.factor(paste0(year, "_", pop)))
write.csv(dapc2plot, file=paste0(PATH, "/results_tables/DAPC_res/", ind.names[[p]], "_K7.csv"), row.names = F)

dapc2plot <- read.csv(paste0(PATH, "/results_tables/DAPC_res/", ind.names[[p]], "_K7.csv")) %>%
  mutate(pop = as.factor(as.numeric(pop)),
         year = fct_relevel(factor(year), c("2014", "2019", "2021")),
         year_pop = as.factor(paste0(year, "_", pop)))

##+ get year_pop group center ----
centroids <- aggregate(cbind(LD1, LD2, LD3, LD4, LD5, LD6) ~ year_pop, dapc2plot, mean) %>%
  left_join(., dapc2plot %>% select(year_pop, pop, year) %>% distinct())

##+ set theme ----
grays <- c("#CCCCCC", "#999999", "#333333")
pca.theme <- list(
  scale_shape_manual(values = c(21, 22, 23)),
  theme(legend.position = "right", 
        panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "gray", linetype = 4),
        axis.title = element_text(size = 23),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key = element_blank(),
        plot.title = element_text(face="bold",size=23)),
  scale_y_continuous(breaks = c(-5, 0, 5), limits = c(-10, 5)) )
##set axis to plot
axis1 <- "LD1"
axis2 <- "LD2"

ggplot() +
  #centroids + path, ellipse for pop
  # stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, group=pop), geom="polygon", type = "t", alpha = 0.4, size=1, level = 0.9, show.legend = F) +
  # # stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], group=Group), color="white", type = "t", alpha = 0.95, linewidth=0.7, level = 0.8) +
  # geom_point(data = centroids %>% filter(pop == 4), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
  # geom_path(data = centroids %>% filter(pop == 4), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
  # geom_point(data = centroids %>% filter(pop == 5), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
  # geom_path(data = centroids %>% filter(pop == 5), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
  # geom_point(data = centroids %>% arrange(pop) %>% filter(!pop %in% c(4, 5)), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
  # geom_path(data = centroids %>% arrange(year) %>% filter(!pop %in% c(4, 5)), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
  # scale_fill_manual(values = aztf.pal, name = "Pop") +
  # scale_color_manual(values = aztf.pal, name = "Pop") +
  #points for every pop
  stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], color=pop), type = "t", level = 0.9, linewidth = 2, show.legend = F) +
  geom_point(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=3, alpha=0.6) +
  scale_fill_manual(values = aztf.pal, name = "Pop") +
  scale_color_manual(values = aztf.pal, name = "Pop") +
  # stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], color=year), type = "t", level = 0.9, linewidth = 2, show.legend = F) +
  # geom_point(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], fill = year, shape = year), size=3, alpha=0.6) +
  # scale_fill_manual(values = grays, name = "") +
  # scale_color_manual(values = grays, name = "") +
  scale_shape_manual(values = c(21, 22, 23), name = "Year") +
  guides(shape = guide_legend(override.aes=list(size = 5))) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 5, alpha = 1)), color="none") +
  xlim(c(-5.1,5.1)) + ylim(-6,4) +
  theme(legend.position = "right", 
        panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "gray", linetype = 4),
        legend.key = element_blank(),
        legend.text = element_text(family = "helvetica", size = 12),
        axis.text = element_text(family = "helvetica", size = 12),
        axis.title = element_text(family = "helvetica", size = 12))
# ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_OverTime.png"), plot = last_plot(), width = 7, height = 6)
ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_OverTime_nopath.png"), plot = last_plot(), width = 7, height = 6)
# ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_YearGrps.png"), plot = last_plot(), width = 7, height = 6)
# ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_YearGrps_ellipse.png"), plot = last_plot(), width = 7, height = 6)
# ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_PopGrps.png"), plot = last_plot(), width = 7, height = 6)
# ggsave(filename = paste0(PATH, "/figures/DAPC_All_K7_indivs.png"), plot = last_plot(), width = 6, height = 6)

#+ scree plots ----
pc.eig <- 100 * cumsum(dapc100$pca.eig) / sum(dapc100$pca.eig)
pc.eig <- data.frame(n_pca = seq(1, length(pc.eig), 1), eig = pc.eig) %>%
  mutate(include = case_when(n_pca < 16 ~ "IN", T ~ "OUT"))
ggplot() +
  geom_col(data = pc.eig, aes(x = n_pca, y = eig, fill = include), show.legend = F) +
  labs(title = "PCA eigenvalues", x = "", y = "") +
  scale_fill_manual(values = c("IN" = "black", "OUT" = "grey")) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste0(PATH, "/figures/DAPC_All_K7_PCAeig.png"), plot = last_plot(), width = 6, height = 5)

da.eig <- data.frame(n_da = seq(1, length(dapc1$eig), 1), eig = dapc1$eig) %>%
  mutate(include = case_when(n_da < 3 ~ "IN", T ~ "OUT"))
ggplot() +
  geom_col(data = da.eig, aes(x = n_da, y = eig, fill = include), show.legend = F) +
  labs(title = "DA eigenvalues", x = "", y = "") +
  scale_fill_manual(values = c("IN" = "#5A5A5A", "OUT" = "lightgrey")) +
  theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste0(PATH, "/figures/DAPC_All_K7_DAeig.png"), plot = last_plot(), width = 6, height = 5)

##+ calculate distance between points ----
#sqrt( (LD1 - LD1)^2 + (LD2 - LD2)^2)
dapc_dist <- centroids %>%
  arrange(pop, year) %>%
  group_by(pop) %>%
  mutate(across(contains("LD"), ~ (.x - lag(.x, default = NA))^2, .names = "{.col}_delt2"), #(deltaLD)^2
         delta_yrs = paste0(lag(year, default = NA), "_", year)) %>%
  select(pop, contains("delt")) %>%
  filter(!is.na(LD1_delt2)) %>%
  ungroup() %>%
  rowwise(pop, delta_yrs) %>%
  summarise(dist = sqrt( sum( c_across( contains("delt2") ) ) ))
long_dist <- centroids %>%
  filter(year != 2019) %>%
  arrange(pop, year) %>%
  group_by(pop) %>%
  mutate(across(contains("LD"), ~ (.x - lag(.x, default = NA))^2, .names = "{.col}_delt2"), #(deltaLD)^2
         delta_yrs = paste0(lag(year, default = NA), "_", year)) %>%
  select(pop, contains("delt")) %>%
  filter(!is.na(LD1_delt2)) %>%
  ungroup() %>%
  rowwise(pop, delta_yrs) %>%
  summarise(dist = sqrt( sum( c_across( contains("delt2") ) ) ))
dapc_dist <- bind_rows(dapc_dist, long_dist) %>% distinct() %>% arrange(pop)
write.csv(dapc_dist, file = paste0(PATH, "/results_tables/DAPC_distance_traveled.csv"), row.names = F)

dapc_dist %>% filter(delta_yrs == "2014_2021") %>% View()

#average change
dapc_dist %>% group_by(pop) %>% mutate(m = mean(dist)) %>% select(pop, m) %>% distinct() %>% View()

#multiple Ks iterations ----
my_k <- c(6, 8)

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(ind.list[[p]], n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(ind.list[[p]], pop = grp_l[[i]]$grp, n.pca = 70, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

##+ plotting ----
###++ make plot list ----
# to compare clusters
dplot_l <- vector(mode = "list", length = length(my_k))

for (i in 1:length(dapc_l)) {

  temp <- as.data.frame(dapc_l[[i]]$ind.coord)
  temp$Group <- dapc_l[[i]]$grp
  
  if (p %in% 3:5) { #single year genind objects
    dapc2plot <- temp %>% 
      mutate(ID = row.names(.)) %>% 
      left_join(aztf, by = "ID") %>% 
      select(ID, year, pop, Group, starts_with("LD")) %>%
      mutate(pop = factor(pop),
             year = factor(year))
  } else { #all year genind objects
    dapc2plot <- temp %>% 
      mutate(ID = row.names(.)) %>% 
      left_join(aztf, by = "ID") %>% 
      select(ID, year, pop, Group, starts_with("LD")) %>%
      mutate(year_pop = paste0(year, "_", pop),
             pop = as.factor(as.numeric(pop)),
             year = fct_relevel(factor(year), c("2014", "2019", "2021")))
    }
  
  write.csv(dapc2plot, file=paste0(PATH, "/results_tables/DAPC_res/", ind.names[[p]], "_K", my_k[[i]],".csv"), row.names = F)
  
  if(my_k[[i]] == 6) {
    centroids <- aggregate(cbind(LD1, LD2, LD3, LD4, LD5) ~ year_pop, dapc2plot, mean) %>%
      left_join(., dapc2plot %>% select(year_pop, pop, year) %>% distinct()) }
  if(my_k[[i]] == 8) {
    centroids <- aggregate(cbind(LD1, LD2, LD3, LD4, LD5, LD6, LD7) ~ year_pop, dapc2plot, mean) %>%
      left_join(., dapc2plot %>% select(year_pop, pop, year) %>% distinct()) }
  
  #calculate distance traveled
  dapc_dist <- centroids %>%
    arrange(pop, year) %>%
    group_by(pop) %>%
    mutate(across(contains("LD"), ~ (.x - lag(.x, default = NA))^2, .names = "{.col}_delt2"), #(deltaLD)^2
           delta_yrs = paste0(lag(year, default = NA), "_", year)) %>%
    select(pop, contains("delt")) %>%
    filter(!is.na(LD1_delt2)) %>%
    ungroup() %>%
    rowwise(pop, delta_yrs) %>%
    summarise(dist = sqrt( sum( c_across( contains("delt2") ) ) ))
  long_dist <- centroids %>%
    filter(year != 2019) %>%
    arrange(pop, year) %>%
    group_by(pop) %>%
    mutate(across(contains("LD"), ~ (.x - lag(.x, default = NA))^2, .names = "{.col}_delt2"), #(deltaLD)^2
           delta_yrs = paste0(lag(year, default = NA), "_", year)) %>%
    select(pop, contains("delt")) %>%
    filter(!is.na(LD1_delt2)) %>%
    ungroup() %>%
    rowwise(pop, delta_yrs) %>%
    summarise(dist = sqrt( sum( c_across( contains("delt2") ) ) ))
  dapc_dist <- bind_rows(dapc_dist, long_dist) %>% distinct() %>% arrange(pop)
  write.csv(dapc_dist, file = paste0(PATH, "/results_tables/DAPC_distance_traveled_K", my_k[[i]], ".csv"), row.names = F)
  
  if ("LD2" %in% colnames(dapc2plot) == FALSE) { 
    
    pca <- ggplot() + 
      geom_density(data = dapc2plot %>% group_by(Group), aes(x=LD1, fill=Group, group=Group), color="black", alpha=0.5, show.legend = F) +
      # geom_tile(data = dapc2plot, aes(x=LD1, y=0, color=pop), size=2, alpha=0.5, height=0.05) +
      geom_segment(data=dapc2plot, aes(x=LD1, xend=LD1, y=-0.01, yend=0.05, color=pop), size=3, alpha=0.8, show.legend = F) +
      labs(x="LD1", y="Density", title = paste0("K = ", my_k[i])) +
      scale_fill_manual(values = grays[c(1,3)]) +
      scale_color_manual(values = aztf.pal, aesthetics = "color") +
      theme(legend.position = "none", 
            panel.background = element_blank(), 
            panel.grid.major = element_line(colour = "gray", linetype = 4))
    
    } else {  

    pca <- ggplot() +
      stat_ellipse(data = dapc2plot, aes(x=LD1, y=LD2, fill=pop, group=pop), geom="polygon", type = "t", alpha = 0.5, size=1, level = 0.9) +
      #stat_ellipse(data = dapc2plot, aes(x=LD1, y=LD2, color=Group, group=Group), type = "t", alpha = 0.95, size=1, level = 0.9) +
      #geom_point(data = dapc2plot, aes(x=LD1, y=LD2, fill=pop, shape=year), color="black", size=2.3, alpha=0.7) +
      # geom_point(data = dapc2plot %>% filter(Group == 2), aes(x=LD1, y=LD2, fill=pop, shape=year), color="black", size=4.5, alpha=0.5) +
      geom_point(data = centroids %>% filter(pop == 4), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
      geom_path(data = centroids %>% filter(pop == 4), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
      geom_point(data = centroids %>% filter(pop == 5), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
      geom_path(data = centroids %>% filter(pop == 5), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
      geom_point(data = centroids %>% arrange(pop) %>% filter(!pop %in% c(4, 5)), aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
      geom_path(data = centroids %>% arrange(year) %>% filter(!pop %in% c(4, 5)), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.08, "inches"))) +
      labs(x="LD1", y="LD2", fill="Pop", color="Group", shape="Year", title=paste0("K = ", my_k[i])) +
      scale_fill_manual(values = aztf.pal, aesthetics = "fill") +
      # scale_color_manual(values = aztf.pal, aesthetics = "color") +
      guides(fill = guide_legend(override.aes=list(shape=21)), color="none") +
      scale_shape_manual(values = c(21, 22, 23)) +
      theme(legend.position = "right", 
            panel.background = element_blank(), 
            panel.grid.major = element_line(colour = "gray", linetype = 4),
            legend.key = element_blank())
  
    }
  
  dplot_l[[i]] <- pca

# pca2 <- ggplot() +
#   stat_ellipse(data = pca.all, aes(x=Axis1, y=Axis2, color=year, group=year), type = "t", alpha = 0.95, size=2, level = 0.9) +
#   geom_point(data = pca.all, aes(x=Axis1, y=Axis2, fill=year, shape=year), color="black", size=4.5, alpha=0.5) +
#   labs(x="PC1 (6.16%)", y="PC2 (4.05%)", fill="Year", color="Year", shape="Year") +
#   scale_fill_manual(values = grays, aesthetics = "fill") +
#   scale_color_manual(values = grays, aesthetics = "color") +
#   pca.theme
# # pca2
}

ggarrange(plotlist=dplot_l, 
          #p1,
          ncol = 2, nrow = round(length(dplot_l)/2), common.legend = TRUE, legend="right")

ggsave(filename = paste0(PATH, "/figures/DAPC_", ind.names[[p]], "_clustercomp.png"), plot = last_plot(), width = 10, height = 7)

# change over time dapc ----
# only relevant to year_pop grouping!
#not using k clustering to assess groups
dapc.o <- dapc(ind.list[[6]], n.da = 100, n.pca = 100)
temp <- optim.a.score(dapc.o)
dapc.o <- dapc(ind.list[[6]], n.da = 100, n.pca = temp$best)

#cross validation pc selection
# temp <- xvalDapc(tab(ind.list[[6]], NA.method = "mean"), pop(ind.list[[6]]))
# dapc.o <- dapc(ind.list[[6]], n.da = 100, n.pca = as.numeric(temp$`Number of PCs Achieving Highest Mean Success`))

#make a dataframe from dapc output
temp <- as.data.frame(dapc.o$ind.coord)
temp$Group <- dapc.o$grp

#make data frame legible 
dapc2plot <- temp %>% 
  mutate(ID = row.names(.)) %>% 
  left_join(aztf, by = "ID") %>% 
  select(ID, year, pop, Group, starts_with("LD")) %>%
  mutate(pop = as.factor(as.numeric(pop)),
         year = fct_relevel(factor(year), c("2014", "2019", "2021")))

#save
write.csv(dapc2plot, file=paste0(PATH, "/results_tables/DAPC_res/Year_Pop_Grp.csv"), row.names = F)
dapc2plot <- read_csv(paste0(PATH, "/results_tables/DAPC_res/Year_Pop_Grp.csv"))

#plot it
##set axis to plot
axis1 <- "LD1"
axis2 <- "LD2"
##get group center
centroids <- as.data.frame(dapc.o$grp.coord) %>% mutate(year_pop = row.names(.)) %>%
  left_join(., dapc2plot %>% select(year, pop, Group), by=c("year_pop" = "Group")) %>% 
  distinct() %>%
  complete(., year, pop, fill = list(V1 = NA, V1 = NA)) %>%
  group_by(pop)

ggplot() +
  stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, group=pop), geom="polygon", type = "t", alpha = 0.4, size=1, level = 0.9, show.legend = F) +
  # stat_ellipse(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], color=pop, group=Group), type = "t", alpha = 0.95, size=0.5, level = 0.9) +
##plot with all individuals, label for year_pop
  # geom_point(data = dapc2plot, aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=2.3, alpha=0.3, show.legend = F) +
  # geom_label(data = centroids, aes(x=V1, y=V2, label=centroids[,1])) +
  # scale_color_manual(values = aztf.pal, aesthetics = "color") +
  # scale_shape_manual(values = c(16, 15, 17)) +
##plot with only year_pop centroid
  geom_point(data = centroids, aes(x=.data[[axis1]], y=.data[[axis2]], fill=pop, shape=year), size=5, alpha=1) +
  geom_path(data = centroids, aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.075, "inches"))) +
  geom_path(data = centroids %>% filter(!(year == "2019" & is.na(year_pop))), aes(x=.data[[axis1]], y=.data[[axis2]], group = pop), arrow = arrow(type = "closed", length=unit(0.075, "inches"))) +
  scale_fill_manual(values = aztf.pal, aesthetics = "fill") +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(fill = guide_legend(override.aes=list(shape=21)), color="none") +
##use in every plot
  labs(x=axis1, y=axis2, fill="Pond", color="Group", shape="Year", title="Population change over time") +
  theme(legend.position = "right", 
        panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "gray", linetype = 4),
        legend.key = element_blank())
ggsave(filename = paste0(PATH, "/figures/DAPC_All_OverTime.png"), plot = last_plot(), width = 7, height = 7)

#first discriminant function only density plot
dens.ld1 <- list(
  labs(x="LD1", y="Density"),
  xlim(-4, 4),
  ylim(0, 0.85),
  scale_fill_manual(values = aztf.pal, aesthetics = "fill"),
  scale_color_manual(values = aztf.pal, aesthetics = "color"),
  theme(legend.position = "none", 
        panel.background = element_blank(), 
        panel.grid.major = element_line(color = "gray", linetype = 4))
)
p14 <- ggplot() + 
  geom_density(data = dapc2plot %>% filter(year == "2014") %>% group_by(pop), 
               aes(x=LD1, fill=pop, group=pop, color=pop), size=1.2, alpha=0.5, show.legend = F) +
  dens.ld1 +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x="", y="")
p19 <- ggplot() + 
  geom_density(data = dapc2plot %>% filter(year == "2019") %>% group_by(pop), 
               aes(x=LD1, fill=pop, group=pop, color=pop), size=1.2, alpha=0.5, show.legend = F) +
  dens.ld1 +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + xlab("")
p21 <- ggplot() + 
  geom_density(data = dapc2plot %>% filter(year == "2021") %>% group_by(pop), 
               aes(x=LD1, fill=pop, group=pop, color=pop), size=1.2, alpha=0.5, show.legend = F) +
  dens.ld1 +
  ylab("")
ggarrange(p14, p19, p21, ncol = 1, labels = c("2014", "2019", "2021"), label.x = 0.04, font.label = list(size=18))
ggsave(filename = paste0(PATH, "/figures/DAPC_All_LD1only.png"), plot=last_plot(), height = 8, width = 8)

# look at clustering another way ----
# my_pal <- RColorBrewer::brewer.pal(n=10, name = "Spectral")
# 
# my_df <- as.data.frame(dapc_l[[ 6 ]]$ind.coord)
# my_df$Group <- dapc_l[[ 6 ]]$grp
# 
# p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) + 
#   geom_point(size = 4, shape = 21) + 
#   theme_bw() + 
#   scale_color_manual(values=c(my_pal)) + 
#   scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
# p2
