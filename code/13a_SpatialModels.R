##TITLE: Spatial Models
##AUTHOR: C. E. Moore
##Updated on 15 SEP 2022

#set up ----
##packages
library(adegenet); library(tidyverse); library(mmod); library(vegan)#; library(ecodist); library(MASS)

##home path
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"

#gene data
source(file = paste0(PATH, "/code/01c_DataLoad.R"))
#add yr_group genind
yr.coord <- coord %>% dplyr::select(year_pop, lon, lat) %>% unique() %>% remove_rownames() %>% column_to_rownames("year_pop")

ind_aztf_yr <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(aztf$year_pop)
ind_aztf_yr$other$xy <- subset(yr.coord, rownames(yr.coord) %in% aztf$year_pop) #assign coordinates

ind.list[[6]] <- ind_aztf_yr
pop.list[[6]] <- genind2genpop(ind_aztf_yr, quiet = T)
ind.names[[6]] <- "Year_Pop"
pop.names[[6]] <- "Year_Pop"

rm(yr.coord, ind_aztf_yr)

#IBD ----

st.d <- read.csv(paste0(PATH, "/results_tables/RiverDist_Matrix.csv")) %>% 
  dplyr::select(-X) %>%
  mutate(pop = case_when(pop == "1" ~ 1,
                         pop == "4" ~ 2,
                         pop == "3" ~ 3,
                         pop == "6" ~ 4, 
                         pop == "7" ~ 5, 
                         pop == "8" ~ 6, 
                         pop == "9" ~ 7, 
                         pop == "10" ~ 8, 
                         pop == "16" ~ 9)) %>%
  arrange(pop) %>%
  column_to_rownames("pop") %>%
  rename_with(.cols = everything(), ~ gsub("X", "", .x)) %>%
  rename_with(.cols = everything(), ~ case_when(.x == "1" ~ 1,
                                                .x == "4" ~ 2,
                                                .x == "3" ~ 3,
                                                .x == "6" ~ 4, 
                                                .x == "7" ~ 5, 
                                                .x == "8" ~ 6, 
                                                .x == "9" ~ 7, 
                                                .x == "10" ~ 8, 
                                                .x == "16" ~ 9)) %>%
  relocate(row.names(.))


##+ GST ----

IBD <- data.frame(group = unlist(ind.names), mantel.r_euc = NA, p.val_euc = NA, mantel.r_str = NA, p.val_str = NA)
for (i in 1:length(ind.list)) {
  
  ##euc distance
  #pairwise distance - gst and spatially
  # pop_dist <- dist(ind.list[[i]]@other$xy, method="euclidean")
  pop_dist <- log(dist(ind.list[[i]]@other$xy, method="euclidean"))
  az.pw_G <- pairwise_Gst_Hedrick(ind.list[[i]], linearized = TRUE) #linearized = Gst / (1-Gst)
  
  daz.pw_G <- round(as.data.frame(as.matrix(az.pw_G)),3)
  daz.pw_G[upper.tri(daz.pw_G)] <- NA
  write.csv(daz.pw_G, file = paste0(PATH, "/results_tables/", ind.names[[i]], "_pairwiseGst.csv"))
  
  #mantel test
  ibd.s <- vegan::mantel(az.pw_G, pop_dist, method = "pearson", 999)
  IBD[IBD$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD[IBD$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD[IBD$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD[IBD$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% dplyr::select(names(az.pw_G)) %>% filter(row.names(st.d) %in% names(az.pw_G)) %>% as.dist()
    st.d2 <- log(st.d2)
    
    #mantel test
    ibd.str <- vegan::mantel(az.pw_G, st.d2, method = "pearson", 999)
    IBD[IBD$group == ind.names[[i]], "mantel.r_str"] <- ibd.str$statistic
    IBD[IBD$group == ind.names[[i]], "p.val_str"] <- ibd.str$signif
    }
}

IBD
# write.csv(IBD, file = paste0(PATH, "/results_tables/IBD_Gst.csv"), row.names = F)
write.csv(IBD, file = paste0(PATH, "/results_tables/IBD_Gst_LogDist.csv"), row.names = F)

#plot to check IBD patterns
dist_plot_gst <- function(genind.obj, title, gen.transform = T, geo.tranform = F) {
  # genind.obj = genind object with xy info
  # gen.measure is metric to use, can be fst, gst, or dps
  # title, is title to put on the plot
  # gen.transform applies the g.dist / (1 - g.dist) transformation
  # geo.transform applies log(distance) transformation
  if (geo.tranform == T) {
    pop_dist <- log(dist(genind.obj@other$xy, method="euclidean"))
    x_title <- "log(Geographic Distance)"
  } else {
    pop_dist <- dist(genind.obj@other$xy, method="euclidean")
    x_title <- "Geographic Distance"
  }
  
  if (gen.transform == T) {
    az.pw_G <- pairwise_Gst_Hedrick(genind.obj, linearized = TRUE)
    y_title <- "Gst / (1 - Gst)"
  } else {
    az.pw_G <- pairwise_Gst_Hedrick(genind.obj, linearized = FALSE)
    y_title <- "Genetic Distance"
  }
  
  df <- data.frame(geo.dist=as.vector(pop_dist), gen.dist=as.vector(az.pw_G))
  p <- ggplot(df, aes(x=geo.dist, y=gen.dist)) +
    stat_density_2d(geom = "raster",
                    aes(fill = after_stat(density)),
                    contour = FALSE,
                    show.legend = FALSE) +
    geom_point() +
    geom_smooth(method = "lm", formula = y~x, se = F, color="black", linewidth=0.5) +
    #geom_smooth(method = "loess", formula = y~x, se = F, color="red", size=0.5) +
    scale_fill_viridis_c() +
    labs(x=x_title, y=y_title, title = title) +
    theme_minimal()
  ggsave(paste0(PATH, "/figures/Dist_x_Gst_", title, ".png"), plot = p, width = 5, height = 4) 
  return(p)
}

dist_plot_gst(ind.list[[4]], "2019 IBD", gen.transform = T, geo.tranform = T)
dist_plot_gst(ind.list[[3]], "2014 IBD")
dist_plot_gst(ind.list[[5]], "2021 IBD", geo.tranform = T)
dist_plot_gst(ind.list[[2]], "All by Year")

#base plot distance figure
# myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
# dens <- MASS::kde2d(as.vector(pop_dist), as.vector(az.pw_G), n=300)
# plot(as.vector(pop_dist), as.vector(az.pw_G), pch=20, cex=0.5,  
#      xlab="Geographic Distance", ylab="Genetic Distance")
# image(dens, col=transp(myPal(300), 0.7), add=TRUE)
# abline(lm(as.vector(az.pw_G) ~ as.vector(pop_dist)))
# lines(loess.smooth(as.vector(pop_dist), as.vector(az.pw_G)), col="red")

##+ FST ----
#bring in dist matrices made in 12_Structure.Rmd
fst_l <- list.files(path = paste0(PATH, "/results_tables"), pattern = "pairwiseFst.csv", full.names = T)
fst_l <- lapply(fst_l, read.csv)

IBD.f <- data.frame(group = unlist(ind.names), mantel.r_euc = NA, p.val_euc = NA, mantel.r_str = NA, p.val_str = NA)
for (i in 1:length(fst_l)) {
  
  pw.f <- fst_l[[i]] %>% 
    column_to_rownames(var = "X") %>%
    rename_with(.cols = everything(), ~ gsub("X", "", .x)) %>% as.dist()
  pop_dist <- dist(ind.list[[i]]@other$xy, method="euclidean")
#if log transforming distance
  pop_dist <- log(pop_dist)
  
  #mantel test
  ibd.s <- vegan::mantel(pw.f, pop_dist, method = "pearson", 999)
  IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD.f[IBD.f$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD.f[IBD.f$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% dplyr::select(names(pw.f)) %>% filter(row.names(st.d) %in% names(pw.f)) %>% as.dist()
#if log transforming distance
    st.d2 <- log(st.d2)
    
    #mantel test
    IBD.f.str <- vegan::mantel(pw.f, st.d2, method = "pearson", 999)
    IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_str"] <- IBD.f.str$statistic
    IBD.f[IBD.f$group == ind.names[[i]], "p.val_str"] <- IBD.f.str$signif
  }
  
}

IBD.f #looks like 2019, 2021 most significant
# write.csv(IBD.f, file = paste0(PATH, "/results_tables/IBD_Fst.csv"), row.names = F)
write.csv(IBD.f, file = paste0(PATH, "/results_tables/IBD_Fst_logDist.csv"), row.names = F)

##+ Dps ----
#bring in dist matrices made in 12_Structure.Rmd
Dps_l <- list.files(path = paste0(PATH, "/results_tables"), pattern = "pairwiseDps.csv", full.names = T)
Dps_l <- lapply(Dps_l, read.csv) 

IBD.d <- data.frame(group = unlist(ind.names), mantel.r_euc = NA, p.val_euc = NA, mantel.r_str = NA, p.val_str = NA)
for (i in 1:length(Dps_l)) {
  
  pw.d <- Dps_l[[i]] %>% 
    column_to_rownames(var = "X") %>%
    rename_with(.cols = everything(), ~ gsub("X", "", .x)) %>% 
    as.dist()
  pop_dist <- dist(ind.list[[i]]@other$xy, method="euclidean")
#if log transforming distance
  pop_dist <- log(pop_dist)
  
  #mantel test
  ibd.s <- vegan::mantel(pw.d, pop_dist, method = "pearson", 999)
  IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD.d[IBD.d$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD.d[IBD.d$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% dplyr::select(names(pw.d)) %>% filter(row.names(st.d) %in% names(pw.d)) %>% as.dist()
#if log transforming distance
    st.d2 <- log(st.d2)
    
    #mantel test
    IBD.d.str <- vegan::mantel(pw.d, st.d2, method = "pearson", 999)
    IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_str"] <- IBD.d.str$statistic
    IBD.d[IBD.d$group == ind.names[[i]], "p.val_str"] <- IBD.d.str$signif
  }
  
}

IBD.d #looks like 2019 is the only one where IBD is a feasible model
# write.csv(IBD.d, file = paste0(PATH, "/results_tables/IBD_Dps.csv"), row.names = F)
write.csv(IBD.d, file = paste0(PATH, "/results_tables/IBD_Dps_LogDist.csv"), row.names = F)

#+ plot IBD ----
dist_plot <- function(genind.obj, obj.name, gen.measure, gen.transform = F, geo.tranform = F) {
  # genind.obj = genind object with xy info
  # obj.name = must grab from same part of genind.obj to match files
  # gen.measure is metric to use, can be fst, gst, or dps
  # title, is title to put on the plot
  # gen.transform applies the g.dist / (1 - g.dist) transformation
  # geo.transform applies log(distance) transformation
  if (geo.tranform == T) {
    pop_dist <- log(dist(genind.obj@other$xy, method="euclidean"))
    x_title <- "log(Geographic Distance)"
  } else {
    pop_dist <- dist(genind.obj@other$xy, method="euclidean")
    x_title <- "Geographic Distance"
  }
  
  if (gen.measure == "Fst") { obj.name <- gsub(" ", "_", obj.name) } #because I named that file differently for some reason
  
  az.pw_G <- read.csv(paste0(PATH, "/results_tables/", obj.name, "_pairwise", gen.measure, ".csv")) #bring in previously calc pw matrix
  az.pw_G <- az.pw_G %>% 
    column_to_rownames(var = "X") %>%
    rename_with(.cols = everything(), ~ gsub("X", "", .x)) %>% 
    as.dist()

  if (gen.transform == T) {
    az.pw_G <- az.pw_G / (1 - az.pw_G)
    y_title <- paste0(gen.measure, " / (1 - ", gen.measure, ")")
    all.title <- paste0(gsub(".*?([0-9]+).*", "\\1", obj.name), " linearized ", gen.measure)
  } else {
    y_title <- gen.measure
    all.title <- paste0(gsub(".*?([0-9]+).*", "\\1", obj.name), " ", gen.measure)
  }
  
  df <- data.frame(geo.dist=as.vector(pop_dist), gen.dist=as.vector(az.pw_G))
  p <- ggplot(df, aes(x=geo.dist, y=gen.dist)) +
    stat_density_2d(geom = "raster",
                    aes(fill = after_stat(density)),
                    contour = FALSE,
                    show.legend = FALSE) +
    geom_point() +
    geom_smooth(method = "lm", formula = y~x, se = F, color="black", size=0.5) +
    geom_smooth(method = "loess", formula = y~x, se = F, color="red", size=0.5) +
    scale_fill_viridis_c() +
    labs(x=x_title, y=y_title, title = all.title) +
    theme_minimal()
  return(p)
}

dist_plot(ind.list[[3]], ind.names[[3]], "Dps")
ggsave(filename = paste0(PATH, "/figures/IBD_Dps_2014.png"), plot = last_plot(), width = 6, height = 6, bg = "white")
dist_plot(ind.list[[4]], ind.names[[4]], "Dps")
ggsave(filename = paste0(PATH, "/figures/IBD_Dps_2019.png"), plot = last_plot(), width = 6, height = 6, bg = "white")
dist_plot(ind.list[[5]], ind.names[[5]], "Dps")
ggsave(filename = paste0(PATH, "/figures/IBD_Dps_2021.png"), plot = last_plot(), width = 6, height = 6, bg = "white")

#GLMs ----


