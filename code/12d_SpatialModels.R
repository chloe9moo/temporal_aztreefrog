##TITLE: Spatial Models
##AUTHOR: C. E. Moore
##Updated on 15 SEP 2022

#set up ----
##packages
library(adegenet); library(tidyverse); library(mmod); library(MASS); library(vegan)

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
  select(-X) %>%
  arrange(pop) %>%
  column_to_rownames("pop") %>%
  rename_with(.cols = everything(), ~ gsub("X", "", .x)) %>%
  relocate(row.names(st.d))

##+ GST ----

IBD <- data.frame(group = unlist(ind.names), mantel.r_euc = NA, p.val_euc = NA, mantel.r_str = NA, p.val_str = NA)
for (i in 1:length(ind.list)) {
  
  ##euc distance
  #pairwise distance - gst and spatially
  pop_dist <- dist(ind.list[[i]]@other$xy, method="euclidean")
  az.pw_G <- pairwise_Gst_Hedrick(ind.list[[i]], linearized = TRUE)
  
  daz.pw_G <- round(as.data.frame(as.matrix(az.pw_G)),3)
  daz.pw_G[upper.tri(daz.pw_G)] <- NA
  write.csv(daz.pw_G, file = paste0(PATH, "/results_tables/", ind.names[[i]], "_pairwiseGst.csv"))
  
  #mantel test
  ibd.s <- mantel(az.pw_G, pop_dist, method = "pearson", 999)
  IBD[IBD$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD[IBD$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD[IBD$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD[IBD$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% select(names(az.pw_G)) %>% filter(row.names(st.d) %in% names(az.pw_G)) %>% as.dist()
    
    #mantel test
    ibd.str <- mantel(az.pw_G, st.d2, method = "pearson", 999)
    IBD[IBD$group == ind.names[[i]], "mantel.r_str"] <- ibd.str$statistic
    IBD[IBD$group == ind.names[[i]], "p.val_str"] <- ibd.str$signif
    }
}

IBD #looks like 2019 is the only one where IBD is a feasible model
write.csv(IBD, file = paste0(PATH, "/results_tables/IBD_Gst.csv"), row.names = F)

#plot to check linearity
pop_dist <- dist(ind.list[[4]]@other$xy, method="euclidean")
az.pw_G <- pairwise_Gst_Hedrick(ind.list[[4]], linearized = TRUE)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
dens <- MASS::kde2d(as.vector(pop_dist), as.vector(az.pw_G), n=300)
plot(as.vector(pop_dist), as.vector(az.pw_G), pch=20, cex=0.5,  
     xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(as.vector(az.pw_G) ~ as.vector(pop_dist)))
lines(loess.smooth(as.vector(pop_dist), as.vector(az.pw_G)), col="red")

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
  
  #mantel test
  ibd.s <- mantel(pw.f, pop_dist, method = "pearson", 999)
  IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD.f[IBD.f$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD.f[IBD.f$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% select(names(pw.f)) %>% filter(row.names(st.d) %in% names(pw.f)) %>% as.dist()
    
    #mantel test
    IBD.f.str <- mantel(pw.f, st.d2, method = "pearson", 999)
    IBD.f[IBD.f$group == ind.names[[i]], "mantel.r_str"] <- IBD.f.str$statistic
    IBD.f[IBD.f$group == ind.names[[i]], "p.val_str"] <- IBD.f.str$signif
  }
  
}

IBD.f #looks like 2019 is the only one where IBD is a feasible model
write.csv(IBD.f, file = paste0(PATH, "/results_tables/IBD_Fst.csv"), row.names = F)

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
  
  #mantel test
  ibd.s <- mantel(pw.d, pop_dist, method = "pearson", 999)
  IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_euc"] <- ibd.s$statistic
  IBD.d[IBD.d$group == ind.names[[i]], "p.val_euc"] <- ibd.s$signif
  
  ##stream distance
  if (i %in% c(2, 6)) { 
    IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_str"] <- NA
    IBD.d[IBD.d$group == ind.names[[i]], "p.val_str"] <- NA
  } else { 
    st.d2 <- st.d %>% select(names(pw.d)) %>% filter(row.names(st.d) %in% names(pw.d)) %>% as.dist()
    
    #mantel test
    IBD.d.str <- mantel(pw.d, st.d2, method = "pearson", 999)
    IBD.d[IBD.d$group == ind.names[[i]], "mantel.r_str"] <- IBD.d.str$statistic
    IBD.d[IBD.d$group == ind.names[[i]], "p.val_str"] <- IBD.d.str$signif
  }
  
}

IBD.d #looks like 2019 is the only one where IBD is a feasible model
write.csv(IBD.d, file = paste0(PATH, "/results_tables/IBD_Dps.csv"), row.names = F)







