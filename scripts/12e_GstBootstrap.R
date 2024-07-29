#significance of pGst
##packages
library(adegenet); library(tidyverse); library(mmod); library(vegan)#; library(ecodist); library(MASS)
library(parallel)

##home path
PATH <- getwd()

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

year_pop <- ind.list[[6]]

rm(yr.coord, ind_aztf_yr, ind.list, ind.names, coord, pop.list, pop.names)

#bootstrap
bs <- chao_bootstrap(year_pop, nreps = 10000)
cat("I've got my boots on.\n")

obs <- pairwise_Gst_Hedrick(bs$obs, linearized = TRUE)

#make empty list
stats <- vector("list", length = length(bs$BS))

#run parallel
# stats <- mclapply(bs$BS[1:5], function(x) pairwise_Gst_Hedrick(x, linearized = TRUE), mc.cores = 5)
stats <- mclapply(bs$BS, function(x) pairwise_Gst_Hedrick(x, linearized = TRUE), mc.cores = 8)

#bind
p.gsts <- do.call(rbind, stats)
write_csv(as.data.frame(p.gsts), paste0(PATH, "/results_tables/pGst_boots_saved.csv"))
cat("The boots are off.\n")

obs.df <- as.data.frame(as.matrix(obs))
obs.df[upper.tri(obs.df)] <- NA
obs.df <- obs.df[-1,-23]

seq.list <- vector("list", length=ncol(obs.df))
for (c in 1:ncol(obs.df)) {
  if(c == 1) { seq.list[[c]] <- as.vector(seq(1, 22)) } else {
    seq.list[[c]] <- as.vector(seq(seq.list[[c-1]][length(seq.list[[c-1]])]+1, seq.list[[c-1]][length(seq.list[[c-1]])]+(23-c)))
  } 
}
mat.name.index <- data.frame(mat.ind = seq(1, ncol(p.gsts)), top = NA, side = NA)
for (c in 1:ncol(obs.df)) {
  mat.name.index[seq.list[[c]],]$top <- names(obs.df)[c]
  mat.name.index[seq.list[[c]],]$side <- row.names(obs.df)[c:22]
}
mat.name.index <- mat.name.index %>% mutate(label = paste0(top, "-", side))
colnames(p.gsts) <- mat.name.index$label

#results
st.dev <- apply(p.gsts, 2, sd)
# mean <- apply(p.gsts, 2, mean)
ul <- rbind(obs) + st.dev*1.96
colnames(ul) <- mat.name.index$label
ul <- as.data.frame(ul) %>%
  pivot_longer(cols = everything(), names_to = "pair", values_to = "upper.limit") %>%
  mutate(top = sub("-[0-9]{4}_[0-9]{1}", "", pair),
         side = sub("[0-9]{4}_[0-9]{1}-", "", pair)) %>%
  select(-pair)

ll <- rbind(obs) - st.dev*1.96
colnames(ll) <- mat.name.index$label
ll <- as.data.frame(ll) %>%
  pivot_longer(cols = everything(), names_to = "pair", values_to = "lower.limit") %>%
  mutate(top = sub("-[0-9]{4}_[0-9]{1}", "", pair),
         side = sub("[0-9]{4}_[0-9]{1}-", "", pair)) %>%
  select(-pair)

ci <- left_join(ul, ll, by=c("top","side"))

write_csv(ci, paste0(PATH, "/results_tables/pGst_CI.csv"))
cat("\nCI saved.\n")

## check outputs
ci <- read_csv(paste0(PATH, "/results_tables/pGst_CI.csv"))
ci <- ci %>% arrange(top, side)

ll <- ci %>% 
  pivot_wider(names_from = top, values_from = lower.limit) %>%
  select(-upper.limit) %>%
  group_by(side) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = T))) %>%
  mutate(across(where(is.numeric), ~na_if(., -Inf))) 
ul <- ci %>% 
  pivot_wider(names_from = top, values_from = upper.limit) %>%
  select(-lower.limit) %>%
  group_by(side) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = T))) %>%
  mutate(across(where(is.numeric), ~na_if(., -Inf))) 
write_csv(ll, paste0(PATH, "/results_tables/pGst_lowerlimit.csv"))
write_csv(ul, paste0(PATH, "/results_tables/pGst_upperlimit.csv"))

p.gsts <- read_csv(paste0(PATH, "/results_tables/pGst_boots_saved.csv"))
