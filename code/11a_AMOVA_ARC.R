#AMOVA ON ARC
library(tidyverse); library(adegenet); library(poppr)

PATH <- "/data"

#load in data
aztf <- read.csv(paste0(PATH, "/final_aztf_loci_nosibs.csv"))
aztf <- aztf %>% mutate(year = ifelse(year == 2018, 2019, year)) %>%
  mutate(year = ifelse(year == 2015, 2014, year))
coord <- read.csv(paste0(PATH, "/pond_coordinates.csv")) %>% arrange(pop)

ind_aztf_pop <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_pop@pop <- as.factor(aztf$pop) #assign populations
ind_aztf_pop$other$xy <- subset(pop.coord, rownames(pop.coord) %in% aztf$pop) #assign coordinates

ind.list[[1]]@other$year <- aztf$year 
ind.list[[1]]@other$pop <- aztf$pop
strata(ind.list[[1]]) <- data.frame(other(ind.list[[1]])$year, other(ind.list[[1]])$pop) #year will be first group, pop is sub group
nameStrata(ind.list[[1]]) <- ~year/pop #rename

#amova
amova <- poppr::poppr.amova(ind.list[[1]], hier = ~year/pop, method = "ade4")

#test for significance
amova.test <- ade4::randtest(amova, nrepet = 10000)

#save
saveRDS(amova.test, file = paste0(PATH, "/AMOVA_All_Test.rds"))


# #check
# test <- readRDS(paste0(PATH, "/AMOVA_All_Test.rds"))
