#AMOVA ON ARC
.libPaths(.libPaths()[3:1])
library(tidyverse); library(adegenet); library(poppr)

PATH <- "/data"
PATH <- "~/Documents/Projects/temporal_aztreefrog"

#load in data
aztf <- read.csv(paste0(PATH, "/final_aztf_loci_nosibs.csv"))
aztf <- aztf %>% mutate(year = ifelse(year == 2018, 2019, year)) %>%
  mutate(year = ifelse(year == 2015, 2014, year))
coord <- read.csv(paste0(PATH, "/pond_coordinates.csv")) %>% arrange(pop)
pop.coord <- coord %>% dplyr::select(pop, lon, lat) %>% unique() %>% remove_rownames() %>% column_to_rownames("pop")

ind_aztf_pop <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_pop@pop <- as.factor(aztf$pop) #assign populations
ind_aztf_pop$other$xy <- subset(pop.coord, rownames(pop.coord) %in% aztf$pop) #assign coordinates

ind_aztf_pop@other$year <- aztf$year
ind_aztf_pop@other$pop <- aztf$pop
strata(ind_aztf_pop) <- data.frame(other(ind_aztf_pop)$year, other(ind_aztf_pop)$pop) #year will be first group, pop is sub group
nameStrata(ind_aztf_pop) <- ~year/pop #rename
           
#amova
amova <- poppr::poppr.amova(ind_aztf_pop, hier = ~year/pop, method = "ade4")

#test for significance
amova.test <- ade4::randtest(amova, nrepet = 10000)

#save
saveRDS(amova.test, file = paste0(PATH, "/AMOVA_All_10K.rds"))


#check
test <- readRDS(paste0(PATH, "/AMOVA_All_10K.rds"))
test
