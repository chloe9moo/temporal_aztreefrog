##TITLE: Source code for loading diff. groupings of the genetic data
##AUTHOR: C. E. Moore
##Updated on 15 SEP 2022

#load in data
aztf <- read.csv(paste0(PATH, "/microsat_data/final_aztf_loci_nosibs.csv"))
aztf <- aztf %>% mutate(year = ifelse(year == 2018, 2019, year)) %>%
  mutate(year = ifelse(year == 2015, 2014, year))
coord <- read.csv(paste0(PATH, "/pond_coordinates.csv"))
pop.coord <- coord %>% dplyr::select(pop, lon, lat) %>% unique() %>% remove_rownames() %>% column_to_rownames("pop")

ind_aztf_pop <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_pop@pop <- as.factor(aztf$pop) #assign populations
ind_aztf_pop$other$xy <- subset(pop.coord, rownames(pop.coord) %in% aztf$pop) #assign coordinates

ind_aztf_yr <- df2genind(aztf[,c(6:22)], sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
ind_aztf_yr@pop <- as.factor(aztf$year)
ind_aztf_yr$other$xy <- data.frame(lon = c(2014, 2019, 2021), lat = c(1,1,1))

o.2015 <- aztf %>% filter(year == '2014')
ind2015 <- df2genind(o.2015[,c(6:22)], sep = "/", ind.names = o.2015$ID, NA.char = "-1/-1")
ind2015@pop <- as.factor(o.2015$pop)
ind2015$other$xy <- subset(pop.coord, rownames(pop.coord) %in% o.2015$pop)

o.2019 <- aztf %>% filter(year == '2019')
ind2019 <- df2genind(o.2019[,c(6:22)], sep = "/", ind.names = o.2019$ID, NA.char = "-1/-1")
ind2019@pop <- as.factor(o.2019$pop)
ind2019$other$xy <- subset(pop.coord, rownames(pop.coord) %in% o.2019$pop)

o.2021 <- aztf %>% filter(year == '2021')
ind2021 <- df2genind(o.2021[,c(6:22)], sep = "/", ind.names = o.2021$ID, NA.char = "-1/-1")
ind2021@pop <- as.factor(o.2021$pop)
ind2021$other$xy <- subset(pop.coord, rownames(pop.coord) %in% o.2021$pop)

#make pop data
pop_aztf_pop <- genind2genpop(ind_aztf_pop, quiet = T) #global
pop_aztf_yr <- genind2genpop(ind_aztf_yr, quiet = T)
pop2015 <- genind2genpop(ind2015, quiet = T)
pop2019 <- genind2genpop(ind2019, quiet = T)
pop2021 <- genind2genpop(ind2021, quiet = T)

#make list
ind.list <- list(ind_aztf_pop, ind_aztf_yr, ind2015, ind2019, ind2021)
ind.names <- list("All Yrs - Indiv x Pop", "All Yrs - Indiv x Year", "Indiv 2014", "Indiv 2019", "Indiv 2021")
pop.list <- list(pop_aztf_pop, pop_aztf_yr, pop2015, pop2019, pop2021)
pop.names <- list("All Yrs - Pop x Pop", "All Yrs - Pop x Year", "Pop 2014", "Pop 2019", "Pop 2021")

#color palette
aztf.pal <- c("1" = "#28231D", 
              "3" = "#5E2D30", 
              "4" = "#008E90", 
              "6" = "#5cc589", 
              "7" = "#1C77A3", 
              "8" = "#C5A387",
              "9" = "#67B8D6", 
              "10" = "#E9D097", 
              "16" = "#283F79")

#remove extras
rm(ind_aztf_pop, ind_aztf_yr, ind2015, ind2019, ind2021, 
   pop_aztf_pop, pop_aztf_yr, pop2015, pop2019, pop2021, 
   o.2015, o.2019, o.2021,
   pop.coord)
