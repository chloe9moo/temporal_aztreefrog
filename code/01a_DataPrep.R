##TITLE: AZ treefrog: Microsat data prep and 2015 comparisons
##AUTHOR: C. E. Moore
##Updated on SEP 29 2022

#-------------------------------#
#### Load in data, set paths ####
#-------------------------------#
library(tidyverse); library(readxl); library(pegas)

PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"
PATH_15 <- paste0(PATH, "/microsat_data/2015")
PATH_18 <- paste0(PATH, "/microsat_data/2018")
PATH_19 <- paste0(PATH, "/microsat_data/2019")
PATH_21 <- paste0(PATH, "/microsat_data/2021")

aztf19 <- lapply(list.files(path = PATH_19, pattern = "2019_Mix*", full.names = TRUE), read.csv)
aztf15 <- lapply(list.files(path = PATH_15, pattern = "*Alleles.csv", full.names = TRUE), read.csv)
aztf18 <- lapply(list.files(path = PATH_18, pattern = "*Alleles.csv", full.names = T), read.csv)
aztf21 <- lapply(list.files(path = PATH_21, pattern = "*Alleles.csv", full.names = T), read.csv)

mims <- read_excel(paste0(PATH_15, "/Mims2016_GenotypeData.xlsx"), sheet = 2)

## remove Mix# / set# from name
aztf19 <- lapply(aztf19, function (x) { cbind(indiv = gsub("([0-9]+)_.*", "\\1", x$Name), x) })
aztf19 <- lapply(aztf19, function (x) { x[!(names(x) %in% c("Name"))] })

aztf15 <- lapply(aztf15, function (x) { cbind(indiv = gsub("([0-9]+)_.*", "\\1", x$Name), x) })
aztf15 <- lapply(aztf15, function (x) { x[!(names(x) %in% c("Name"))] })

aztf21 <- lapply(aztf21, function (x) { cbind(indiv = gsub("([0-9]+)_.*", "\\1", x$Name), x) })
aztf21 <- lapply(aztf21, function (x) { x[!(names(x) %in% c("Name"))] })

aztf18 <- lapply(aztf18, function (x) { cbind(indiv = gsub("([0-9]+)_.*", "\\1", x$Name), x) })
aztf18 <- lapply(aztf18, function (x) { x[!(names(x) %in% c("Name"))] })

#-----------------------------#
#### 2015 data comparisons ####
#-----------------------------#

#remove weird individual duplicates -- already did this in geneious
# aztf15[[1]] <- aztf15[[1]][-c(24,26,28,30,32,34,36,38,39,40,43,44,46),]
# aztf15[[2]] <- aztf15[[2]][-c(6,32,34,35,37,44),]
# aztf15[[3]] <- aztf15[[3]][-c(26,28,30,32,34,36,38,40,41,42,44,46),]

#reformat and consolidate
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("4269")) %>% unique()
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("29495")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("1422")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("30215")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("3318")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[1]] %>% select(indiv, contains("2932")) %>% unique() %>% right_join(., aztf2015)

aztf2015 <- aztf15[[2]] %>% select(indiv, contains("1316")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[2]] %>% select(indiv, contains("16672")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[2]] %>% select(indiv, contains("34484")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[2]] %>% select(indiv, contains("2688")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[2]] %>% select(indiv, contains("4093")) %>% unique() %>% right_join(., aztf2015)

aztf2015 <- aztf15[[3]] %>% select(indiv, contains("23452")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[3]] %>% select(indiv, contains("20812")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[3]] %>% select(indiv, contains("4370")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[3]] %>% select(indiv, contains("12115")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[3]] %>% select(indiv, contains("30594")) %>% unique() %>% right_join(., aztf2015)
aztf2015 <- aztf15[[3]] %>% select(indiv, contains("10374")) %>% unique() %>% right_join(., aztf2015)

#subtract original data from new calls
mims <- mims[c(1:59), c(1,4:37)] #remove rows and cols not needed
mims <- mims %>% rename(indiv = Sample)
mims <- mims[order(mims$indiv),]

aztf2015 <- aztf2015[order(aztf2015$indiv),]
colnames(aztf2015) <- gsub("\\.\\.\\.", "_", colnames(aztf2015))
colnames(aztf2015) <- gsub("Hwri", "", colnames(aztf2015))
indiv <- aztf2015[,1]
aztf2015 <- aztf2015[,-1]
mims <- mims[,-1]

cols <- sort(intersect(names(aztf2015), names(mims)))
call.diff <- aztf2015[cols] - mims[cols]
call.diff$indiv <- indiv
call.diff <- call.diff[,c(35, 1:34)]
mims$indiv <- indiv
mims <- mims[,c(35,1:34)]
aztf2015$indiv <- indiv
aztf2015 <- aztf2015[,c(35,1:34)]
call.diff[is.na(call.diff)] <- 0

write.csv(call.diff, file = paste0(PATH, "/subset_differences2015.csv"), row.names = FALSE)

# // summarize differences ----

colnames(call.diff) <- paste0("hwri", colnames(call.diff))
call.diff <- call.diff %>% rename(indiv = hwriindiv)
sum.diff <- call.diff %>% replace(is.na(.), 0) %>%
  transmute(indiv = indiv,
            hwri10374 = hwri10374_1+hwri10374_2,
            hwri12115 = hwri12115_1+hwri12115_2,
            hwri1316 = hwri1316_1+hwri1316_2,
            hwri1422 = hwri1422_1+hwri1422_2,
            hwri16672 = hwri16672_1+hwri16672_2,
            hwri20812 = hwri20812_1+hwri20812_2,
            hwri23452 = hwri23452_1+hwri23452_2,
            hwri2688 = hwri2688_1+hwri2688_2,
            hwri2932 = hwri2932_1+hwri2932_2,
            hwri29495 = hwri29495_1+hwri29495_2,
            hwri30215 = hwri30215_1+hwri30215_2,
            hwri30594 = hwri30594_1+hwri30594_2,
            hwri3318 = hwri3318_1+hwri3318_2,
            hwri34484 = hwri34484_1+hwri34484_2,
            hwri4093 = hwri4093_1+hwri4093_2,
            hwri4269 = hwri4269_1+hwri4269_2,
            hwri4370 = hwri4370_1+hwri4370_2)
write.csv(sum.diff, file = paste0(PATH, "/sum_diff_2015.csv"), row.names = FALSE)

summary <- sum.diff %>% mutate_if(is.numeric, ~1 * (. != 0))
summary <- colMeans(summary[sapply(summary, is.numeric)])
write.csv(summary, file = paste0(PATH, "/avg_incorrect_loci.csv"))


# // adjust + check adjustments ----

colnames(mims) <- paste0("hwri", colnames(mims))
mims.2 <- mims %>%
  mutate(hwri16672_1 = if_else(hwri16672_1 == 177, 177, hwri16672_1 - 4), #keep bin 177 the same, subtract 4 from all others
         hwri16672_2 = if_else(hwri16672_2 == 177, 177, hwri16672_2 - 4),
         hwri20812_1 = case_when(hwri20812_1 == 296 ~ 296, hwri20812_1 == 372 ~ 380, TRUE ~ as.numeric(hwri20812_1)+4), #bin 296 no change, bin 372 + 8, all others + 4
         hwri20812_2 = case_when(hwri20812_2 == 296 ~ 296, hwri20812_2 == 372 ~ 380, TRUE ~ as.numeric(hwri20812_2)+4),
         hwri23452_1 = hwri23452_1 + 2, #add 2 to all indivs
         hwri23452_2 = hwri23452_2 + 2,
         hwri34484_1 = case_when(hwri34484_1 == 163 ~ 167, hwri34484_1 == 155 ~ 159, TRUE ~ as.numeric(hwri34484_1)), #163 becomes 167, 155 becomes 159
         hwri34484_2 = case_when(hwri34484_2 == 163 ~ 167, hwri34484_2 == 155 ~ 159, TRUE ~ as.numeric(hwri34484_2)),
         hwri4370_1 = hwri4370_1 - 2, #subtract 2 from all indivs
         hwri4370_2 = hwri4370_2 - 2) 
colnames(mims.2) <- gsub("hwri", "", colnames(mims.2))
indiv <- mims.2[,1]
mims.2 <- mims.2[,-1]
aztf2015 <- aztf2015[,-1]

cols <- sort(intersect(names(aztf2015), names(mims.2)))
call.diff.2 <- aztf2015[cols] - mims.2[cols]
call.diff.2$indiv <- indiv
call.diff.2 <- call.diff.2[,c(35, 1:34)]
mims.2$indiv <- indiv
mims.2 <- mims.2[,c(35,1:34)]
aztf2015$indiv <- indiv
aztf2015 <- aztf2015[,c(35,1:34)]
call.diff.2[is.na(call.diff.2)] <- 0

write.csv(call.diff.2, file = paste0(PATH, "/subset_differences2015_updated.csv"), row.names = FALSE)

colnames(call.diff.2) <- paste0("hwri", colnames(call.diff.2))
call.diff.2 <- call.diff.2 %>% rename(indiv = hwriindiv)
sum.diff <- call.diff.2 %>% replace(is.na(.), 0) %>%
  transmute(indiv = indiv,
            hwri10374 = hwri10374_1+hwri10374_2,
            hwri12115 = hwri12115_1+hwri12115_2,
            hwri1316 = hwri1316_1+hwri1316_2,
            hwri1422 = hwri1422_1+hwri1422_2,
            hwri16672 = hwri16672_1+hwri16672_2,
            hwri20812 = hwri20812_1+hwri20812_2,
            hwri23452 = hwri23452_1+hwri23452_2,
            hwri2688 = hwri2688_1+hwri2688_2,
            hwri2932 = hwri2932_1+hwri2932_2,
            hwri29495 = hwri29495_1+hwri29495_2,
            hwri30215 = hwri30215_1+hwri30215_2,
            hwri30594 = hwri30594_1+hwri30594_2,
            hwri3318 = hwri3318_1+hwri3318_2,
            hwri34484 = hwri34484_1+hwri34484_2,
            hwri4093 = hwri4093_1+hwri4093_2,
            hwri4269 = hwri4269_1+hwri4269_2,
            hwri4370 = hwri4370_1+hwri4370_2)
write.csv(sum.diff, file = paste0(PATH, "/sum_diff_2015_updated.csv"), row.names = FALSE)

summary <- sum.diff %>% mutate_if(is.numeric, ~1 * (. != 0))
summary <- colMeans(summary[sapply(summary, is.numeric)])
write.csv(summary, file = paste0(PATH, "/avg_incorrect_loci_updated.csv"))

# // update 2015 data and save ----

aztf2015.v2 <- read_excel(paste0(PATH_15, "/Mims2016_GenotypeData.xlsx"), sheet = 2)
colnames(aztf2015.v2) <- paste0("hwri", colnames(aztf2015.v2))
aztf2015.v2 <- aztf2015.v2 %>% rename(indiv = hwriSample, pop = hwriSiteNum, age = 'hwriLife Stage')
aztf2015.v2 <- aztf2015.v2 %>%
  mutate(hwri16672_1 = if_else(hwri16672_1 == 177, 177, hwri16672_1 - 4), #keep bin 177 the same, subtract 4 from all others
         hwri16672_2 = if_else(hwri16672_2 == 177, 177, hwri16672_2 - 4),
         hwri20812_1 = case_when(hwri20812_1 == 296 ~ 296, hwri20812_1 == 372 ~ 380, TRUE ~ as.numeric(hwri20812_1)+4), #bin 296 no change, bin 372 + 8, all others + 4
         hwri20812_2 = case_when(hwri20812_2 == 296 ~ 296, hwri20812_2 == 372 ~ 380, TRUE ~ as.numeric(hwri20812_2)+4),
         hwri23452_1 = hwri23452_1 + 2, #add 2 to all indivs
         hwri23452_2 = hwri23452_2 + 2,
         hwri34484_1 = case_when(hwri34484_1 == 163 ~ 167, hwri34484_1 == 155 ~ 159, TRUE ~ as.numeric(hwri34484_1)), #163 becomes 167, 155 becomes 159
         hwri34484_2 = case_when(hwri34484_2 == 163 ~ 167, hwri34484_2 == 155 ~ 159, TRUE ~ as.numeric(hwri34484_2)),
         hwri4370_1 = hwri4370_1 - 2, #subtract 2 from all indivs
         hwri4370_2 = hwri4370_2 - 2)
aztf2015.v2$year <- 2015
aztf2015.v2 <- aztf2015.v2[,c(1:3,38,4:37)]
aztf2015.v2$pop <- as.character(aztf2015.v2$pop)

write.csv(aztf2015.v2, file = paste0(PATH_15, "/Mims2016_GenotypeData_Updated.csv"))

#--------------------------#
#### reformat data 2019 ####
#--------------------------#

#remove negatives
aztf19 <- lapply(aztf19, function (x) { x[!grepl("neg", x$indiv),] })
ml <- lapply(aztf19, function (x) colnames(x)) %>% #get loci names
  unlist() %>% as.data.frame() %>% rename(loci = '.') %>% #create dataframe
  mutate(loci = gsub("Hwri(.+)\\.\\.\\.[1-2]", "\\1", loci)) %>% unique() %>% #get only the numbers of the loci
  filter(loci != 'indiv')

#set dataframe to populate with microsat data
aztf2019 <- lapply(aztf19, function (x) x %>% select(indiv)) %>% bind_rows() %>% unique()

#consolidate alleles into single dataframe
for (i in 1:3) {
  for(y in 1:nrow(ml)) {
    
    l <- ml[y,1]
    aztf2019 <- aztf19[[i]] %>% select(indiv, contains(l)) %>% unique() %>% na.omit() %>% right_join(., aztf2019)
    
  }
}

rr <- aztf19[[4]]
aztf2019c <- right_join(rr, aztf2019) #doesn't seem to be any differences in the rerun plate

#add in population info
pops <- read.csv(paste0(PATH, "/multiyear_hywr_tracking.csv"))
pops <- pops %>% select(Individual.ID, Pond.Number, Type) %>% rename(indiv = Individual.ID, pop = Pond.Number, age = Type)
pops[pops[, c("indiv")] == "19-074*",]$indiv <- "19-074"
pops <- pops %>% 
  mutate(age = case_when(age == 'Buccal' ~ 'Adult', age == 'Tail' ~ 'Larvae', age == 'Toe' ~ 'Adult',
                         age == 'tail' ~ 'Larvae', age == 'toe' ~ 'Adult',  age == 'toe ' ~ 'Adult', TRUE ~ age),
         # pop = case_when(pop == 'T9' ~ '7',
         #                 pop == 'T8' ~ '6',
         #                 pop == 'T11' ~ '9',
         #                 pop == 'FS4636' ~ '4',
         #                 pop == 'UGCP' ~ '8',
         #                 pop == 'T4' ~ '11',
         #                 pop == 'A5' ~ '12',
         #                 pop == 'Unknown' ~ '15',
         #                 pop == 'Double Tank Up-Canelo Hills' ~ '14',
         #                 pop == 'T10' ~ '13',
         #                 T ~ pop)
         ) %>%
  filter(age %in% c('Adult', 'Larvae'))

aztf2019 <- merge(aztf2019, pops, by = "indiv")
aztf2019$year <- 2019
aztf2019 <- aztf2019[,c(1,36,37,38,2:35)]
colnames(aztf2019) <- sub("\\.\\.\\.", "_", colnames(aztf2019)) #matching names to 2015 data
names(aztf2019) <- tolower(names(aztf2019))
aztf2019[is.na(aztf2019)] <- '-1' #match nas to 2015

write.csv(aztf2019, file = paste0(PATH_19, "/AZTF19_GenotypeData.csv"), row.names = FALSE)

#--------------------------#
#### reformat data 2021 ####
#--------------------------#

#remove negatives
aztf21 <- lapply(aztf21, function (x) { x[!grepl("neg", x$indiv),] })
ml <- lapply(aztf21, function (x) colnames(x)) %>% #get loci names
  unlist() %>% as.data.frame() %>% rename(loci = '.') %>% #create dataframe
  mutate(loci = gsub("Hwri(.+)\\.\\.\\.[1-2]", "\\1", loci)) %>% unique() %>% #get only the numbers of the loci
  filter(loci != 'indiv')

#set dataframe to populate with microsat data
aztf2021 <- lapply(aztf21, function (x) x %>% select(indiv)) %>% bind_rows() %>% unique()

#consolidate alleles into single dataframe
for (i in 1:3) {
  for(y in 1:nrow(ml)) {
    
    l <- ml[y,1]
    aztf2021 <- aztf21[[i]] %>% select(indiv, contains(l)) %>% unique() %>% na.omit() %>% right_join(., aztf2021)
    
  }
}


#add in population info
# pops <- read.csv(paste0(PATH, "/updated_multiyear_hywr_tracksheet.csv"))
# pops <- pops %>% 
#   filter(Year.Sampled == "2021") %>%
#   mutate(indiv = str_extract(Individual.ID, "21_.*")) %>%
#   rename(year = Year.Sampled, pop = Pond.Number, age = Type) %>%
#   dplyr::select(indiv, age, pop, year) %>%
#   mutate(indiv = str_replace(indiv, "_", "-"),
#          age = case_when(age == 'toe' ~ 'Adult', age == 'toe ' ~ 'Adult', age == 'tail' ~ 'Larvae', TRUE ~ age)) %>%
#   filter(indiv %in% aztf2021$indiv)
aztf2021 <- merge(aztf2021, pops, by = "indiv")
aztf2021$year <- 2021
aztf2021 <- aztf2021[,c(1,36,37,38,2:35)]

colnames(aztf2021) <- sub("\\.\\.\\.", "_", colnames(aztf2021)) #matching names to 2015 data
names(aztf2021) <- tolower(names(aztf2021))
aztf2021[is.na(aztf2021)] <- '-1' #match nas to 2015

write.csv(aztf2021, file = paste0(PATH_21, "/AZTF21_GenotypeData.csv"), row.names = FALSE)

#--------------------------#
#### reformat data 2018 ####
#--------------------------#

#remove negatives
aztf18 <- lapply(aztf18, function (x) { x[!grepl("neg", x$indiv),] })
ml <- lapply(aztf18, function (x) colnames(x)) %>% #get loci names
  unlist() %>% as.data.frame() %>% rename(loci = '.') %>% #create dataframe
  mutate(loci = gsub("Hwri(.+)\\.\\.\\.[1-2]", "\\1", loci)) %>% unique() %>% #get only the numbers of the loci
  filter(loci != 'indiv')

#set dataframe to populate with microsat data
aztf2018 <- lapply(aztf18, function (x) x %>% select(indiv)) %>% bind_rows() %>% unique()

#consolidate alleles into single dataframe
for (i in 1:3) {
  for(y in 1:nrow(ml)) {
    
    l <- ml[y,1]
    aztf2018 <- aztf18[[i]] %>% select(indiv, contains(l)) %>% unique() %>% na.omit() %>% right_join(., aztf2018)
    
  }
}


#add in population info
# pops <- read.csv(paste0(PATH, "/updated_multiyear_hywr_tracksheet.csv"))
# pops <- pops %>% 
#   filter(Year.Sampled == "2018") %>%
#   mutate(indiv = str_extract(Individual.ID, "18_.*")) %>%
#   rename(year = Year.Sampled, pop = Pond.Number, age = Type) %>%
#   dplyr::select(indiv, age, pop, year) %>%
#   mutate(indiv = str_replace(indiv, "_", "-"),
#          age = case_when(age == 'toe' ~ 'Adult', age == 'toe ' ~ 'Adult', age == 'tail' ~ 'Larvae', TRUE ~ age)) %>%
#   filter(indiv %in% aztf2018$indiv)
aztf2018 <- merge(aztf2018, pops, by = "indiv")
aztf2018$year <- 2018
aztf2018 <- aztf2018[,c(1,36,37,38,2:35)]

colnames(aztf2018) <- sub("\\.\\.\\.", "_", colnames(aztf2018)) #matching names to 2015 data
names(aztf2018) <- tolower(names(aztf2018))
aztf2018[is.na(aztf2018)] <- '-1' #match nas to 2015

write.csv(aztf2018, file = paste0(PATH_18, "/AZTF18_GenotypeData.csv"), row.names = FALSE)


#-----------------------------------------#
#### combine 2015, 2018, 2019, 2021 - format ####
#-----------------------------------------#
#read in data
aztf2019 <- read.csv(paste0(PATH_19, "/AZTF19_GenotypeData.csv")) %>% mutate(pop = as.character(pop))
aztf2015.v2 <- read.csv(paste0(PATH_15, "/Mims2016_GenotypeData_Updated.csv")) %>% select(-X) %>%
  mutate(old.pop = as.character(pop)) %>%
  mutate(pop = case_when(old.pop == '6'~'9', old.pop == '9'~'6', T ~ old.pop)) %>%
  select(-old.pop)
aztf2015.v2 <- aztf2015.v2 %>% 
  mutate(across(starts_with("hwri"), ~ ifelse(.x < 10, -1, .x))) #change the incorrect neg. values
aztf2021 <- read.csv(paste0(PATH_21, "/AZTF21_GenotypeData.csv")) %>% mutate(pop = as.character(pop))
aztf2018 <- read.csv(paste0(PATH_18, "/AZTF18_GenotypeData.csv")) %>% mutate(pop = as.character(pop))
# aztf2019[,c(5:38)] <- apply(aztf2019[,c(5:38)], 2, as.numeric)

aztf_loci <- aztf2015.v2 %>% full_join(aztf2019) %>% full_join(aztf2021) %>% full_join(aztf2018)

write.csv(aztf_loci, file = paste0(PATH, "/microsat_data/combined_aztf_genotypes.csv"), row.names = FALSE)

#format two column allele format to single column with / separating start and end
##note that alleles2loci seems to incorrectly order alleles when the first one is <100 bp
microsat <- aztf_loci[,c(5:38)] #separate out loci
microsat <- pegas::alleles2loci(microsat) #convert
colnames(microsat) <- sub("_1", "", colnames(microsat))
microsat[c("ID", "pop", "age", "year")] <- aztf_loci[,c(1:4)]
microsat <- microsat %>% mutate(year_pop = paste0(year, "_", pop)) #add in another grouping option
microsat[microsat == "2018_16"] <- "2019_16"
microsat <- microsat[,c(18:22,1:17)]

write.csv(microsat, file=paste0(PATH, "/microsat_data/combined_aztf_loci.csv"), row.names = F)

# remove siblings ----
aztf <- read.csv(paste0(PATH, "/microsat_data/combined_aztf_loci.csv"))
sibs <- read.csv(paste0(PATH, "/multiyear_hywr_tracking.csv"))

sibs <- sibs %>% select(Individual.ID, Sibling.Action)
k.sibs <- sibs %>% filter(Sibling.Action == "keep") %>% select(Individual.ID) %>% rename(ID = Individual.ID)

aztf.s <- aztf %>% filter(ID %in% k.sibs$ID)

write.csv(aztf.s, file=paste0(PATH, "/microsat_data/combined_aztf_loci_nosibs.csv"), row.names = F)

#check which sibs weren't removed that Meryl removed
# mims <- read_excel(paste0(PATH_15, "/Mims2016_GenotypeData.xlsx"), sheet = 3)
# aztf.f <- read.csv(paste0(PATH, "/microsat_data/filtered_aztf_loci.csv"))
# mims <- mims %>% select(Sample) %>% rename(ID = Sample)
# aztf.s %>% filter(year == '2015' & !ID %in% mims$ID) %>% select(ID) #samples I kept, Meryl didn't
# mims %>% filter(!ID %in% aztf.s$ID)

# NeEstimator prep ----
library(adegenet)
# aztf <- read.csv(paste0(PATH,"/microsat_data/final_aztf_loci_sibs.csv"))
aztf <- read.csv(paste0(PATH,"/microsat_data/final_aztf_loci_nosibs.csv"))
microsat <- aztf %>% select(contains("hwri"))
# microsat <- df2genind(microsat, sep = "/", ind.names = aztf$ID, NA.char = "-1/-1")
microsat <- microsat %>% 
  mutate(across(everything(), ~gsub("/", "", .x))) %>% 
  mutate(across(everything(), ~gsub("-1-1", "000000", .x)))
microsat["pop"] <- aztf$year_pop
microsat <- microsat[,c(18,1:17)]
microsat <- arrange(microsat, pop)

# write.table(microsat, file=paste0(PATH,"/Ne_Estimator/Ne_loci_original_sibs.txt"), sep = '\t', row.names = F)
write.table(microsat, file=paste0(PATH,"/Ne_Estimator/Ne_loci_original_nosibs.txt"), sep = '\t', row.names = F)


# STRUCTURE prep ----
aztf <- read.csv(paste0(PATH,"/microsat_data/final_aztf_loci_nosibs.csv"))
pop.key <- aztf %>% select(year_pop) %>% distinct() %>% mutate(group = row_number())
write.csv(pop.key, file = paste0(PATH, "/microsat_data/year_pop_number_key4structure.csv"), row.names = F)

str <- aztf %>% select(-pop, -year, -age) %>% left_join(., pop.key) %>% 
  mutate(Label = row_number()) %>%
  relocate(ID, Label, group) %>% 
  select(-year_pop)
str %>% select(ID, Label) %>% write.csv(., file=paste0(PATH, "/microsat_data/ind_key4structure.csv"), row.names = F)
str <- str %>% select(-ID)

ind_aztf_pop <- df2genind(str[,c(3:19)], sep = "/", ind.names = str$Label, NA.char = "-1/-1")
ind_aztf_pop@pop <- as.factor(str$group) #assign populations
str <- genind2df(ind_aztf_pop, sep = "/", usepop = T, oneColPerAll = T)
str <- str %>% mutate(Label = row.names(str), Pop = pop) %>% select(-pop) %>% relocate(Label, Pop)
# str <- str %>% select(-Label)
names(str) <- gsub("\\.1", "", names(str))
names(str) <- gsub("\\.2", "", names(str))
str[str=="NA"] <- -1

write.table(str, file=paste0(PATH,"/microsat_data/structure_loci.txt"), sep = '\t', row.names = F, quote = F)




