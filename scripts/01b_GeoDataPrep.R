##TITLE: Prepping spatial data for downstream analyses
##AUTHOR: C. E. Moore
##Updated on 15 SEP 2022

library(tidyverse); library(sf); library(mapview); library(adegenet)

# set up ----
PATH <- "/home/chloe9mo/Documents/Projects/temporal_aztreefrog"
source(file = paste0(PATH, "/code/01c_DataLoad.R"))

dat <- read.csv(paste0(PATH, "/multiyear_hywr_tracking.csv")) %>% 
  select("Individual.ID", "Type","Year.Sampled","Pond.Number", "Pond.Name","UTME","UTMN") %>%
  filter(Individual.ID != "19-015")

dat.sf <- st_as_sf(dat, coords = c("UTME","UTMN"), crs=32612)
pop.dat <- dat.sf %>% select("Year.Sampled","Pond.Number") %>% unique() %>% 
  mutate(Year.Sampled = fct_relevel(as.factor(Year.Sampled), "2014", "2018", "2019", "2021"))

# pal <-  mapviewPalette("mapviewSpectralColors")

mapview(pop.dat, zcol="Year.Sampled", legend=T)

#adding pop numbers
# dat[dat$Pond.Name == "Riata Tank", c("Pond.Number")] <- "17"
# dat[dat$Pond.Name == "Turkey Creek Crossing (83)", c("Pond.Number")] <- "1"
# dat[dat$Pond.Name == "Anderson Dam", c("Pond.Number")] <- "18"
# dat[dat$Pond.Name == "Perch Tank", c("Pond.Number")] <- "19"
# dat[dat$Pond.Name == "Parker Tank", c("Pond.Number")] <- "20"
# dat[dat$Pond.Name == "Mountain Tank", c("Pond.Number")] <- "21"
# dat[dat$Pond.Name == "Garden Tank", c("Pond.Number")] <- "22"
# dat[dat$Pond.Name == "Midden Tank", c("Pond.Number")] <- "23"
# dat[dat$Pond.Name == "Cemetery Tank", c("Pond.Number")] <- "24"
# dat[dat$Pond.Name == "Double Tank (east)", c("Pond.Number")] <- "10"
# dat[dat$Pond.Name == "Duncan Tank", c("Pond.Number")] <- "25"


#summarizing ----
sumdat <- read.csv(paste0(PATH, "/microsat_data/final_aztf_loci_nosibs.csv")) %>% select(year, pop) %>%
  mutate(year = case_when(year == 2018 ~ 2019,
                          year == 2015 ~ 2014,
                          T ~ as.numeric(year)),
         pop = case_when(pop == "1" ~ 1,
                         pop == "4" ~ 2,
                         pop == "3" ~ 3,
                         pop == "6" ~ 4,
                         pop == "7" ~ 5,
                         pop == "8" ~ 6,
                         pop == "9" ~ 7,
                         pop == "10" ~ 8,
                         pop == "16" ~ 9)) %>%
  count(year, pop) %>%
  pivot_wider(id_cols=pop, names_from = year, values_from=n) %>%
  replace(is.na(.), 0) %>%
  mutate(TotSamp = rowSums(across(.cols = c('2014','2019','2021'))))
write.csv(sumdat, file = paste0(PATH,"/HYWR_multiyear_popsumm_nosib.csv"), row.names = F)

# dat <- dat %>% mutate(UTME = round(as.numeric(UTME), 0), UTMN = round(as.numeric(UTMN), 0)) %>%
#   arrange()
# 
# pop.key <- dat %>% select(Pond.Number, Year.Sampled, UTME, UTMN) %>% unique()
# pop.key <- pop.key %>% 
#   arrange(Pond.Number, desc(Year.Sampled == '2014')) %>% 
#   distinct(Pond.Number, .keep_all=T) %>%
#   select(Pond.Number, UTME, UTMN)
# 
# dat.fix <- dat %>% select(-UTME, -UTMN) %>% left_join(., pop.key, by="Pond.Number")
# 
# write.csv(dat.fix, file=paste0(PATH, "/updated_multiyear_hywr_tracksheet.csv"), row.names = F)

#plot new coords
# dat.fix.sf <- st_as_sf(dat.fix, coords = c("UTME","UTMN"), crs=32612)
# dat.fix.sf %>% select(Pond.Number) %>% unique() %>% mapview(., legend=F)

#pond coords ----
# xy <- dat %>% dplyr::select(Pond.Number, UTME, UTMN) %>% unique()
# xy <- st_as_sf(xy, coords = c("UTME", "UTMN"), crs = 32612) #crs is WGS 84, UTM zone 12N
# #checking it worked
# # mapview(xy)
# #convert to lat long
# xy.ll <- st_transform(xy, crs = "+proj=longlat +datum=WGS84")
# #checking it worked
# # mapview(list(xy, xy.ll), col.regions=list("red", "blue"))
# #save as dataframe
# xy.ll <- xy.ll %>% 
#   mutate(lon = sf::st_coordinates(.)[,1],
#          lat = sf::st_coordinates(.)[,2]) %>% 
#   st_drop_geometry()
# loci <- read.csv(paste0(PATH, "/microsat_data/combined_aztf_loci.csv"))
# loci <- loci %>% mutate(pop = as.character(pop))
# 
# xy.ll <- loci %>% select(pop, year_pop) %>% left_join(., xy.ll, by = c("pop" = "Pond.Number")) %>% unique()
# write.csv(xy.ll, file=paste0(PATH, "/pond_coordinates.csv"), row.names = F)

#sample extent maps ----
# library(ggmap); library(cowplot)

sumdat <- read.csv(paste0(PATH,"/HYWR_multiyear_popsumm_nosib.csv")) %>% rename_with(.cols = everything(), ~ gsub("X", "", .x))

coord <- coord %>%
  filter(pop %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")) %>%
  left_join(sumdat, by="pop") %>%
  mutate(pop = fct_relevel(as.factor(pop), c("1", "2", "3", "4", "5", "6", "7", "8", "9"))) %>%
  rename(y2014 = '2014',
         y2019 = '2019',
         y2021 = '2021') %>%
  st_as_sf(., coords=c("lon", "lat"), crs=st_crs(4269), remove=F)

e <- as.data.frame(as.matrix(st_bbox(coord %>% st_transform(26912)))) %>% rownames_to_column(var = "name") %>% pivot_wider(names_from = name, values_from = V1)
e[,1] <- e[,1] - 2000 #honestly not sure what the units are
e[,2] <- e[,2] - 10000 
e[,3] <- e[,3] + 2000 
e[,4] <- e[,4] + 2000 

#+ clip background layers ----
library(terra); library(tidyterra); library(ggnewscale)

#get hillshade
elev <- rast(paste0(PATH, "/ElevationData/Combined_ElevationRast.tif"))
slope <- terrain(elev, "slope", unit = "radians")
aspect <- terrain(elev, "aspect", unit = "radians")
hill <- shade(slope, aspect, 25, 270)
#get extent and crop
e.r <- ext(c(xmin=e$xmin, xmax=e$xmax, ymin=e$ymin, ymax=e$ymax))
hill.c <- crop(hill, e.r)

pal_greys <- hcl.colors(1500, "Grays", rev = T)
pal_greys <- pal_greys[c(1:1000)]
# pal_greys <- pal_greys[c(100:600,900:1400)]

#test background
# ggplot() +
#   geom_spatraster(data=hill.c, show.legend = T) +
#   geom_sf(data = coord) +
#   scale_fill_gradientn(colors = pal_greys, na.value = NA)

#+ sample extent map for manuscript ----
library(maps); library(scatterpie); library(ggrepel); library(ggspatial)

#full N. America inset
na <- st_as_sf(map("world", regions = c("usa(?!:Hawaii)(?!:Alaska)", "mexico"), plot = FALSE, fill = TRUE)) %>%
  st_transform(st_crs("+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"))
#az shape
az <- st_as_sf(map("state", "arizona", plot = FALSE, fill = TRUE)) %>% st_transform(st_crs("+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"))
#aztf range polygon
range <- read_sf(paste0(PATH, "/spatial_data/range_map/Approx_Range_Map.shp")) %>% st_transform(st_crs("+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"))

#usa + mex inset map
p.full <- ggplot() +
  geom_sf(data = na, fill="grey", color="black", size=0.3) +
  geom_sf(data = az, fill="#444444", color="black") +
  geom_sf(data = range, fill="red", color="red", alpha = 0.5, size=0.2) +
  geom_sf(data = coord %>% filter(pop == 5), fill=NA, color="yellow", shape=22, size=5, alpha=0.6, stroke=1) +
  labs(x="Longitude", y="Latitude") +
  # coord_sf(xlim=c(-2500000, 0), ylim = c(-1800000, 1200000)) + #see more of the usa/mex
  coord_sf(xlim=c(-2100000, -800000), ylim = c(-1700000, 100)) +
  theme(rect = element_rect(fill="transparent"),
        panel.grid.major = element_line(color="gray", linetype = "dashed"),
        panel.background = element_rect(fill="transparent"),
        #axis.text = element_blank(),
        panel.border = element_rect(color = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12))
p.full
ggsave(filename=paste0(PATH, "/figures/Full_StudyExtent_Map.png"), plot=p.full, width = 4, height = 6, bg="transparent")

#reformat data for plotting
coord.base <- coord %>% select(-year, -year_pop) %>% distinct(.) %>%
  mutate(radius = TotSamp/7000, 
         across(.cols = starts_with("l"), ~as.numeric(.x)),
         across(.cols = starts_with("y20"), ~as.numeric(.x)))
pie.coord <- coord.base %>% st_drop_geometry()

#match projection
hill.c <- project(hill.c, crs(coord.base))

#set colors
az.colors <- c("#F4D1AE", "#12263A", "#00C6B8")

#plot
p.zoom <- ggplot() +
  geom_spatraster(data = hill.c, show.legend = F) +
  scale_fill_gradientn(colors = pal_greys, na.value = NA) +
  new_scale_fill() +
  geom_text_repel(data = coord.base %>% filter(TotSamp < 37), aes(x=lon, y=lat, label=pop), nudge_y = 0.009, nudge_x = 0.009, size=7, fontface="bold") +
  geom_text_repel(data = coord.base %>% filter(TotSamp < 80 & TotSamp > 50), aes(x=lon, y=lat, label=pop), nudge_y = 0.013, nudge_x = 0.013, size=7, fontface="bold") +
  geom_text_repel(data = coord.base %>% filter(TotSamp > 80), aes(x=lon, y=lat, label=pop), nudge_y = 0.015, nudge_x = 0.014, size=7, fontface="bold") +
  geom_scatterpie(data = pie.coord, aes(x=lon, y=lat, r=radius, fill=type), cols = c("y2014", "y2019", "y2021"), size = 0.2) +
  # coord_equal() +
  # geom_scatterpie_legend(coord.base$radius, x=-110.508, y=31.382, labeller = function(x) x*7000, n=3) +
  # geom_text(x=-110.51, y=31.405, label="Arizona treefrog \nindividuals collected", size=5, hjust = 0.5) +
  scale_fill_manual(values = az.colors, aesthetics = "fill",
                    guide = guide_legend(title = "", title.position = "top", direction = "vertical"),
                    labels = c('2014', '2019', '2021')) +
  labs(x="Longitude", y="Latitude") +
  annotation_scale(location = "bl", width_hint = 0.24, pad_y = unit(0.17, "inches"), pad_x = unit(0.2, "inches"), style = "ticks", text_cex = 1.1, text_pad = unit(0.3, "cm"), tick_height = 1, line_width = 1.5) +
  annotation_north_arrow(location = "tr", which_north = "grid", width = unit(0.35, "inches"), height = unit(0.35, "inches")) +
  scale_x_continuous(expand = c(-0.001,-0.001)) +
  scale_y_continuous(expand = c(-0.001,-0.001), position = "right") +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill=NA, color = "black"),
        legend.position = c(-0.15, -0.016),
        # legend.position = "left",
        legend.justification = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(1.5, "line"),
        legend.background = element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 15),
        plot.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.y = element_text(angle = 90, hjust = 0.6),
        # axis.title = element_text(size = 15),
        # axis.title.y = element_text(hjust = 0.1),
        axis.title = element_blank())
# p.zoom
ggsave(filename = paste0(PATH, "/figures/Full_StudyExtent_Inset.png"), plot = p.zoom, dpi = 300, width = 8.5, height = 5.5, bg="transparent")

transp <- hcl.colors(1000, "viridis", alpha = 0)
ggplot() + 
  geom_spatraster(data=hill.c, show.legend = F) +
  scale_fill_gradientn(colors = transp, na.value = NA) +
  geom_scatterpie_legend(coord.base$radius, x=-110.4, y=31.45, labeller = function(x) x*7000, n=3) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid.major = element_blank())
  # geom_scatterpie(data = pie.coord, aes(x=lon, y=lat, r=radius, fill=type), cols = c("y2014", "y2019", "y2021"), size = 0.2)
ggsave(paste0(PATH, "/figures/Full_StudyExtent_PieLegend.png"), plot = last_plot(), width = 8, height = 5, bg="transparent")

#+ pop number map ----
coord <- coord %>% st_as_sf(coords=c("lon", "lat"), crs=st_crs(32612))
pm <- ggplot() +
  geom_spatraster(data=hill.c, show.legend = F) +
  scale_fill_gradientn(colors = pal_greys, na.value = NA) +
  new_scale_fill() +
  # geom_sf(data = coord, size=2, shape=4, inherit.aes = F) +
  # geom_sf(data = coord %>% filter(pop %in% c("1", "3", "4", "6", "7", "8", "9", "10", "16")), aes(fill=pop), color="black", shape=21, inherit.aes = F) +
  geom_sf_label(data = coord , aes(fill=pop, label=pop), color="white", size=7, inherit.aes = F) +
  coord_sf(datum = sf::st_crs(coord)) +
  scale_color_manual(values = aztf.pal, aesthetics = "color") +
  scale_fill_manual(values = aztf.pal, aesthetics = "fill") +
  scale_x_continuous(expand = c(-0.001,-0.001)) +
  scale_y_continuous(expand = c(-0.001,-0.001)) +
  labs(x = "", y = "", title = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA, color = "black"),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 18),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))
pm
ggsave(filename = paste0(PATH,"/figures/dapc_aztf_samplemap.png"), plot = pm, width = 7, height = 7)

# #+ pop size maps ----
# nsampmap.graphics <- list(
#   coord_sf(datum = sf::st_crs(coord)),
#     # scale_color_manual(values = aztf.pal, aesthetics = "color"),
#     scale_fill_manual(values = aztf.pal, aesthetics = "fill", guide = "none", labels = coord$pop),
#     scale_size_continuous(range = c(4, 14), guide = "legend",
#                           limits=c(0, 45),
#                           breaks = c(5, 15, 25, 35, 45)),
#     labs(x = "Longitude", y = "Latitude", size="# of \nindividuals"),
#     theme(panel.background = element_blank(),
#           panel.border = element_rect(fill=NA, color = "black"),
#           legend.key = element_blank(),
#           legend.text = element_text(size=18),
#           legend.title = element_text(size = 18),
#           plot.title = element_text(face = "bold", size = 18),
#           axis.title = element_text(size=15),
#           axis.text = element_text(size=13))
# )
# 
# pm2 <- ggmap(az.base) +
#   geom_sf(data = coord %>% filter(y2014 == 0), size=7, shape=13, inherit.aes = F) +
#   geom_sf(data = coord %>% filter(y2014 > 0), aes(fill=pop, size=y2014), color="black", shape=21, inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("4", "6", "7", "8", "9", "10")), aes(label=pop), color="black", inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("1", "3")), aes(label=pop), color="white", inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("16")), aes(label=pop), color="black", nudge_y = 0.008, inherit.aes = F) +
#   ggtitle("2014 sampled populations") +
#   nsampmap.graphics +
#   theme(legend.position = "bottom",
#         legend.title.align = 1) +
#   guides(size = guide_legend(label.position = "bottom"))
# pm2
# pm3 <- ggmap(az.base) +
#   geom_sf(data = coord %>% filter(y2019 == 0), size=7, shape=13, inherit.aes = F) +
#   geom_sf(data = coord %>% filter(y2019 > 0), aes(fill=pop, size=y2019), color="black", shape=21, inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(y2019 > 0), aes(label=pop), color="black", inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(y2019 == 0), aes(label=pop), color="black", nudge_y = 0.008, inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("16")), aes(label=pop), color="white", inherit.aes = F) +
#   ggtitle("2019 sampled populations") +
#   nsampmap.graphics +
#   theme(legend.position = "bottom", 
#         axis.title.y = element_blank(), 
#         axis.text.y = element_blank(),
#         legend.title.align = 1) +
#   guides(size = guide_legend(label.position = "bottom"))
# pm3
# pm4 <- ggmap(az.base) +
#   geom_sf(data = coord %>% filter(y2021 == 0), size=7, shape=13, inherit.aes = F) +
#   geom_sf(data = coord %>% filter(y2021 > 0), aes(fill=pop, size=y2021), color="black", shape=21, inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("4", "6", "7", "8", "9")), aes(label=pop), color="black", inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("1", "3", "16")), aes(label=pop), color="white", inherit.aes = F) +
#   # geom_sf_text(data = coord %>% filter(pop %in% c("10")), aes(label=pop), color="black", nudge_y = 0.008, inherit.aes = F) +
#   ggtitle("2021 sampled populations") +
#   nsampmap.graphics +
#   theme(legend.position = "bottom", 
#         axis.title.y = element_blank(), 
#         axis.text.y = element_blank(),
#         legend.title.align = 1) +
#   guides(size = guide_legend(label.position = "bottom"))
# pm4
# 
# pm5 <- ggmap(az.base) +
#   geom_sf(data = coord, aes(fill=pop, size=TotSamp), color="black", shape=21, inherit.aes = F) +
#   geom_sf_text(data = coord %>% filter(pop %in% c("4", "6", "7", "8", "9", "1")), aes(label=pop), color="black", inherit.aes = F) +
#   geom_sf_text(data = coord %>% filter(pop %in% c("10", "3", "16")), aes(label=pop), color="white", inherit.aes = F) +
#   ggtitle("Total Samples") +
#   scale_size_continuous(range = c(2, 12), guide = "legend",
#                         limits=c(10, 105),
#                         breaks = c(20, 40, 60, 80, 100)) +
#   coord_sf(datum = sf::st_crs(coord)) +
#   scale_fill_manual(values = aztf.pal, aesthetics = "fill", guide = "none", labels = coord$pop) +
#   labs(x = "Longitude", y = "Latitude", size="# of \nindividuals") +
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(fill=NA, color = "black"),
#         legend.key = element_blank(),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size = 15),
#         plot.title = element_text(face = "bold", size = 18),
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=13))
#   
# pm5
# 
# ggsave(filename = paste0(PATH,"/figures/aztf_samplemap_2014.png"), plot = pm2, width = 6.8, height = 6.5)
# ggsave(filename = paste0(PATH,"/figures/aztf_samplemap_2019.png"), plot = pm3, width = 6.3, height = 6.5)
# ggsave(filename = paste0(PATH,"/figures/aztf_samplemap_2021.png"), plot = pm4, width = 6.3, height = 6.5)
# ggsave(filename = paste0(PATH,"/figures/aztf_samplemap_total.png"), plot = pm5, width = 7, height = 7)


#kmz to point for network creation ----
sp.coord <- st_as_sf(coord, coords = c("lon", "lat"), crs = 4326) #pond coordinates w/ any samples
net <- lapply(list.files(paste0(PATH, "/spatial_data"), pattern = "Huachuca", full.names = T), st_read) %>% #huachuca ponds id'ed by traci
  bind_rows() %>%
  st_zm() %>%
  st_transform(st_crs(sp.coord))

rd.net <- st_filter(net, sp.coord %>%
                      st_buffer(7000) %>% #7 km like parsley et al 2020
                      st_union() %>%
                      st_make_valid()) %>%
  select(Name, geometry)
st_write(rd.net, paste0(PATH, "/spatial_data/huachuca_network_ponds.shp"), delete_dsn = T)

rm(net)

#stream dist between sites ----
library(riverdist)

#decided to go into qgis and select the important flowlines instead
#the full nhd is ~20 G of mem btw
# nhd <- st_read("~/Documents/Projects/PNW_fishes/NHD/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb", layer = "NHDFlowline_Network")
# nhd <- st_filter(nhd, sp.coord %>%
#                    st_buffer(100000) %>% #7 km like parsley et al 2020
#                    st_union() %>%
#                    st_make_valid() %>%
#                    st_transform(st_crs(nhd))) %>% 
#   st_zm()
# st_write(nhd, paste0(PATH, "/spatial_data/huachuca_nhdflowline.shp"), delete_dsn=T)

##load in for package
nhd <- line2network(path=paste0(PATH, "/spatial_data"), layer="aztf_streamnetwork",
                    reproject = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
nhd_fixed <- cleanup(nhd)
saveRDS(nhd_fixed, file=paste0(PATH, "/spatial_data/aztf_streamnetwork_fixed.rds"))

rd.pt.c <- sp.coord %>% 
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  mutate(lon = st_coordinates(.)[,"X"],
         lat = st_coordinates(.)[,"Y"]) %>%
  select(-year_pop) %>%
  filter(pop %in% aztf$pop) %>%
  filter(!pop %in% c(7, 9)) %>% #these are closer to a different river network, probably need to find distance between them tho
  distinct()
rd.pt <- xy2segvert(x=rd.pt.c$lon, y=rd.pt.c$lat, rivers = nhd_fixed)
rd.pt$id <- rd.pt.c$pop

#plot to check
points(rd.pt.c$lon, rd.pt.c$lat, pch=16, col="red")
riverpoints(seg=rd.pt$seg, vert=rd.pt$vert, rivers=nhd_fixed, pch=15, col="blue")

#dist matrix
dmat <- riverdistancemat(rd.pt$seg, rd.pt$vert, nhd_fixed, ID=rd.pt$id)
write.csv(dmat, file=paste0(PATH, "/results_tables/RiverDist_Matrix.csv"), row.names = T)

#add back in other pops
md <- max(dmat) * 10 #one order of magnitude greater than furthest stream distance for pops 7+9 and all others
nhd_fixed2 <- trimriver(trimto = c(283, 284, 285), rivers=nhd) #only the relevant segments
nhd_fixed2 <- cleanup(nhd_fixed2)
rd.pt.c <- sp.coord %>% 
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  mutate(lon = st_coordinates(.)[,"X"],
         lat = st_coordinates(.)[,"Y"]) %>%
  select(-year_pop) %>%
  filter(pop %in% c(7, 9)) %>% #these are closer to a different river network, probably need to find distance between them tho
  distinct()
rd.pt <- xy2segvert(x=rd.pt.c$lon, y=rd.pt.c$lat, rivers = nhd_fixed2)
rd.pt$id <- rd.pt.c$pop
d2mat <- riverdistancemat(rd.pt$seg, rd.pt$vert, nhd_fixed, ID=rd.pt$id) %>% as.data.frame()


dmat <- read.csv(paste0(PATH, "/results_tables/RiverDist_Matrix.csv"))
d2mat[names(dmat)[names(dmat) != 'X']] <- md
d2mat <- d2mat %>% 
  rename_with(.cols = -contains("X"), ~ paste0("X", .)) %>%
  rownames_to_column(var="pop")
dmat <- dmat %>%
  rename(pop = X) %>%
  mutate('X7' = md, 'X9' = md, pop = as.character(pop))

dmat <- full_join(dmat, d2mat)
write.csv(dmat, file=paste0(PATH, "/results_tables/RiverDist_Matrix.csv"), row.names = T)

#prep precip data ----
library(terra); library(exactextractr)
sp.coord <- st_as_sf(coord, coords = c("lon", "lat"), crs = 4326) #pond coordinates w/ any samples
rd.net <- st_read(paste0(PATH, "/spatial_data/huachuca_network_ponds.shp")) #read in pond network

pri.l <- list.files(path = paste0(PATH, "/spatial_data/prism"), pattern = "*.bil$", full.names = T)
c.buff <- sp.coord %>% 
  st_transform(st_crs(rast(pri.l[[1]]))) %>%
  st_buffer(7000)
s.buff <- sp.coord %>%
  st_transform(st_crs(rast(pri.l[[1]]))) %>%
  st_buffer(7000) %>%
  st_union()

for (i in 1:length(pri.l)) {
  
  pri <- rast(pri.l[[i]]) #get layer
  yr <- paste0("ppt_", gsub("PRISM_ppt_stable_4kmM3_", "", names(pri)) %>% gsub("_bil", "", .)) #set year
  ra <- exact_extract(pri, c.buff, "mean", weights = "area") #get average for each ind. site
  c.buff <- c.buff %>% mutate(precip = ra) #add to df
  names(c.buff)[names(c.buff) == "precip"] <- yr #set col name
  
  ra <- exact_extract(pri, s.buff, "mean", weights = "area") #get avg for all sites
  sd <- exact_extract(pri, s.buff, "stdev", weights = "area") #get sd across all sites
  
  c.buff <- c.buff %>% mutate(precip = ra, stdev = sd) #add to df
  names(c.buff)[names(c.buff) == "precip"] <- paste0(yr, "_all") #set col name
  names(c.buff)[names(c.buff) == "stdev"] <- paste0(yr, "_all_sd") #set col name
  
}

c.buff <- c.buff %>% 
  rowwise() %>%
  mutate(d21 = sd(c(ppt_2021, ppt_2020, ppt_2019)),
        d19 = sd(c(ppt_2019, ppt_2018, ppt_2017)),
        d14 = sd(c(ppt_2014, ppt_2013, ppt_2012))) %>% 
  st_drop_geometry()
write.csv(c.buff, paste0(PATH, "/spatial_data/Precip4Sites.csv"), row.names = F)

#look at precip change over time
precip.change <- data.frame(year = seq(2012, 2021, 1), ppt.a = NA, ppt.sd = NA)

for (i in 1:nrow(precip.change)) {
  yr <- precip.change[i,1]
  tryCatch({ #if an error is thrown when trying to select the column, fill that cell with na
    
    c.buff %>% select(paste0("ppt_", yr, "_all"))
    precip.change[i,2] <- (c.buff %>% select(paste0("ppt_", yr, "_all")) %>% st_drop_geometry() %>% distinct())[1,1]
  
    }, error=function(e) precip.change[i,2] <- NA)
  tryCatch({ #if an error is thrown when trying to select the column, fill that cell with na
    
    c.buff %>% select(paste0("ppt_", yr, "_all_sd"))
    precip.change[i,3] <- (c.buff %>% select(paste0("ppt_", yr, "_all_sd")) %>% st_drop_geometry() %>% distinct())[1,1]
    
  }, error=function(e) precip.change[i,3] <- NA)
  
}

ggplot(data=na.omit(precip.change), aes(x=year, y=ppt.a)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=ppt.a-ppt.sd, ymax=ppt.a+ppt.sd), size=0.1) + theme_classic()
