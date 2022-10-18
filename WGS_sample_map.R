###making a map of the samples used in the WGS project

### loady loady
# set options
options(stringsAsFactors = F)         # no automatic data transformation
options("scipen" = 100, "digits" = 4) # suppress math annotation
op <- options(gvis.plot.tag='chart')  # set gViz options

# load packages
library(OpenStreetMap)
library(DT)
library(RColorBrewer)
library(mapproj)
library(sf)
library(RgoogleMaps)
library(scales)
library(rworldmap)
library(maps)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(maptools)
library(leaflet)
library(sf)
library(tmap)
library(here)
library(rgdal)
library(scales)
library(flextable)
library(ggplot2)
library(cowplot)
library(rcartocolor)
library(ggmap)
# activate klippy for copy-to-clipboard button
klippy::klippy()



# set working directory here
setwd("/Users/libbynatola/Documents/UBC/Bioinformatics/wgs/map")

#load coordinate/sample data
wgs_samples <- read.csv("../wgs_sample_info.csv")


# make world map                  
world <- ne_countries(scale = "medium", returnclass = "sf")

# plot it
pdf("wgs_sample_map.pdf")
ggplot() +
  geom_sf(data = world, fill="gray97") + 
  labs( x = "Longitude", y = "Latitude", color = "Transect") +
  coord_sf(xlim = c(-155,-67), ylim = c(65,35), expand = T) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(legend.position=c(0.1, 0.3)) +
  geom_point(shape = 20,alpha = 0.7, size = 3.5, data=wgs_samples, aes(x=lon, y=lat, color=Group)) + scale_colour_manual(values = c("#bd7777", "#bd5757", "#b26ec4", "#394ca2", "#2f5734", "grey", "black", "#e6e045", "#c6924d"), labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), name=NULL) + labs(x="Longitude", y="Latitude")
dev.off()

pdf("wgs_sample_map_shapes.pdf")
ggplot() +
  geom_sf(data = world, fill="gray97") + 
  labs( x = "Longitude", y = "Latitude", color = "Transect") +
  coord_sf(xlim = c(-155,-67), ylim = c(65,35), expand = T) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(legend.position=c(0.1, 0.3)) +
  geom_point(alpha = 0.7, size = 3.5, data=wgs_samples, aes(x=lon, y=lat, color=Group, shape = Group)) + scale_colour_manual(name = "Species", values = c("#bd7777", "#bd5757", "#b26ec4", "#394ca2", "#2f5734", "grey25", "black", "#e6e045", "#c6924d"), labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber")) + 
  scale_shape_manual(name = "Species", labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), values = c(16, 19, 5, 15, 6, 13, 18, 17, 0)) + labs(x="Longitude", y="Latitude")
dev.off()

#scale_shape_manual(values = c(16, 16, 5, 15, 6, 4, 13, 17, 5)
#set coordinate limits to include everyone
base = get_map(location=c(-155,33,-60,67), zoom=5, maptype="terrain")

#I makea the map
map1 = ggmap(base)
map1

# add samples
map2 <- map1 + geom_point(shape = 1, data=wgs_samples, aes(x=lon, y=lat, color=Group, pch = 1))
map2

# logical colors and axes labels
map3 <- map2 + scale_colour_manual(values = c("red4", "red", "purple", "blue", "green", "grey", "black", "yellow", "orange")) + labs(x="Longitude", y="Latitude") 
map3

# latin names in scale
map4 <- map2 + scale_colour_manual(values = c("red4", "red", "purple", "blue", "green", "grey", "black", "yellow", "orange"), labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), name=NULL) + labs(x="Longitude", y="Latitude")
map4

# pdf(file = "WGS_sample_map.pdf")
map4
#dev.off()


# make world map                  
world <- ne_countries(scale = "medium", returnclass = "sf")

# jitter lat and long
wgs_samples$jitter_lat <- jitter(wgs_samples$lat, factor = 500)
wgs_samples$jitter_lon <- jitter(wgs_samples$lon, factor = 500)

# plot it
pdf("wgs_map.pdf")
ggplot() +
  geom_sf(data = world, fill="gray98") + 
  labs( x = "Longitude", y = "Latitude", color = "Transect") +
  coord_sf(xlim = c(-147,-67), ylim = c(36,66), expand = T) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(text = element_text(size = 18)) +
  geom_point(size = 4, data=wgs_samples, aes(x=jitter_lon, y=jitter_lat, color = Group, shape = Group)) + scale_colour_manual(name = "Species", labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), values = c("red4", "red", "purple", "blue", "green", "grey40", "black", "gold", "orange")) +
  scale_shape_manual(name = "Species", labels = c("Red-breasted daggetti", "Red-breasted ruber", "Red-naped x Red-breasted ruber", "Red-naped", "Red-naped x Yellow-bellied", "Red-breasted x Red-naped x Yellow-bellied", "Williamson's", "Yellow-bellied", "Yellow-bellied x Red-breasted ruber"), values = c(6, 2, 14, 0, 10, 8, 3, 1, 13))
dev.off()

