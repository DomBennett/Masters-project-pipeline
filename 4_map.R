## MRes Project 2013
## Stage 4: 
## In: 0_data, 1_results | Out: 2_figures
## 21/08/2013

## Libraries
source(file.path('functions','EcoDataTools.R'))
require(ggplot2)
require(maps)
require(mapdata)

## Dirs
data.dir <- "0_data"
figures.dir <- "2_figures"
results.dir <- "1_results"

## 0 Input + data manipulation
all.data <- read.csv(file.path(data.dir, "urban_data.csv"), stringsAsFactors = FALSE)
pdmpd <- read.csv(file.path(results.dir, "pdmpd.csv"))

## 1. Make a map!
studies <-paste0(all.data$Source_ID, all.data$Study_number)
ustudies <- unique(studies)
long <- all.data$Longitude[match(ustudies, studies)]
lat <- all.data$Latitude[match(ustudies, studies)]
phyla <- all.data$Phylum[match(ustudies, studies)]
present <- as.numeric(ustudies %in% unique(pdmpd$studies))
mapobj <- data.frame(study = ustudies, lat, long, Phylum = phyla)
mapobj <- mapobj[as.logical(present),]
worldmap <- map_data(map = "worldHires", ylim=c(-60,90))
worldmap <- geom_path(aes(x = long, y = lat, group = group),
                      data = worldmap)
map <- ggplot() + worldmap +
  geom_point(aes(x = long, y = lat, colour = Phylum), size = 3,
             alpha = 0.8, data = mapobj) +
  theme_bw() + xlab("Latitude") + ylab("Longitude")
pdf(file.path(figures.dir, "map.pdf"), height = 6, width = 10.5)
print (map)
dev.off()