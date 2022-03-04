## Sampling sites
### World map
Given it is a "global diversity" study, worth having a world map with sampling sites, distinction between ancient and modern samples, and the fact that some some from humans, animals, and the environment (ancient).

```R
setwd("/nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/MAP")

# load libraries
library(ggplot2)
library(dplyr)
require(maps)
library(ggrepel)
library(patchwork)
library(ggsci)

# load world data
world_map <- map_data("world")

# load metadata
data <- read.delim("map_metadata.txt", sep="\t", header=T)

# make a map
ggplot() +
     geom_polygon(data = world_map, aes(x = world_map$long, y = world_map$lat, group = world_map$group), fill="grey90") +
     geom_point(data = data, aes(x = LONGITUDE, y = LATITUDE, colour = REGION, shape = SAMPLE_AGE), size=3) +
     geom_text_repel(data = data, aes(x = LONGITUDE, y = LATITUDE, label = paste0(COUNTRY," (",POPULATION_ID,"); n = ", SAMPLE_N)), size=3, max.overlaps = Inf) +
     theme_void() +
     ylim(-55,85) +
     labs(title="A", colour="", shape="") +
     scale_colour_npg()

# save it
ggsave("worldmap_samplingsites.png", height=5, width=12)
ggsave("worldmap_samplingsites.pdf", height=5, width=12, useDingbats=FALSE)

```

- will use this as Figure 1A
![worldmap_samplingsites](../04_analysis/worldmap_samplingsites.png)

```R

library(ggplot2)
library(dplyr)
require(maps)
library(ggrepel)
library(patchwork)
library(ggsci)

# load world data
world_map <- map_data("world")


# code for the scale bar - need to cut and paste some functions from here: https://egallic.fr/en/scale-bar-and-north-arrow-on-a-ggplot2-map/

data <- read.delim("map_metadata_ancient.txt", sep="\t", header=T)
fst <- read.delim("map_fst_data.txt", sep="\t", header=T)

ggplot() +
     geom_polygon(data = world_map, aes(x = world_map$long, y = world_map$lat, group = world_map$group), fill="grey90",col="white") +
     geom_segment(data=fst, aes(x=POP1_LONG , y=POP1_LAT ,xend=POP2_LONG , yend=POP2_LAT, size=1-FST ),linetype="dashed")+
     geom_point(data = data, aes(x = LONGITUDE, y = LATITUDE, colour = COUNTRY), size=3) +
     geom_text_repel(data = data, aes(x = LONGITUDE, y = LATITUDE, label = POPULATION_ID), size=3)+
     coord_cartesian(xlim = c(0,28), ylim = c(51,58))+
     scale_colour_npg()+
     scale_size_continuous(range = c(0, 1), guide = 'none')+
     theme_bw()+ labs(x="Longitude", y="Latitude") #+
     #scale_bar(lon = 25, lat = 51,
     #distance_lon = 100, distance_lat = 30, distance_legend = 54,
     #dist_unit = "km", orientation = FALSE)
```

x
### Sampling timepoints

```R
library(ggplot2)

data <- read.delim("ancient_times.txt",header=F,sep="\t")

ggplot(data, aes(x=V11,xend=V12,y=reorder(paste0(V1," (",V4,")"),V11,FUN=mean),yend=paste0(V1," (",V4,")"), colour=V10)) +
     geom_segment(size=5) +
     xlim(1000,2020) +
     labs(x = "Estimated age of sampling site (AD)", y = "", colour = "Sample site") +
     scale_y_discrete(limits=rev) +
     theme_bw() + theme(legend.position="bottom")

ggsave("samplingsites_time.png", height=5, width=7)
ggsave("samplingsites_time.pdf", height=5, width=7, useDingbats=FALSE)

```

Figure: [map](../04_analysis/samplingsites_time.pdf)
- will use this in the supplementary data

![samplingsites_time](../04_analysis/samplingsites_time.png)
