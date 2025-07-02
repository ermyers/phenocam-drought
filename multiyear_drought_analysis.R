# Analysis of multi-year droughts

library(sf)
library(terra)
library(tidyverse)
library(patchwork)
library(lubridate)

# Load continuous drought and ecoregion data
load("outputs/drought_statistics.RData")
load("outputs/phen_with_ecoregion.RData")

# Filter drought statistics to only include primary ROIs
unique_phenocam_rois_annotated <- read.csv("data/unique_phenocam_rois_annotated_04-30-2025.csv")
primary_rois <- filter(unique_phenocam_rois_annotated, action=="keep" & rank==1)
drought_statistics_filtered <- filter(drought_statistics, Phenocam %in% primary_rois$Phenocam)
rm(drought_statistics,phen_usdm_modified,unique_phenocam_rois_annotated)

# Add primary veg type and lat/lon to drought statistics
drought_statistics_filtered <- left_join(drought_statistics_filtered, select(primary_rois, Phenocam, Veg_Type, Primary_Veg_Type, lat, lon), by=join_by(Phenocam))

# Add ecoregion and aridity index to drought statistics
drought_statistics_filtered <- left_join(drought_statistics_filtered, select(phen_with_ecoregion_df, site, NA_L1CODE, NA_L1NAME, aridity_index), by=join_by(Phenocam==site))
rm(phen_with_ecoregion_df)

drought_statistics_filtered <- drought_statistics_filtered %>%
  mutate(aridity_label = case_when(aridity_index<0.03 ~ "Hyper Arid",
                                   aridity_index<0.2 ~ "Arid",
                                   aridity_index<0.5 ~ "Semi-Arid",
                                   aridity_index<0.65 ~ "Dry Sub-Humid",
                                   aridity_index>=0.65 ~ "Humid")) %>%
  mutate(aridity_label = factor(aridity_label, levels=c("Hyper Arid", "Arid", "Semi-Arid", "Dry Sub-Humid", "Humid")))

# Select multiyear droughts
multiyear_droughts <- filter(drought_statistics_filtered, Number_Weeks>=104)

# Multiyear drought statistics
veg_type_count <- multiyear_droughts %>% count(Veg_Type)
primary_veg_type_count <- multiyear_droughts %>% count(Primary_Veg_Type)
aridity_index_count <- multiyear_droughts %>% count(aridity_label)
ecoregion_count <- multiyear_droughts %>% count(NA_L1NAME)

#######################
# Visual representation
#######################

# Plot number of weeks vs. severity
p1 <- ggplot() +
  geom_point(data = drought_statistics_filtered, aes(x=Sum_USDM, y=Number_Weeks, color=aridity_label)) +
  labs(color = "Aridity Index") +
  xlab("Cumulative DM") +
  ylab("Drought Length (weeks)") +
  ggtitle("Continuous droughts at PhenoCam locations") +
  theme_bw()
p1

# Load vector data (USA states and phenocam locations)
#type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
                                     usa_states_vect$NAME != "Hawaii" &
                                     usa_states_vect$NAME != "Puerto Rico",]

# Convert phenocam shapefile into dataframe with coordinates
#type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to drought statistics
#drought_statistics_df <- left_join(drought_statistics_filtered, select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert site-year and drought statistics into vectors
drought_statistics_vect <- vect(drought_statistics_filtered, geom=c("lon", "lat"), crs="epsg:4326", keepgeom=FALSE)

# Only look at sites in the contiguous United States
drought_statistics_vect <- terra::intersect(drought_statistics_vect,project(usa_states_vect,crs(drought_statistics_vect)))

# Select droughts longer than 104 weeks
multiyear_droughts_vect <- drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]

# # Plot locations of droughts longer than 104 weeks (color = number of weeks)
# p2 <- ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]), aes(color=Number_Weeks)) +
#   labs(color = "Number of Weeks in Drought") +
#   scale_color_viridis_c(option = 'viridis') +
#   ggtitle("Longest droughts (>=104 weeks) at phenocam locations")
# p2
# 
# # Plot locations of droughts longer than 104 weeks (color = aridity index)
# p3 <- ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]), aes(color=aridity_label)) +
#   labs(color = "Aridity Index") +
#   #scale_color_viridis_c(option = 'viridis') +
#   ggtitle("Longest droughts (>=104 weeks) at phenocam locations")
# p3
# 
# # Plot locations of droughts longer than 104 weeks (color = number of weeks, shape = aridity label)
# p4 <- ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]), aes(color=Number_Weeks, shape=aridity_label), size=2.5) +
#   labs(color = "Number of Weeks in Drought", shape = "Aridity Index") +
#   scale_color_viridis_c(option = 'viridis') +
#   ggtitle("Longest droughts (>=104 weeks) at phenocam locations")
# p4

# Plot locations of droughts longer than 104 weeks (color = number of weeks, shape = aridity label)
cols <- RColorBrewer::brewer.pal(9, "YlOrRd")
cols <- cols[4:9]
breaks <- seq(0,200, by=20)
p4 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]), aes(color=Number_Weeks, shape=aridity_label), size=2.5) +
  labs(color = "Number of Weeks in Drought", shape = "Aridity Index") +
  scale_color_stepsn(colours = cols,
                    breaks = seq(0,200, by=20),
                    name = "Number of Weeks in Drought") +
  ggtitle("Longest droughts (>=104 weeks)") +
  theme_bw()
p4

# # Plot droughts with >150 cumulative severity
# p5 <- ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Sum_USDM>=150,]), aes(color=Sum_USDM)) +
#   labs(color = "Cumulative Drought Severity") +
#   scale_color_viridis_c(option = 'viridis') +
#   ggtitle("Most severe droughts (>=150 cumulative DM) at phenocam locations")
# p5

# Plot droughts with >150 cumulative severity (color = severity, shape = aridity label)
cols <- RColorBrewer::brewer.pal(9, "YlOrRd")
cols <- cols[4:9]
p6 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Sum_USDM>=150,]), aes(color=Sum_USDM, shape=aridity_label), size=2.5) +
  labs(color = "Cumulative Drought Severity", shape = "Aridity Index") +
  scale_color_stepsn(colours = cols,
                     breaks = seq(0,500, by=100),
                     name = "Cumulative Drought Severity") +
  ggtitle("Most severe droughts (>=150 cumulative DM)") +
  theme_bw()
p6

# Current size for copying drought plots: 700x300
