# Code for pulling RAP data at phenocam locations

library(sf)
# remotes::install_github(repo = "landscape-data-commons/trex")
library(trex)
library(geojsonsf)
library(terra)
library(tidyverse)
library(patchwork)
library(lubridate)

# Load data
load("outputs/drought_statistics.RData")
unique_phenocam_rois_annotated <- read.csv("data/unique_phenocam_rois_annotated_04-30-2025.csv")
primary_rois <- filter(unique_phenocam_rois_annotated, action=="keep" & rank==1)
drought_statistics_filtered <- filter(drought_statistics, Phenocam %in% primary_rois$Phenocam)
rm(drought_statistics,phen_usdm_modified,unique_phenocam_rois_annotated)

type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
                                     usa_states_vect$NAME != "Hawaii" &
                                     usa_states_vect$NAME != "Puerto Rico",]

# Convert phenocam shapefile into dataframe with coordinates
type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to drought statistics
drought_statistics_df <- left_join(drought_statistics_filtered, select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert site-year and drought statistics into vectors
drought_statistics_vect <- vect(drought_statistics_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Only look at sites in the contiguous United States
drought_statistics_vect <- terra::intersect(drought_statistics_vect,project(usa_states_vect,crs(drought_statistics_vect)))

# Select droughts longer than 104 weeks
#multiyear_droughts <- sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,])
multiyear_droughts <- drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]

# Add a 30-m buffer around each point
multiyear_droughts <- buffer(multiyear_droughts, width=30)

# Convert to shapefile
multiyear_droughts <- sf::st_as_sf(multiyear_droughts)

plot(multiyear_droughts["Number_Weeks"])

# Loop through each point and fetch 16-day RAP data
rap_list <- list()

for (i in seq_len(nrow(multiyear_droughts))) {
  roi <- multiyear_droughts[i,]
  message("Fetching RAP data for phenocam: ", roi$Phenocam)
  
  # Fetch RAP 16-day production data (all available years by default)
  rap_data <- trex::fetch_rap(polygons=roi, data_type="production16day", mask=FALSE)
  
  # Add Phenocam name for tracking and store in list
  rap_data$Phenocam <- roi$Phenocam
  rap_list[[i]] <- rap_data
}

# Combine into a single data frame
rap_16day_all_areas <- bind_rows(rap_list)

# Loop through each point and fetch yearly cover data
rap_list <- list()

for (i in seq_len(nrow(multiyear_droughts))) {
  roi <- multiyear_droughts[i,]
  message("Fetching RAP data for phenocam: ", roi$Phenocam)
  
  # Fetch RAP yearly cover data (all available years by default)
  rap_data <- trex::fetch_rap(polygons=roi, data_type="cover", mask=FALSE)
  
  # Add Phenocam name for tracking and store in list
  rap_data$Phenocam <- roi$Phenocam
  rap_list[[i]] <- rap_data
}

# Combine into a single data frame
rap_yearly_cover_all_areas <- bind_rows(rap_list)

# Load SPEI data
spei_1month_df <- data.frame()
spei_3month_df <- data.frame()
spei_6month_df <- data.frame()
spei_12month_df <- data.frame()

# Read in SPEI data
for (filename in list.files("data/spei_1mn_por40_ref30", full.names=TRUE)){
  temp_df <- read.csv(filename)
  spei_1month_df <- rbind(spei_1month_df, temp_df)
}

for (filename in list.files("data/spei_3mn_por40_ref30", full.names=TRUE)){
  temp_df <- read.csv(filename)
  spei_3month_df <- rbind(spei_3month_df, temp_df)
}

for (filename in list.files("data/spei_6mn_por40_ref30", full.names=TRUE)){
  temp_df <- read.csv(filename)
  spei_6month_df <- rbind(spei_6month_df, temp_df)
}

for (filename in list.files("data/spei_12mn_por60_ref30", full.names=TRUE)){
  temp_df <- read.csv(filename)
  spei_12month_df <- rbind(spei_12month_df, temp_df)
}

rm(temp_df, filename)

# Rename variables and remove invalid values
spei_1month_df <- spei_1month_df %>%
  rename(SPEI_1month = values, Phenocam = site, Date = date) %>%
  mutate(Date = as.Date(Date), Year=year(Date), DOY=yday(Date)) %>%
  mutate(SPEI_1month = replace(SPEI_1month, SPEI_1month==-9999, NA))

spei_3month_df <- spei_3month_df %>%
  rename(SPEI_3month = values, Phenocam = site, Date = date) %>%
  mutate(Date = as.Date(Date), Year=year(Date), DOY=yday(Date)) %>%
  mutate(SPEI_3month = replace(SPEI_3month, SPEI_3month==-9999, NA))

spei_6month_df <- spei_6month_df %>%
  rename(SPEI_6month = values, Phenocam = site, Date = date) %>%
  mutate(Date = as.Date(Date), Year=year(Date), DOY=yday(Date)) %>%
  mutate(SPEI_6month = replace(SPEI_6month, SPEI_6month==-9999, NA))

spei_12month_df <- spei_12month_df %>%
  rename(SPEI_12month = values, Phenocam = site, Date = date) %>%
  mutate(Date = as.Date(Date), Year=year(Date), DOY=yday(Date)) %>%
  mutate(SPEI_12month = replace(SPEI_12month, SPEI_12month==-9999, NA))

# Join dataframes into one
spei_df <- left_join(spei_1month_df,spei_3month_df, by=join_by(Phenocam,lat,lon,elev,Date,Year,DOY))
spei_df <- left_join(spei_df,spei_6month_df, by=join_by(Phenocam,lat,lon,elev,Date,Year,DOY))
spei_df <- left_join(spei_df,spei_12month_df, by=join_by(Phenocam,lat,lon,elev,Date,Year,DOY))
spei_df <- select(spei_df, Phenocam, lat, lon, elev, Date, Year, DOY, SPEI_1month, SPEI_3month, SPEI_6month, SPEI_12month)

rm(spei_1month_df, spei_3month_df, spei_6month_df, spei_12month_df)

# Plot production over time

test_data <- filter(rap_16day_all_areas, Phenocam=="ibp", year>=2014, year<2025)
test_data <- pivot_longer(test_data,
                          cols=c("AFG","PFG","HER"),
                          names_to="Plant_Functional_Group",
                          values_to="Production")
test_plot <- ggplot() +
  geom_line(data=test_data, aes(x=date, y=Production, color=Plant_Functional_Group)) +
  geom_vline(xintercept = filter(multiyear_droughts, Phenocam=="ibp")$Start_Date) +
  geom_vline(xintercept = filter(multiyear_droughts, Phenocam=="ibp")$End_Date) +
  #legend(title="Plant Functional Group") +
  ylab("16-day production (lbs/acre)") +
  ggtitle("IBP 16-day production data")
test_plot

# Plot yearly cover over time
test_data_2 <- rap_yearly_cover_all_areas %>%
  filter(Phenocam=="ibp", year>=2014, year<2025) %>%
  pivot_longer(cols=c("AFG", "PFG", "SHR", "TRE", "LTR", "BGR"),
               names_to="Plant_Functional_Group",
               values_to="Percent_Cover")

test_plot_2 <- ggplot() +
  geom_point(data=test_data_2, aes(x=year, y=Percent_Cover, color=Plant_Functional_Group)) +
  geom_vline(xintercept = filter(multiyear_droughts, Phenocam=="ibp")$Start_Year) +
  geom_vline(xintercept = filter(multiyear_droughts, Phenocam=="ibp")$End_Year) +
  ylab("Yearly Percent Cover") +
  ggtitle("IBP yearly cover data")
test_plot_2

# Plot production over time for all 35 sites
for (i in seq_len(nrow(multiyear_droughts))) {
  phenocam <- multiyear_droughts[i,]$Phenocam
  primary_veg <- filter(primary_rois, Phenocam==phenocam)$Primary_Veg_Type
  temp_data <- rap_16day_all_areas %>%
    filter(Phenocam==phenocam, year>=2017, year<2025) %>%
    pivot_longer(cols=c("AFG","PFG","HER"),
                 names_to="Plant_Functional_Group",
                 values_to="Production") %>%
    mutate(Production = replace(Production, Production<0, NA))
  temp_plot <- ggplot() +
    geom_line(data=temp_data, aes(x=date, y=Production, color=Plant_Functional_Group)) +
    geom_vline(xintercept = multiyear_droughts[i,]$Start_Date) +
    geom_vline(xintercept = multiyear_droughts[i,]$End_Date) +
    #legend(title="Plant Functional Group") +
    ylab("Production (lbs/acre)") +
    #ggtitle(paste(phenocam,"16-day production data"))
    ggtitle(paste(phenocam,", ",primary_veg,sep=''))
  
  temp_data_2 <- rap_yearly_cover_all_areas %>%
    filter(Phenocam==phenocam, year>=2017, year<2025) %>%
    pivot_longer(cols=c("AFG", "PFG", "SHR", "TRE"),
                 names_to="Plant_Functional_Group",
                 values_to="Percent_Cover")
  temp_data_2 <- left_join(unique(select(temp_data, date,year,doy)),temp_data_2, by=c("year"), relationship='many-to-many')
  temp_plot_2 <- ggplot() +
    geom_line(data=temp_data_2, aes(x=date, y=Percent_Cover, color=Plant_Functional_Group)) +
    geom_vline(xintercept = multiyear_droughts[i,]$Start_Date) +
    geom_vline(xintercept = multiyear_droughts[i,]$End_Date) +
    ylab("Percent Cover")
    #ylab("Percent Cover") +
    #ggtitle(paste(phenocam,"yearly cover data"))
  
  temp_data_3 <- spei_df %>%
    filter(Phenocam==phenocam, Year>=2017, Year<2025) %>%
    pivot_longer(cols=c("SPEI_3month","SPEI_12month"),
                 names_to="SPEI_timescale",
                 values_to="SPEI_value")
  temp_plot_3 <- ggplot() +
    geom_line(data=temp_data_3, aes(x=Date, y=SPEI_value, color=SPEI_timescale)) +
    ylab("SPEI")
  
  plot <- temp_plot / temp_plot_3
  print(plot)
}

