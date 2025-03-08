# Generate list of unique PhenoCam ROIs and save as a CSV for manual annotation
# Use 2016-2024 data

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(dplyr)

# Load growing season data
load('outputs/growing_seasons_phen_2016_to_2024.RData')

# Load phenocam site information (for lat/lon)
load("outputs/type1_sites.RData")

# Make a list of all phenocam ROIs with at least 4 growing seasons between 2016 and 2024
unique_phenocam_rois <- growing_seasons_phen_ave %>% filter(Number_Years >= 4) %>% distinct(Phenocam,Veg_Type,ROI)

# For each phenocam ROI, append lat and lon
unique_phenocam_rois <- cbind(unique_phenocam_rois, lat=NA, lon=NA)
for (i in 1:nrow(unique_phenocam_rois)){
  lat <- filter(type1_sites, site==unique_phenocam_rois$Phenocam[i])$lat
  lon <- filter(type1_sites, site==unique_phenocam_rois$Phenocam[i])$lon
  unique_phenocam_rois$lat[i] <- lat
  unique_phenocam_rois$lon[i] <- lon
}
rm(lat,lon,i)

# For each phenocam ROI, calculate start and end year
unique_phenocam_rois <- cbind(unique_phenocam_rois, start_year=NA, end_year=NA)
for (i in 1:nrow(unique_phenocam_rois)){
  temp_df <- filter(growing_seasons_phen_number,
                    Phenocam==unique_phenocam_rois$Phenocam[i],
                    Veg_Type==unique_phenocam_rois$Veg_Type[i],
                    ROI==unique_phenocam_rois$ROI[i])
  start_year <- min(temp_df$Year, na.rm = TRUE)
  end_year <- max(temp_df$Year, na.rm = TRUE)
  unique_phenocam_rois$start_year[i] <- start_year
  unique_phenocam_rois$end_year[i] <- end_year
}

rm(temp_df, i, start_year, end_year)

# Save outputs to CSV
write.csv(unique_phenocam_rois,"outputs/unique_phenocam_rois.csv", row.names = FALSE)