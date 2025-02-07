# Created 11/20/2024

library(dplyr)

# Load GCC data
load('outputs/phen_gcc.RData')

# Load phenocam site information (for lat/lon)
type1_sites_usa <- read.csv("outputs/type1_sites_usa.csv")

# Make a list of all phenocam ROIs
unique_phenocam_rois <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)

# For each phenocam ROI, append lat and lon
unique_phenocam_rois <- cbind(unique_phenocam_rois, lat=NA, lon=NA)
for (i in 1:nrow(unique_phenocam_rois)){
  lat <- filter(type1_sites_usa, site==unique_phenocam_rois$Phenocam[i])$y
  lon <- filter(type1_sites_usa, site==unique_phenocam_rois$Phenocam[i])$x
  unique_phenocam_rois$lat[i] <- lat
  unique_phenocam_rois$lon[i] <- lon
}
rm(lat,lon,i)

# For each phenocam ROI, calculate start and end year
unique_phenocam_rois <- cbind(unique_phenocam_rois, start_year=NA, end_year=NA)
for (i in 1:nrow(unique_phenocam_rois)){
  filtered_gcc <- filter(phen_gcc,
                         Phenocam==unique_phenocam_rois$Phenocam[i],
                         Veg_Type==unique_phenocam_rois$Veg_Type[i],
                         ROI==unique_phenocam_rois$ROI[i],
                         is.na(gcc_90)==FALSE)
  start_year <- min(filtered_gcc$Year, na.rm = TRUE)
  end_year <- max(filtered_gcc$Year, na.rm = TRUE)
  unique_phenocam_rois$start_year[i] <- start_year
  unique_phenocam_rois$end_year[i] <- end_year
}

rm(filtered_gcc, i, start_year, end_year)

# Make a list of all phenocams with total number of ROIs
phenocam_roi_count <- count(unique_phenocam_rois,Phenocam)

# Make a list of all phenocams with total number of veg types
phenocam_veg_count <- count(distinct(phen_gcc,Phenocam,Veg_Type),Phenocam)

# Save outputs to CSV
write.csv(unique_phenocam_rois,"outputs/unique_phenocam_rois.csv", row.names = FALSE)

# # Calculate phenocams with multiple ROIs
# unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)
# phenocam_roi_count <- count(unique_phenocams,Phenocam)
# multiple_roi_phenocams <- filter(unique_phenocams,Phenocam %in% filter(phenocam_roi_count,n>1)$Phenocam)
# multiple_roi_phenocams_condensed <- multiple_roi_phenocams %>% select(Phenocam,Veg_Type) %>% distinct(Phenocam,Veg_Type)
# multiple_roi_phenocams_list <- unique(multiple_roi_phenocams$Phenocam)
# write.csv(multiple_roi_phenocams,"outputs/multiple_roi_phenocams.csv")
# write.csv(multiple_roi_phenocams_condensed,"outputs/multiple_roi_phenocams_condensed.csv")
# 
# # Calculate phenocams with multiple veg types
# unique_phenocams_veg <- phen_gcc %>% distinct(Phenocam,Veg_Type)
# phenocam_veg_count <- count(unique_phenocams_veg,Phenocam)
# multiple_veg_phenocams <- filter(unique_phenocams_veg,Phenocam %in% filter(phenocam_veg_count,n>1)$Phenocam)
# multiple_veg_phenocams_list <- unique(multiple_veg_phenocams$Phenocam)
# write.csv(multiple_veg_phenocams,"outputs/multiple_veg_phenocams.csv")
