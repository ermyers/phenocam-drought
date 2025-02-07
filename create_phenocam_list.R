# Come up with list (or multiple lists) of valid phenocams for study, based on:
# 1) Type 1
# 2) Multiple years of data (3 or 5 to start)
# 3) In the USA (for USDM analysis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(phenocamapi)
# installed via following code lines:
# if(!require(devtools)) install.packages('devtools')
# devtools::install_github('bnasr/phenocamapi')
library(dplyr)
library(terra)

# Download list of phenocam sites
site_metadata <- get_phenos()

# Filter list of phenocam sites
# All Type I sites
type1_sites <- dplyr::filter(site_metadata,site_type=="I")

# Type I sites with 3 and 5+ years of data
type1_sites_3years <- dplyr::filter(type1_sites,as.Date(date_last)-as.Date(date_first)>=1095)
type1_sites_5years <- dplyr::filter(type1_sites,as.Date(date_last)-as.Date(date_first)>=1826)

# To do USA-only, we will need to convert our dataframes into vectors
type1_sites_vect <- vect(type1_sites, geom=c("lon","lat"), crs="epsg:4326")
type1_sites_3years_vect <- vect(type1_sites_3years, geom=c("lon","lat"), crs="epsg:4326")
type1_sites_5years_vect <- vect(type1_sites_5years, geom=c("lon","lat"), crs="epsg:4326")

# Load in USA county data...
usa_counties_2023 <- vect("data/shapefiles/tl_2023_us_county.shp")
# ...and reproject it into the same CRS as the phenocams
usa_counties_2023 <- project(usa_counties_2023,type1_sites_vect)

# Subset phenocams using USA counties
type1_sites_usa_vect <- type1_sites_vect[usa_counties_2023]
type1_sites_3years_usa_vect <- type1_sites_3years_vect[usa_counties_2023]
type1_sites_5years_usa_vect <- type1_sites_5years_vect[usa_counties_2023]

# Convert vectors back to dataframes
type1_sites_usa <- as.data.frame(type1_sites_usa_vect)
type1_sites_3years_usa <- as.data.frame(type1_sites_3years_usa_vect)
type1_sites_5years_usa <- as.data.frame(type1_sites_5years_usa_vect)

# Get list of site names from dataframes
type1_sitenames_usa <- type1_sites_usa$site
type1_sitenames_3years_usa <- type1_sites_3years_usa$site
type1_sitenames_5years_usa <- type1_sites_5years_usa$site

# Split list of sitenames to download in subsets
# (right now looking at 3-year data, but could do with any list)
list1 <- type1_sitenames_3years_usa[1:110]
list2 <- type1_sitenames_3years_usa[111:220]
list3 <- type1_sitenames_3years_usa[221:330]
list4 <- type1_sitenames_3years_usa[331:441]

# Save lists
save(list=c("type1_sitenames_usa","type1_sitenames_3years_usa","type1_sitenames_5years_usa","list1","list2","list3","list4"), file="outputs/phenocam_lists_split.RData")
save(list=c("type1_sitenames_usa","type1_sitenames_3years_usa","type1_sitenames_5years_usa"), file="outputs/phenocam_lists.RData")

# Write and save shapefiles of phenocam and county data
# phenocam
writeVector(type1_sites_usa_vect,"outputs/shapefiles/type1_sites_usa.shp")
writeVector(type1_sites_3years_usa_vect,"outputs/shapefiles/type1_sites_3years_usa.shp")
writeVector(type1_sites_5years_usa_vect,"outputs/shapefiles/type1_sites_5years_usa.shp")
# county
usa_counties_type1 <- usa_counties_2023[type1_sites_usa_vect]
writeVector(usa_counties_type1, "outputs/shapefiles/usa_counties_type1.shp")
usa_counties_3years_type1 <- usa_counties_2023[type1_sites_3years_usa_vect]
writeVector(usa_counties_2023[type1_sites_3years_usa_vect], "outputs/shapefiles/usa_counties_3years_type1.shp")