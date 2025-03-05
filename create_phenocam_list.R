# Come up with an initial list of study phenocams, based on:
# 1) Type 1
# 2) Multiple years of data (3 or 5 to start)
# 3) Not NEON understory or spruce experimental manipulation plot
# 4) In the USA (for USDM analysis)

# The list will be further refined manually later on in the process.

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

# Remove NEON DP1.00042 understory
type1_sites_no_understory <- dplyr::filter(type1_sites,!grepl('DP1.00042', site))

# Remove spruce experimental plots
type1_sites_no_spruce <- dplyr::filter(type1_sites_no_understory,!grepl('spruceT', site))

# Type I sites with 3, 4, and 5+ years of data
type1_sites_3years <- dplyr::filter(type1_sites_no_spruce,as.Date(date_last)-as.Date(date_first)>=1095)
type1_sites_4years <- dplyr::filter(type1_sites_no_spruce,as.Date(date_last)-as.Date(date_first)>=1460)
type1_sites_5years <- dplyr::filter(type1_sites_no_spruce,as.Date(date_last)-as.Date(date_first)>=1826)

# To do USA-only, we will need to convert our dataframes into vectors
type1_sites_vect <- vect(type1_sites, geom=c("lon","lat"), crs="epsg:4326")
type1_sites_3years_vect <- vect(type1_sites_3years, geom=c("lon","lat"), crs="epsg:4326")
type1_sites_4years_vect <- vect(type1_sites_4years, geom=c("lon","lat"), crs="epsg:4326")
type1_sites_5years_vect <- vect(type1_sites_5years, geom=c("lon","lat"), crs="epsg:4326")

# Load in USA county data...
usa_counties_2023 <- vect("data/shapefiles/tl_2023_us_county.shp")
# ...and reproject it into the same CRS as the phenocams
usa_counties_2023 <- project(usa_counties_2023,type1_sites_vect)

# Subset phenocams using USA counties
type1_sites_usa_vect <- type1_sites_vect[usa_counties_2023]
type1_sites_3years_usa_vect <- type1_sites_3years_vect[usa_counties_2023]
type1_sites_4years_usa_vect <- type1_sites_4years_vect[usa_counties_2023]
type1_sites_5years_usa_vect <- type1_sites_5years_vect[usa_counties_2023]

# Convert vectors back to dataframes
type1_sites_usa <- as.data.frame(type1_sites_usa_vect)
type1_sites_3years_usa <- as.data.frame(type1_sites_3years_usa_vect)
type1_sites_4years_usa <- as.data.frame(type1_sites_4years_usa_vect)
type1_sites_5years_usa <- as.data.frame(type1_sites_5years_usa_vect)

# Get list of site names from dataframes
type1_sitenames_usa <- type1_sites_usa$site
type1_sitenames_3years_usa <- type1_sites_3years_usa$site
type1_sitenames_4years_usa <- type1_sites_4years_usa$site
type1_sitenames_5years_usa <- type1_sites_5years_usa$site

# Save lists
save(list=c("type1_sitenames_usa","type1_sitenames_3years_usa","type1_sitenames_4years_usa","type1_sitenames_5years_usa"), file="outputs/phenocam_lists.RData")

# Write and save shapefiles of phenocam data
# phenocam
writeVector(type1_sites_usa_vect,"outputs/shapefiles/type1_sites_usa.shp", overwrite=TRUE)
writeVector(type1_sites_3years_usa_vect,"outputs/shapefiles/type1_sites_3years_usa.shp", overwrite=TRUE)
writeVector(type1_sites_4years_usa_vect,"outputs/shapefiles/type1_sites_4years_usa.shp", overwrite=TRUE)
writeVector(type1_sites_5years_usa_vect,"outputs/shapefiles/type1_sites_5years_usa.shp", overwrite=TRUE)

# # county
# usa_counties_type1 <- usa_counties_2023[type1_sites_usa_vect]
# writeVector(usa_counties_type1, "outputs/shapefiles/usa_counties_type1.shp", overwrite=TRUE)
# usa_counties_3years_type1 <- usa_counties_2023[type1_sites_3years_usa_vect]
# writeVector(usa_counties_2023[type1_sites_3years_usa_vect], "outputs/shapefiles/usa_counties_3years_type1.shp")