# Extract Level 1 ecoregions at phenocam locations

library(terra)

# Load in the phenocam locations
type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
crs(type1_sites_usa_vect)
type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect)

# Load in the ecoregions
level_1_ecoregions_vect <- vect("data/ecoregions/NA_CEC_Eco_Level1.shp")
crs(level_1_ecoregions_vect)
level_1_ecoregions_vect <- project(level_1_ecoregions_vect,type1_sites_usa_vect)
crs(level_1_ecoregions_vect)
level_1_ecoregions_df <- as.data.frame(level_1_ecoregions_vect)

# Extract ecoregion values at phenocam locations
phen_with_ecoregion <- terra::intersect(type1_sites_usa_vect,level_1_ecoregions_vect)
phen_with_ecoregion_df <- as.data.frame(phen_with_ecoregion)

# Plot
plot(phen_with_ecoregion,39)

# Save
save(phen_with_ecoregion_df,file="outputs/phen_with_ecoregion.RData")