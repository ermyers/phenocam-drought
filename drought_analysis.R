# Phenocam USDM comparison analysis

# Note 10/29/24: Currently sites with multiple growing seasons (e.g. waahila)
# are causing problems with calculation of typical growing season

# Load packages
library(dplyr)
library(ggplot2)
library(patchwork)
library(terra)
library(sf)
library(lubridate)

# Load data
load("outputs/growing_seasons_phen.RData")
#load("outputs/phen_usdm.RData")
load("outputs/drought_statistics.RData")

###################################
# Filter data according to ROI list
###################################

# Load list of manually filtered PhenoCam ROIs
unique_phenocam_rois_annotated <- read.csv("data/unique_phenocam_rois_annotated.csv")

# Filter list to include only primary ROIs
primary_rois <- filter(unique_phenocam_rois_annotated, action=="keep" & rank==1)

# Filter growing season and drought statistics to include only primary ROIs
growing_seasons_phen_filtered <- data.frame()
growing_seasons_phen_ave_filtered <- data.frame()
growing_seasons_phen_number_filtered <- data.frame()
#drought_statistics_filtered <- data.frame()
#phen_usdm_modified_filtered <- data.frame()

for (i in 1:nrow(primary_rois)){
  phenocam <- primary_rois$Phenocam[i]
  veg_type <- primary_rois$Veg_Type[i]
  roi <- primary_rois$ROI[i]
  
  growing_seasons_phen_filtered <- rbind(growing_seasons_phen_filtered,
                                         filter(growing_seasons_phen,
                                                Phenocam==phenocam,
                                                Veg_Type==veg_type,
                                                ROI==roi))
  
  growing_seasons_phen_ave_filtered <- rbind(growing_seasons_phen_ave_filtered,
                                             filter(growing_seasons_phen_ave,
                                                    Phenocam==phenocam,
                                                    Veg_Type==veg_type,
                                                    ROI==roi))
  
  growing_seasons_phen_number_filtered <- rbind(growing_seasons_phen_number_filtered,
                                                filter(growing_seasons_phen_number,
                                                       Phenocam==phenocam,
                                                       Veg_Type==veg_type,
                                                       ROI==roi))
}

rm(i, phenocam, roi, veg_type)

drought_statistics_filtered <- filter(drought_statistics, Phenocam %in% primary_rois$Phenocam)
phen_usdm_modified_filtered <- filter(phen_usdm_modified, Phenocam %in% primary_rois$Phenocam)

# Remove unfiltered data
rm(drought_statistics, growing_seasons_phen, growing_seasons_phen_ave,
   growing_seasons_phen_number, phen_usdm_modified)

###########################################
# Calculate drought statistics by site-year
###########################################

growing_seasons_phen_usdm_filtered <- cbind(growing_seasons_phen_filtered,Cumulative_DM=NA, Cumulative_DM_with_D0=NA,Drought_Weeks=NA)

# for(i in 1:nrow(growing_seasons_phen_usdm_filtered)){
#   phenocam <- growing_seasons_phen_usdm_filtered$Phenocam[i]
#   year <- growing_seasons_phen_usdm_filtered$Year[i]
#   eos_25 <- growing_seasons_phen_usdm_filtered$EOS_25[i]
#   if(year(eos_25)!=year){
#     # will not calculate cumulative DM if EOS occurs into next calendar year
#   } else{
#     temp_phen_usdm <- filter(phen_usdm_modified_filtered,Phenocam==phenocam,Year==year,Date<=eos_25,DM>=1)
#     cumulative_dm <- sum(temp_phen_usdm$DM,na.rm=TRUE)
#     drought_weeks <- nrow(temp_phen_usdm)
#     growing_seasons_phen_usdm_filtered$Cumulative_DM[i] <- cumulative_dm
#     growing_seasons_phen_usdm_filtered$Drought_Weeks[i] <- drought_weeks
#   }
# }

for(i in 1:nrow(growing_seasons_phen_usdm_filtered)){
  phenocam <- growing_seasons_phen_usdm_filtered$Phenocam[i]
  year <- growing_seasons_phen_usdm_filtered$Year[i]
  temp_phen_usdm <- filter(phen_usdm_modified_filtered,Phenocam==phenocam,Year==year)
  cumulative_dm_with_d0 <- (nrow(filter(temp_phen_usdm,DM==0)) + 2*nrow(filter(temp_phen_usdm,DM==1)) +
    3*nrow(filter(temp_phen_usdm,DM==2)) + 4*nrow(filter(temp_phen_usdm,DM==3)) +
    5*nrow(filter(temp_phen_usdm,DM==4)))
  cumulative_dm_with_d0 <- cumulative_dm_with_d0/nrow(temp_phen_usdm)
  temp_phen_usdm <- filter(phen_usdm_modified_filtered,Phenocam==phenocam,Year==year,DM>=1)
  cumulative_dm <- sum(temp_phen_usdm$DM,na.rm=TRUE)
  drought_weeks <- nrow(temp_phen_usdm)
  growing_seasons_phen_usdm_filtered$Cumulative_DM[i] <- cumulative_dm
  growing_seasons_phen_usdm_filtered$Cumulative_DM_with_D0[i] <- cumulative_dm_with_d0
  growing_seasons_phen_usdm_filtered$Drought_Weeks[i] <- drought_weeks
}

rm(i,phenocam,year,temp_phen_usdm,cumulative_dm,cumulative_dm_with_d0,drought_weeks)

# Plot cumulative DM vs cumulative and peak GCC by veg type
temp_data <- mutate(growing_seasons_phen_usdm_filtered,GSL_25 = EOS_25-SOS_25, SOS_25_DOY = yday(SOS_25))

facet_peak <- ggplot(data=filter(temp_data, Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH"), aes(x=Cumulative_DM,y=Peak_GCC)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("Peak GCC vs. Cumulative DM (not including D0)")
facet_peak

facet_peak <- ggplot(data=filter(temp_data, Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH"), aes(x=Cumulative_DM_with_D0,y=Peak_GCC)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("Peak GCC vs. Weighted Average DM (including D0)")
facet_peak

# Plot residuals for the above by veg type
linear_model_1 <- temp_data %>%
  filter(Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH") %>%
  group_by(Veg_Type) %>%
  do(model = lm(Peak_GCC ~ Cumulative_DM, data = .))

linear_model_2 <- temp_data %>%
  filter(Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH") %>%
  group_by(Veg_Type) %>%
  do(model = lm(Peak_GCC ~ Cumulative_DM_with_D0, data = .))

residuals_df_1 <- data.frame()
residuals_df_2 <- data.frame()

for (i in 1:5){
  veg_type <- linear_model_1$Veg_Type[i]
  cum_dm_1 <- linear_model_1[[2]][[i]][["model"]][["Cumulative_DM"]]
  peak_gcc <- linear_model_1[[2]][[i]][["model"]][["Peak_GCC"]]
  resid_1 <- linear_model_1[[2]][[i]][["residuals"]]
  temp_df <- data.frame("Veg_Type" = veg_type,
                        "Cumulative_DM" = cum_dm_1,
                        "Peak_GCC" = peak_gcc,
                        "Residuals_no_D0" = resid_1)
  residuals_df_1 <- rbind(residuals_df_1,temp_df)
}

for (i in 1:5){
  veg_type <- linear_model_2$Veg_Type[i]
  peak_gcc <- linear_model_2[[2]][[i]][["model"]][["Peak_GCC"]]
  cum_dm_2 <- linear_model_2[[2]][[i]][["model"]][["Cumulative_DM_with_D0"]]
  resid_2 <- linear_model_2[[2]][[i]][["residuals"]]
  temp_df <- data.frame("Veg_Type" = veg_type,
                        "Cumulative_DM_with_D0" = cum_dm_2,
                        "Peak_GCC" = peak_gcc,
                        "Residuals_with_D0" = resid_2)
  residuals_df_2 <- rbind(residuals_df_2,temp_df)
}

rm(cum_dm_1,cum_dm_2,i,peak_gcc,resid_1,resid_2,veg_type)

facet_resid <- ggplot(data=residuals_df_1, aes(x=Cumulative_DM,y=Residuals_no_D0)) +
  geom_point() +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("Residuals for Peak GCC vs. Cumulative DM (not including D0)")
facet_resid

facet_resid <- ggplot(data=residuals_df_2, aes(x=Cumulative_DM_with_D0,y=Residuals_with_D0)) +
  geom_point() +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("Residuals for Peak GCC vs. Weighted Average DM (including D0)")
facet_resid

# Linear model summaries (cumulative DM, no D0)
# AG
summary(lm(Peak_GCC ~ Cumulative_DM, filter(temp_data, Veg_Type=="AG")))
# DB
summary(lm(Peak_GCC ~ Cumulative_DM, filter(temp_data, Veg_Type=="DB")))
# EN
summary(lm(Peak_GCC ~ Cumulative_DM, filter(temp_data, Veg_Type=="EN")))
# GR
summary(lm(Peak_GCC ~ Cumulative_DM, filter(temp_data, Veg_Type=="GR")))
# SH
summary(lm(Peak_GCC ~ Cumulative_DM, filter(temp_data, Veg_Type=="SH")))

# Linear model summaries (cumulative DM, including D0)
# AG
summary(lm(Peak_GCC ~ Cumulative_DM_with_D0, filter(temp_data, Veg_Type=="AG")))
# DB
summary(lm(Peak_GCC ~ Cumulative_DM_with_D0, filter(temp_data, Veg_Type=="DB")))
# EN
summary(lm(Peak_GCC ~ Cumulative_DM_with_D0, filter(temp_data, Veg_Type=="EN")))
# GR
summary(lm(Peak_GCC ~ Cumulative_DM_with_D0, filter(temp_data, Veg_Type=="GR")))
# SH
summary(lm(Peak_GCC ~ Cumulative_DM_with_D0, filter(temp_data, Veg_Type=="SH")))

facet_gsl <- ggplot(data=filter(temp_data, Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH"), aes(x=Cumulative_DM,y=GSL_25)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("GSL vs. Cumulative DM (not including D0)")
facet_gsl

facet_gsl <- ggplot(data=filter(temp_data, Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH"), aes(x=Cumulative_DM_with_D0,y=GSL_25)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Veg_Type,ncol=2) +
  ggtitle("GSL vs. Weighted Average DM (including D0)")
facet_gsl

for(veg_type in unique(temp_data$Veg_Type)){
  p1 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=Cumulative_GCC_yearly)) +
    geom_point() +
    ggtitle(paste("Veg type =",veg_type))
  
  p2 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=Peak_GCC)) +
    geom_point()
  
  p3 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=GSL_25)) +
    geom_point()
  
  p4 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=SOS_25_DOY)) +
    geom_point()
  
  print((p1 + p2)/(p3 + p4))
}

# Plot cumulative DM vs cumulative and peak GCC by ecoregion
temp_data <- left_join(temp_data,distinct(select(phen_usdm_modified_filtered,Phenocam,NA_L1NAME)))

for(ecoregion in unique(temp_data$NA_L1NAME)){
  p1 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=Cumulative_GCC_yearly)) +
    geom_point() +
    ggtitle(paste("Ecoregion =",ecoregion))
  
  p2 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=Peak_GCC)) +
    geom_point()
  
  p3 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=GSL_25)) +
    geom_point()
  
  p4 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=SOS_25_DOY)) +
    geom_point()
  
  print((p1 + p2)/(p3 + p4))
}

# Plot cumulative DM and peak GCC by veg type, faceted

# ################################
# # Phenocam Drought Indicator #1:
# # No Growing Season
# ################################
# 
# # No growing season site-years
# no_gs <- filter(growing_seasons_phen_number, Number_Growing_Seasons==0)
# 
# # Add average growing season information
# no_gs <- left_join(no_gs, growing_seasons_phen_ave, join_by(Phenocam,Veg_Type,ROI))
# 
# # Remove entries with fewer than 5 years of data
# no_gs <- filter(no_gs, Number_Years>=5)
# 
# # Add columns for USDM statistics
# no_gs <- cbind(no_gs, DM_total_yearly=0, Drought_Weeks_yearly=0,
#                DM_total_gs=0, Drought_Weeks_gs=0, Weeks_gs=0)
# 
# # Calculate USDM statistics for calendar year and for expected growing season
# for(i in 1:nrow(no_gs)){
#   phenocam <- no_gs$Phenocam[i]
#   year <- no_gs$Year[i]
#   temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
#   temp_gs_df <- filter(temp_df, DOY >= no_gs$SOS_25_mean[i] & DOY <= no_gs$EOS_25_mean[i])
#   
#   no_gs$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
#   no_gs$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
#   no_gs$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
#   no_gs$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
#   no_gs$Weeks_gs[i] <- nrow(temp_gs_df)
# }
# 
# rm(i,phenocam,year,temp_df,temp_gs_df)
# 
# # Calculate USDM statistics for same sites as above, other years
# gs_df <- right_join(growing_seasons_phen, distinct(select(no_gs,Phenocam,Veg_Type,ROI)), join_by(Phenocam,Veg_Type,ROI))
# 
# # Add columns for USDM statistics
# gs_df <- cbind(gs_df, DM_total_yearly=0, Drought_Weeks_yearly=0,
#                DM_total_gs=0, Drought_Weeks_gs=0, Weeks_gs=0)
# 
# # Calculate USDM statistics for calendar year and for calculated growing season
# for(i in 1:nrow(gs_df)){
#   phenocam <- gs_df$Phenocam[i]
#   year <- gs_df$Year[i]
#   temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
#   temp_gs_df <- filter(temp_df, Date >= gs_df$SOS_25[i] & Date <= gs_df$EOS_25[i])
#   
#   gs_df$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
#   gs_df$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
#   gs_df$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
#   gs_df$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
#   gs_df$Weeks_gs[i] <- nrow(temp_gs_df)
# }
# 
# rm(i,phenocam,year,temp_df,temp_gs_df)

################################
# Phenocam Drought Indicator #2:
# Low Growing Season GCC
################################

# Add columns to growing_seasons_phen_number
growing_seasons_phen_number <- cbind(growing_seasons_phen_number,
                                     Peak_GCC_Percent_of_Ave=0, Peak_GCC_STD_from_Ave = 0)

# Calculate growing season peak GCC as a percentage of the average
for(i in 1:nrow(growing_seasons_phen_number)){
  phenocam <- growing_seasons_phen_number$Phenocam[i]
  veg_type <- growing_seasons_phen_number$Veg_Type[i]
  roi <- growing_seasons_phen_number$ROI[i]
  peak_gcc_mean <- filter(growing_seasons_phen_ave, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi)$Peak_GCC_mean
  peak_gcc_std <- filter(growing_seasons_phen_ave, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi)$Peak_GCC_std
  
  growing_seasons_phen_number$Peak_GCC_Percent_of_Ave[i] <- growing_seasons_phen_number$Peak_GCC[i]/peak_gcc_mean
  growing_seasons_phen_number$Peak_GCC_STD_from_Ave[i] <- (growing_seasons_phen_number$Peak_GCC[i]-peak_gcc_mean)/peak_gcc_std
}

rm(i, phenocam, veg_type, roi, peak_gcc_mean, peak_gcc_std)

# Low peak GCC site-years
low_gcc <- filter(growing_seasons_phen_number, Peak_GCC_STD_from_Ave <= -1)

# Add columns for USDM statistics
low_gcc <- cbind(low_gcc, DM_total_yearly=0, Drought_Weeks_yearly=0,
                 DM_total_gs=NA, Drought_Weeks_gs=NA, Weeks_gs=NA)

# Calculate USDM statistics for calendar year
for(i in 1:nrow(low_gcc)){
  phenocam <- low_gcc$Phenocam[i]
  veg_type <- low_gcc$Veg_Type[i]
  roi <- low_gcc$ROI[i]
  year <- low_gcc$Year[i]
  temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
  low_gcc$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
  low_gcc$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
  if(low_gcc$Number_Growing_Seasons[i]==1){
    sos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$SOS_25
    eos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$EOS_25
    temp_gs_df <- filter(temp_df, Date >= sos_25 & Date <= eos_25)
    
    low_gcc$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
    low_gcc$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
    low_gcc$Weeks_gs[i] <- nrow(temp_gs_df)
  }
}

rm(i, phenocam, roi, veg_type, year, temp_df, temp_gs_df, sos_25, eos_25)

# Calculate USDM statistics for same sites as above, other years
normal_gcc <- filter(growing_seasons_phen_number, Peak_GCC_STD_from_Ave > -1)
normal_gcc <- right_join(normal_gcc, distinct(select(low_gcc,Phenocam,Veg_Type,ROI)), join_by(Phenocam,Veg_Type,ROI))

# Add columns for USDM statistics
normal_gcc <- cbind(normal_gcc, DM_total_yearly=0, Drought_Weeks_yearly=0,
                    DM_total_gs=NA, Drought_Weeks_gs=NA, Weeks_gs=NA)

# Calculate USDM statistics for calendar year
for(i in 1:nrow(normal_gcc)){
  phenocam <- normal_gcc$Phenocam[i]
  veg_type <- normal_gcc$Veg_Type[i]
  roi <- normal_gcc$ROI[i]
  year <- normal_gcc$Year[i]
  temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
  normal_gcc$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
  normal_gcc$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
  
  if(normal_gcc$Number_Growing_Seasons[i]==1){
    sos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$SOS_25
    eos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$EOS_25
    temp_gs_df <- filter(temp_df, Date >= sos_25 & Date <= eos_25)
    
    normal_gcc$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
    normal_gcc$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
    normal_gcc$Weeks_gs[i] <- nrow(temp_gs_df)
  }
}

rm(i, phenocam, veg_type, roi, year, temp_df, temp_gs_df, sos_25, eos_25)

#########################
# Plot PhenoCam Locations
#########################

# Load shapefiles
type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
                                     usa_states_vect$NAME != "Hawaii" &
                                     usa_states_vect$NAME != "Puerto Rico",]

# Convert phenocam shapefile into dataframe with coordinates
type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to primary ROIs
primary_rois_df <- left_join(primary_rois,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert primary ROIs into vector
primary_rois_vect <- vect(primary_rois_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Convert primary ROIs into sf
primary_rois_sf <- sf::st_as_sf(primary_rois_vect)

# Plot primary ROIs (terra)
plot(usa_states_vect)
plot(primary_rois_vect,2,add=TRUE)

# Plot primary ROIs (ggplot)
ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = primary_rois_sf, aes(color=Veg_Type)) +
  #scale_color_viridis_c(option = 'viridis') +
  #scale_color_stepsn(colors=c("purple","blue","green","yellow","orange","red","brown"),
  #                   breaks=c(0,10,20,30,40,50,100,200),
  #                   guide = "colorsteps") +
  labs(color = 'Veg Type') +
  ggtitle("PhenoCam locations")

########################
# Plot Drought Locations
########################

# # Load shapefiles
# type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
# usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
# usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
#                                      usa_states_vect$NAME != "Hawaii" &
#                                      usa_states_vect$NAME != "Puerto Rico",]
# 
# # Convert phenocam shapefile into dataframe with coordinates
# type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to drought statistics
drought_statistics_df <- left_join(drought_statistics_filtered,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert drought statistics into vector
drought_statistics_vect <- vect(drought_statistics_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Convert drought statistics into sf
drought_statistics_sf <- sf::st_as_sf(drought_statistics_vect)

# Plot drought statistics (terra)
plot(usa_states_vect)
plot(drought_statistics_vect,6,add=TRUE) # Number of weeks of drought
plot(usa_states_vect)
plot(drought_statistics_vect,7,add=TRUE) # Sum of drought monitor index

# Plot drought statistics (ggplot)
ggplot(data = drought_statistics_sf, aes(color=Number_Weeks)) +
  #scale_color_distiller(palette = 'YlOrRd', direction = 1) +
  scale_color_viridis_c(option = 'viridis') +
  geom_sf() +
  labs(color = 'Drought length (weeks)') +
  ggtitle("Number of Weeks of Drought")

ggplot(data = drought_statistics_sf, aes(color=Sum_USDM)) +
  #scale_color_distiller(palette = 'YlOrRd', direction = 1) +
  scale_color_viridis_c(option = 'viridis') +
  geom_sf() +
  labs(color = 'Sum of drought index') +
  ggtitle("Cumulative drought severity")

# Plot drought statistics (ggplot, contiguous USA)
drought_statistics_contig_vect <- terra::intersect(drought_statistics_vect,project(usa_states_vect,crs(drought_statistics_vect)))
drought_statistics_contig_sf <- sf::st_as_sf(drought_statistics_contig_vect)
drought_statistics_contig_sf <- mutate(drought_statistics_contig_sf, Number_Weeks_Binned = cut(Number_Weeks, breaks=c(0,10,20,40,60,100,200)))
drought_statistics_contig_sf <- mutate(drought_statistics_contig_sf, Sum_USDM_Binned = cut(Sum_USDM, breaks=c(0,25,50,100,200,400,600)))

ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = drought_statistics_contig_sf, aes(color=Number_Weeks_Binned)) +
  #scale_color_viridis_c(option = 'viridis') +
  #scale_color_stepsn(colors=c("purple","blue","green","yellow","orange","red","brown"),
  #                   breaks=c(0,10,20,30,40,50,100,200),
  #                   guide = "colorsteps") +
  labs(color = 'Drought length (weeks)') +
  ggtitle("Number of Weeks of Continuous Drought")

ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = drought_statistics_contig_sf, aes(color=Sum_USDM_Binned)) +
  #scale_color_viridis_c(option = 'viridis') +
  labs(color = 'Sum of drought index') +
  ggtitle("Cumulative Severity (DM) of Continuous Drought")

#######################
# Plot no-gs site-years
#######################

# Load shapefiles
# type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
# usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
# usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska"]

# Convert phenocam shapefile into dataframe with coordinates
# type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to no_gs
no_gs_df <- left_join(no_gs,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert drought statistics into vector
no_gs_vect <- vect(no_gs_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Plot no gs site-years
plot(usa_states_vect)
plot(no_gs_vect,2,add=TRUE) # Veg Type
plot(no_gs_vect,2) # Veg Type

##########################################
# Plot SOS vs. EOS with standard deviation
##########################################

# plot of SOS vs EOS, colored by veg type
ggplot(data = growing_seasons_phen_ave_filtered, aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  ggtitle("All Primary ROIs")

# plot without error bars, filtered so that max stdev is no more than 20 days
ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=20 & EOS_25_std<=20), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  ggtitle("Phenocam ROIs with SOS and EOS stdev <20 days")

# plot of SOS vs EOS, with error bars
ggplot(data = growing_seasons_phen_ave_filtered, aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
  geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
  ggtitle("All Phenocam ROIs")

# plot with error bars, but filtered so that max stdev is no more than 30 days
ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=30 & EOS_25_std<=30), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
  geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
  ggtitle("SOS and EOS stdev <30 days")

# plot with error bars, but filtered so that max stdev is no more than 20 days
ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=20 & EOS_25_std<=20), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
  geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
  ggtitle("SOS and EOS stdev <20 days")

# plot with error bars, but filtered so that max stdev is no more than 10 days
ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=10 & EOS_25_std<=10), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
  geom_point() +
  geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
  geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
  ggtitle("SOS and EOS stdev <10 days")

####################
# Preliminary plots
####################

drought_indicator_1 <- right_join(select(no_gs, Phenocam, Veg_Type, ROI, Year, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
                                  select(gs_df, Phenocam, Veg_Type, ROI, Year, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
                                  join_by(Phenocam,Veg_Type,ROI))

p1 <- ggplot(data = drought_indicator_1, aes(x = DM_total_yearly.x, y = DM_total_yearly.y)) +
  geom_point() +
  xlab("Total DM for no growing season year") +
  ylab("Total DM for growing season years, same site")
p2 <- ggplot(data = drought_indicator_1, aes(x = Drought_Weeks_yearly.x, y = Drought_Weeks_yearly.y)) +
  geom_point() +
  xlab("# Weeks in Drought for no growing season year") +
  ylab("# Weeks in Drought for growing season years, same site")
p1
p2

drought_indicator_2 <- right_join(select(low_gcc, Phenocam, Veg_Type, ROI, Year, Peak_GCC_Percent_of_Ave, Peak_GCC_STD_from_Ave, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
                                  select(normal_gcc, Phenocam, Veg_Type, ROI, Year, Peak_GCC_Percent_of_Ave, Peak_GCC_STD_from_Ave, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
                                  join_by(Phenocam,Veg_Type,ROI), relationship = "many-to-many")

p3 <- ggplot(data = drought_indicator_2, aes(x = DM_total_yearly.x, y = DM_total_yearly.y)) +
  geom_point() +
  xlab("Total DM for low GCC year") +
  ylab("Total DM for normal GCC years, same site")
p3

p4 <- ggplot() +
  geom_point(data = low_gcc, aes(x = Peak_GCC_STD_from_Ave, y = DM_total_yearly), color = "red") +
  geom_point(data = normal_gcc, aes(x = Peak_GCC_STD_from_Ave, y = DM_total_yearly), color = "black") +
  xlim(-10,10) +
  xlab("Standard deviations of peak GCC away from mean") +
  ylab("Total DM for site-year")
p4

p5 <- ggplot() +
  geom_point(data = low_gcc, aes(x = Peak_GCC_STD_from_Ave, y = Drought_Weeks_yearly), color = "red") +
  geom_point(data = normal_gcc, aes(x = Peak_GCC_STD_from_Ave, y = Drought_Weeks_yearly), color = "black") +
  xlim(-10,10) +
  xlab("Standard deviations of peak GCC away from mean") +
  ylab("# Weeks in drought for site-year")
p5

###
#
#