# Phenocam USDM comparison analysis

# Note 04/30/25: Redid primary ROIs to remove Alaska sites.

# Note 04/23/25: Redid baseline subtraction to remove problem with "vaira".

# Note 04/21/25: Trying code with baseline-subtracted GCC instead of raw values.
# Noticed that at least one site ("vaira") had a problem with baseline
# subtraction due to location of peak.

# Note 04/09/25: Designated "primary" growing seasons (highest peak) to deal
# with issues (skewed SOS and EOS, skewed typical growing season) relating to
# multiple growing season site-years. Statistical calculations are currently
# performed on primary growing seasons only.

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
load("outputs/growing_seasons_phen_2016_to_2024_baseline_subtracted.RData")
#load("outputs/drought_statistics.RData")
load("outputs/usdm_statistics.RData")
load("outputs/vpd_statistics.RData")
load("outputs/spei_statistics.RData")
load("outputs/phen_with_ecoregion.RData")

###################################
# Filter data according to ROI list
###################################

# Load list of manually filtered PhenoCam ROIs
unique_phenocam_rois_annotated <- read.csv("data/unique_phenocam_rois_annotated_04-30-2025.csv")

# Filter list to include only primary ROIs
primary_rois <- filter(unique_phenocam_rois_annotated, action=="keep" & rank==1)

# Filter growing season and drought statistics to include only primary ROIs
growing_seasons_phen_filtered <- data.frame()
growing_seasons_phen_ave_filtered <- data.frame()
site_year_statistics_filtered <- data.frame()

# For growing season statistics, need to filter by Phenocam, Veg Type, and ROI
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
  
  site_year_statistics_filtered <- rbind(site_year_statistics_filtered,
                                         filter(site_year_statistics,
                                                Phenocam==phenocam,
                                                Veg_Type==veg_type,
                                                ROI==roi))
}

rm(i, phenocam, roi, veg_type)

# For drought metrics, only need to filter by Phenocam
#drought_statistics_filtered <- filter(drought_statistics, Phenocam %in% primary_rois$Phenocam)
vpd_statistics_filtered <- filter(vpd_statistics, Phenocam %in% primary_rois$Phenocam)
spei_statistics_filtered <- filter(spei_statistics, Phenocam %in% primary_rois$Phenocam)
#phen_usdm_modified_filtered <- filter(phen_usdm_modified, Phenocam %in% primary_rois$Phenocam)
usdm_statistics_filtered <- filter(usdm_statistics, Phenocam %in% primary_rois$Phenocam)

# Remove unfiltered data
# rm(drought_statistics, growing_seasons_phen, growing_seasons_phen_ave,
#    phen_usdm_modified, site_year_statistics, spei_statistics, vpd_statistics)
rm(growing_seasons_phen, growing_seasons_phen_ave, site_year_statistics,
   spei_statistics, usdm_statistics, vpd_statistics)

# Add Primary_Veg_Type field to filtered dataframes
# Primary_Veg_Type is the same as Veg_Type for most ROIs, but assigns a veg type to XX ROIs and a few other ROIs
growing_seasons_phen_ave_filtered <- left_join(growing_seasons_phen_ave_filtered,select(primary_rois,Phenocam,Veg_Type,Primary_Veg_Type),by=join_by(Phenocam,Veg_Type))
growing_seasons_phen_filtered <- left_join(growing_seasons_phen_filtered,select(primary_rois,Phenocam,Veg_Type,Primary_Veg_Type),by=join_by(Phenocam,Veg_Type))
site_year_statistics_filtered <- left_join(site_year_statistics_filtered,select(primary_rois,Phenocam,Veg_Type,Primary_Veg_Type),by=join_by(Phenocam,Veg_Type))

# Add difference from mean to growing_seasons_phen_filtered
# Add mean values
site_year_statistics_filtered <- left_join(site_year_statistics_filtered, select(growing_seasons_phen_ave_filtered, Phenocam,Veg_Type,ROI,Peak_GCC_mean,Peak_GCC_baseline_subtracted_mean,Amplitude_GCC_yearly_mean,Cumulative_GCC_yearly_mean,SOS_25_mean),by=join_by(Phenocam,Veg_Type,ROI))
# Calculate difference (percent of mean for GCC metrics, number of days difference for SOS)
site_year_statistics_filtered <- site_year_statistics_filtered %>% mutate(Peak_GCC_percent_of_mean = (Peak_GCC-Peak_GCC_mean)/Peak_GCC_mean,
                                                                          Amplitude_GCC_yearly_percent_of_mean = (Amplitude_GCC_yearly-Amplitude_GCC_yearly_mean)/Amplitude_GCC_yearly_mean,
                                                                          Cumulative_GCC_yearly_percent_of_mean = (Cumulative_GCC_yearly-Cumulative_GCC_yearly_mean)/Cumulative_GCC_yearly_mean,
                                                                          SOS_25_difference_from_mean = as.numeric(SOS_25 - as.Date(paste(Year-1,"-12-31",sep="")))-SOS_25_mean)

# ###########################################
# # Calculate drought statistics by site-year
# ###########################################
# 
# site_year_statistics_drought_filtered <- cbind(site_year_statistics_filtered,Cumulative_DM=NA, Cumulative_DM_with_D0=NA,Drought_Weeks=NA)
# 
# for(i in 1:nrow(site_year_statistics_drought_filtered)){
#   phenocam <- site_year_statistics_drought_filtered$Phenocam[i]
#   year <- site_year_statistics_drought_filtered$Year[i]
#   temp_phen_usdm <- filter(phen_usdm_modified_filtered,Phenocam==phenocam,Year==year)
#   cumulative_dm_with_d0 <- (nrow(filter(temp_phen_usdm,DM==0)) + 2*nrow(filter(temp_phen_usdm,DM==1)) +
#     3*nrow(filter(temp_phen_usdm,DM==2)) + 4*nrow(filter(temp_phen_usdm,DM==3)) +
#     5*nrow(filter(temp_phen_usdm,DM==4)))
#   cumulative_dm_with_d0 <- cumulative_dm_with_d0/nrow(temp_phen_usdm)
#   temp_phen_usdm <- filter(phen_usdm_modified_filtered,Phenocam==phenocam,Year==year,DM>=1)
#   cumulative_dm <- sum(temp_phen_usdm$DM,na.rm=TRUE)
#   drought_weeks <- nrow(temp_phen_usdm)
#   site_year_statistics_drought_filtered$Cumulative_DM[i] <- cumulative_dm
#   site_year_statistics_drought_filtered$Cumulative_DM_with_D0[i] <- cumulative_dm_with_d0
#   site_year_statistics_drought_filtered$Drought_Weeks[i] <- drought_weeks
# }
# 
# rm(i,phenocam,year,temp_phen_usdm,cumulative_dm,cumulative_dm_with_d0,drought_weeks)

######################
# Add USDM statistics 
######################

site_year_statistics_drought_filtered <- left_join(site_year_statistics_filtered, usdm_statistics_filtered, by=join_by(Phenocam,Year))
rm(usdm_statistics_filterd)

#####################
# Add VPD statistics 
#####################

site_year_statistics_drought_filtered <- left_join(site_year_statistics_drought_filtered, vpd_statistics_filtered, by=join_by(Phenocam,Year))
rm(vpd_statistics_filtered)

######################
# Add SPEI statistics
######################

site_year_statistics_drought_filtered <- left_join(site_year_statistics_drought_filtered, spei_statistics_filtered, by=join_by(Phenocam,Year))
rm(spei_statistics_filtered)

##################################
# Add ecoregion and aridity index
##################################

site_year_statistics_drought_filtered <- left_join(site_year_statistics_drought_filtered, select(phen_with_ecoregion_df, site, NA_L1CODE, NA_L1NAME, aridity_index), by=join_by(Phenocam==site))
rm(phen_with_ecoregion_df)

site_year_statistics_drought_filtered <- site_year_statistics_drought_filtered %>% mutate(aridity_label = case_when(aridity_index<0.03 ~ "Hyper Arid",
                                                                                                                    aridity_index<0.2 ~ "Arid",
                                                                                                                    aridity_index<0.5 ~ "Semi-Arid",
                                                                                                                    aridity_index<0.65 ~ "Dry Sub-Humid",
                                                                                                                    aridity_index>=0.65 ~ "Humid"))

######################################################
# Calculate site numbers by primary veg type and year
######################################################

site_numbers <- site_year_statistics_drought_filtered %>% count(Primary_Veg_Type,Year, sort=TRUE)
veg_type_ecoregions <- site_year_statistics_drought_filtered %>% count(Primary_Veg_Type, NA_L1NAME)

##########################################################
# Filter by primary veg type for correlation and plotting
##########################################################

temp_data <- site_year_statistics_drought_filtered %>%
  filter(Primary_Veg_Type=="AG" | Primary_Veg_Type=="EN" | Primary_Veg_Type=="DB" | Primary_Veg_Type=="GR" | Primary_Veg_Type=="SH")

#########################
# Linear Model Summaries
#########################

# Results matrices
veg_type_list <- c("AG","DB","EN","GR","SH")

usdm_correlation_summaries <- data.frame()

####################
# Weighted Ave USDM
####################

###########################
# USDM and peak GCC - start
###########################

# Peak GCC vs. Weighted Ave USDM (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Peak GCC vs. Weighted Ave USDM (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Peak GCC vs. Weighted Ave USDM (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Weighted_Ave_Oct_Feb_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_Oct_Feb_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Peak GCC vs. Weighted Ave USDM (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Peak GCC vs. Weighted Ave USDM (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

#########################
# USDM and peak GCC - end
#########################

################################
# USDM and amplitude GCC - start
################################

# GCC Amplitude vs. Weighted Ave USDM (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GCC Amplitude vs. Weighted Ave USDM (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GCC Amplitude vs. Weighted Ave USDM (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Weighted_Ave_Oct_Feb_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_Oct_Feb_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GCC Amplitude vs. Weighted Ave USDM (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GCC Amplitude vs. Weighted Ave USDM (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

##############################
# USDM and amplitude GCC - end
##############################

#################################
# USDM and cumulative GCC - start
#################################

# Cumulative GCC vs. Weighted Ave USDM (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Cumulative GCC vs. Weighted Ave USDM (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Cumulative GCC vs. Weighted Ave USDM (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Weighted_Ave_Oct_Feb_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_Oct_Feb_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Cumulative GCC vs. Weighted Ave USDM (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# Cumulative GCC vs. Weighted Ave USDM (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

###############################
# USDM and cumulative GCC - end
###############################

##############################
# USDM and SOS anomaly - start
##############################

# # SOS difference from mean vs. Weighted Ave USDM (Current year)
# for (veg_type in veg_type_list){
#   temp_model <- lm(SOS_25_difference_from_mean ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
#   new_row <- data.frame("Veg_Type" = veg_type,
#                         "Phen_Metric" = "SOS_25_difference_from_mean",
#                         "Drought_Index" = "USDM",
#                         "Drought_Metric" = "Weighted_Ave_USDM",
#                         "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
#                         "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
#                         "Adjusted_R2" = summary(temp_model)$adj.r.squared,
#                         "F_Statistic" = summary(temp_model)$fstatistic[[1]],
#                         "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
#   usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
# }

# SOS difference from mean vs. Weighted Ave USDM (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# SOS difference from mean vs. Weighted Ave USDM (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Weighted_Ave_Oct_Feb_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_Oct_Feb_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# SOS difference from mean vs. Weighted Ave USDM (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# SOS difference from mean vs. Weighted Ave USDM (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

############################
# USDM and SOS anomaly - end
############################

######################
# USDM and GSL - start
######################

# GSL difference from mean vs. Weighted Ave USDM (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GSL vs. Weighted Ave USDM (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GSL vs. Weighted Ave USDM (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Weighted_Ave_Oct_Feb_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Weighted_Ave_Oct_Feb_USDM",
                        "Slope" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Weighted_Ave_Oct_Feb_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GSL vs. Weighted Ave USDM (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Weighted_Ave_May_Sept_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_May_Sept_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_May_Sept_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

# GSL vs. Weighted Ave USDM (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "USDM",
                        "Drought_Metric" = "Prev_Year_Weighted_Ave_USDM",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Weighted_Ave_USDM","Pr(>|t|)"])
  usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
}

####################
# USDM and GSL - end
####################

# # Relative Yearly GCC Amplitude vs. Cumulative DM, including D0
# for (veg_type in veg_type_list){
#   temp_model <- lm(Amplitude_GCC_yearly_percent_of_mean ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
#   new_row <- data.frame("Veg_Type" = veg_type,
#                         "Phen_Metric" = "Amplitude_GCC_yearly_percent_of_mean",
#                         "Drought_Index" = "USDM",
#                         "Drought_Metric" = "Weighted_Average_DM",
#                         "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
#                         "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
#                         "Adjusted_R2" = summary(temp_model)$adj.r.squared,
#                         "F_Statistic" = summary(temp_model)$fstatistic[[1]],
#                         "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
#   usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
# }
# 
# # Relative Yearly Cumulative GCC vs. Cumulative DM, including D0
# for (veg_type in veg_type_list){
#   temp_model <- lm(Cumulative_GCC_yearly_percent_of_mean ~ Weighted_Ave_USDM, filter(temp_data, Primary_Veg_Type==veg_type))
#   new_row <- data.frame("Veg_Type" = veg_type,
#                         "Phen_Metric" = "Cumulative_GCC_yearly_percent_of_mean",
#                         "Drought_Index" = "USDM",
#                         "Drought_Metric" = "Weighted_Average_DM",
#                         "Slope" = coef(summary(temp_model))["Weighted_Ave_USDM","Estimate"],
#                         "Error" = coef(summary(temp_model))["Weighted_Ave_USDM","Std. Error"],
#                         "Adjusted_R2" = summary(temp_model)$adj.r.squared,
#                         "F_Statistic" = summary(temp_model)$fstatistic[[1]],
#                         "P_Value" = coef(summary(temp_model))["Weighted_Ave_USDM","Pr(>|t|)"])
#   usdm_correlation_summaries <- rbind(usdm_correlation_summaries, new_row)
# }

######
# VPD
######

vpd_correlation_summaries <- data.frame()

##########################
# VPD and peak GCC - start
##########################

# Peak GCC vs. Mean VPD Anomaly (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Peak GCC vs. Mean VPD Anomaly (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Peak GCC vs. Mean VPD Anomaly (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_Oct_Feb_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_Oct_Feb_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Peak GCC vs. Mean VPD Anomaly (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Peak GCC vs. Mean VPD Anomaly (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

########################
# VPD and peak GCC - end
########################

###############################
# VPD and amplitude GCC - start
###############################

# GCC Amplitude vs. Mean VPD Anomaly (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean VPD Anomaly (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean VPD Anomaly (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_Oct_Feb_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_Oct_Feb_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean VPD Anomaly (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean VPD Anomaly (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

#############################
# VPD and amplitude GCC - end
#############################

################################
# VPD and cumulative GCC - start
################################

# Cumulative GCC vs. Mean VPD Anomaly (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean VPD Anomaly (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean VPD Anomaly (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_Oct_Feb_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_Oct_Feb_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean VPD Anomaly (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean VPD Anomaly (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

##############################
# VPD and cumulative GCC - end
##############################

#############################
# VPD and SOS anomaly - start
#############################

# SOS difference from mean vs. Mean VPD Anomaly (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean VPD Anomaly (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean VPD Anomaly (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_Oct_Feb_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_Oct_Feb_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean VPD Anomaly (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean VPD Anomaly (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

###########################
# VPD and SOS anomaly - end
###########################

#####################
# VPD and GSL - start
#####################

# GSL vs. Mean VPD Anomaly (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GSL vs. Mean VPD Anomaly (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GSL vs. Mean VPD Anomaly (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_Oct_Feb_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Mean_Oct_Feb_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GSL vs. Mean VPD Anomaly (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_May_Sept_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

# GSL vs. Mean VPD Anomaly (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_VPD_Anomaly, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "VPD",
                        "Drought_Metric" = "Prev_Year_Mean_VPD_Anomaly",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_VPD_Anomaly","Pr(>|t|)"])
  vpd_correlation_summaries <- rbind(vpd_correlation_summaries, new_row)
}

###################
# VPD and GSL - end
###################

#######
# SPEI
#######

spei_correlation_summaries <- data.frame()

###################################
# 3-month SPEI and peak GCC - start
###################################

# Peak GCC vs. Mean SPEI 3-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 3-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 3-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_Oct_Feb_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 3-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 3-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#################################
# 3-month SPEI and peak GCC - end
#################################

###################################
# 6-month SPEI and peak GCC - start
###################################

# Peak GCC vs. Mean SPEI 6-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 6-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 6-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_Oct_Feb_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 6-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 6-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#################################
# 6-month SPEI and peak GCC - end
#################################

####################################
# 12-month SPEI and peak GCC - start
####################################

# Peak GCC vs. Mean SPEI 12-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 12-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 12-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Mean_Oct_Feb_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 12-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Peak GCC vs. Mean SPEI 12-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Peak_GCC ~ Prev_Year_Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Peak_GCC",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

##################################
# 12-month SPEI and peak GCC - end
##################################

########################################
# 3-month SPEI and amplitude GCC - start
########################################

# GCC Amplitude vs. Mean SPEI 3-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 3-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 3-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_Oct_Feb_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 3-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 3-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

######################################
# 3-month SPEI and amplitude GCC - end
######################################

########################################
# 6-month SPEI and amplitude GCC - start
########################################

# GCC Amplitude vs. Mean SPEI 6-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 6-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 6-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_Oct_Feb_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 6-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 6-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

######################################
# 6-month SPEI and amplitude GCC - end
######################################

#########################################
# 12-month SPEI and amplitude GCC - start
#########################################

# GCC Amplitude vs. Mean SPEI 12-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 12-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 12-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Mean_Oct_Feb_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 12-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GCC Amplitude vs. Mean SPEI 12-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Amplitude_GCC_yearly ~ Prev_Year_Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Amplitude_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#######################################
# 12-month SPEI and amplitude GCC - end
#######################################

#########################################
# 3-month SPEI and cumulative GCC - start
#########################################

# Cumulative GCC vs. Mean SPEI 3-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 3-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 3-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_Oct_Feb_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 3-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 3-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#######################################
# 3-month SPEI and cumulative GCC - end
#######################################

#########################################
# 6-month SPEI and cumulative GCC - start
#########################################

# Cumulative GCC vs. Mean SPEI 6-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 6-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 6-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_Oct_Feb_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 6-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 6-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#######################################
# 6-month SPEI and cumulative GCC - end
#######################################

##########################################
# 12-month SPEI and cumulative GCC - start
##########################################

# Cumulative GCC vs. Mean SPEI 12-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 12-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 12-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Mean_Oct_Feb_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 12-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# Cumulative GCC vs. Mean SPEI 12-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(Cumulative_GCC_yearly ~ Prev_Year_Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "Cumulative_GCC_yearly",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

########################################
# 12-month SPEI and cumulative GCC - end
########################################

######################################
# 3-month SPEI and SOS anomaly - start
######################################

# SOS difference from mean vs. Mean SPEI 3-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 3-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 3-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_Oct_Feb_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 3-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 3-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

####################################
# 3-month SPEI and SOS anomaly - end
####################################

######################################
# 6-month SPEI and SOS anomaly - start
######################################

# SOS difference from mean vs. Mean SPEI 6-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 6-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 6-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_Oct_Feb_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 6-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 6-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

####################################
# 6-month SPEI and SOS anomaly - end
####################################

#######################################
# 12-month SPEI and SOS anomaly - start
#######################################

# SOS difference from mean vs. Mean SPEI 12-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 12-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 12-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Mean_Oct_Feb_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 12-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# SOS difference from mean vs. Mean SPEI 12-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(SOS_25_difference_from_mean ~ Prev_Year_Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "SOS_25_difference_from_mean",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#####################################
# 12-month SPEI and SOS anomaly - end
#####################################

##############################
# 3-month SPEI and GSL - start
##############################

# GSL vs. Mean SPEI 3-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 3-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 3-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_Oct_Feb_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 3-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_May_Sept_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 3-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_SPEI_3month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_3month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_3month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

############################
# 3-month SPEI and GSL - end
############################

##############################
# 6-month SPEI and GSL - start
##############################

# GSL vs. Mean SPEI 6-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 6-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 6-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_Oct_Feb_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 6-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_May_Sept_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 6-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_SPEI_6month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_6month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_6month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

############################
# 6-month SPEI and GSL - end
############################

###############################
# 12-month SPEI and GSL - start
###############################

# GSL vs. Mean SPEI 12-month (Current year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 12-month (Current year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 12-month (Antecedent cooler season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Mean_Oct_Feb_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Mean_Oct_Feb_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Mean_Oct_Feb_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 12-month (Prev year growing season)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_May_Sept_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_May_Sept_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_May_Sept_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

# GSL vs. Mean SPEI 12-month (Prev year)
for (veg_type in veg_type_list){
  temp_model <- lm(GSL_25 ~ Prev_Year_Mean_SPEI_12month, filter(temp_data, Primary_Veg_Type==veg_type))
  new_row <- data.frame("Veg_Type" = veg_type,
                        "Phen_Metric" = "GSL_25",
                        "Drought_Index" = "SPEI",
                        "Drought_Metric" = "Prev_Year_Mean_SPEI_12month",
                        "Slope" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Estimate"],
                        "Error" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Std. Error"],
                        "Adjusted_R2" = summary(temp_model)$adj.r.squared,
                        "F_Statistic" = summary(temp_model)$fstatistic[[1]],
                        "P_Value" = coef(summary(temp_model))["Prev_Year_Mean_SPEI_12month","Pr(>|t|)"])
  spei_correlation_summaries <- rbind(spei_correlation_summaries, new_row)
}

#############################
# 12-month SPEI and GSL - end
#############################

# Add column for significant results
usdm_correlation_summaries <- usdm_correlation_summaries %>% mutate(significant = case_when(P_Value <= 0.05 ~ "yes"))
vpd_correlation_summaries <- vpd_correlation_summaries %>% mutate(significant = case_when(P_Value <= 0.05 ~ "yes"))
spei_correlation_summaries <- spei_correlation_summaries %>% mutate(significant = case_when(P_Value <= 0.05 ~ "yes"))

# Remove superfluous variables
rm(new_row, temp_model, veg_type, veg_type_list)

# Save results as CSV
write.csv(usdm_correlation_summaries,"outputs/usdm_correlation_summaries.csv",row.names=FALSE)
write.csv(vpd_correlation_summaries,"outputs/vpd_correlation_summaries.csv",row.names=FALSE)
write.csv(spei_correlation_summaries,"outputs/spei_correlation_summaries.csv",row.names=FALSE)

###################################################################
# Scatter plots (Phenocam metrics vs. drought metrics by veg type)
###################################################################

######################
# Weighted Average DM
######################

# Peak GCC vs. Weighted Average Yearly USDM
facet_peak <- ggplot(data=temp_data, aes(x=Weighted_Ave_USDM,y=Peak_GCC)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Primary_Veg_Type,ncol=2) +
  ggtitle("Peak GCC vs. Weighted Average Yearly USDM")
facet_peak

# Cumulative GCC vs. Weighted Average Yearly USDM
facet_peak <- ggplot(data=temp_data, aes(x=Weighted_Ave_USDM,y=Cumulative_GCC_yearly)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~Primary_Veg_Type,ncol=2) +
  ggtitle("Cumulative Yearly GCC vs. Weighted Average USDM")
facet_peak

#######################
# Yearly drought plots
#######################

# Load shapefiles
type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
                                     usa_states_vect$NAME != "Hawaii" &
                                     usa_states_vect$NAME != "Puerto Rico",]
level_1_ecoregions_vect <- vect("data/ecoregions/NA_CEC_Eco_Level1.shp")

# Convert phenocam shapefile into dataframe with coordinates
type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")

# Add XY coordinates to site-year and drought statistics
site_year_df <- left_join(temp_data, select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))
drought_statistics_df <- left_join(drought_statistics_filtered, select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))

# Convert site-year and drought statistics into vectors
site_year_vect <- vect(site_year_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
drought_statistics_vect <- vect(drought_statistics_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Only look at sites in the contiguous United States
site_year_vect <- terra::intersect(site_year_vect,project(usa_states_vect,crs(site_year_vect)))
drought_statistics_vect <- terra::intersect(drought_statistics_vect,project(usa_states_vect,crs(drought_statistics_vect)))

# Plot droughts longer than 104 weeks
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Number_Weeks>=104,]), aes(color=Number_Weeks)) +
  labs(color = "Number of Weeks in Drought") +
  scale_color_viridis_c(option = 'viridis') +
  ggtitle("Longest droughts (>=104 weeks) at phenocam locations")
p1

# Plot droughts with >150 cumulative severity
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(drought_statistics_vect[drought_statistics_vect$Sum_USDM>=150,]), aes(color=Sum_USDM)) +
  labs(color = "Cumulative Drought Severity") +
  scale_color_viridis_c(option = 'viridis') +
  ggtitle("Most severe droughts (>=150 cumulative DM) at phenocam locations")
p2

p1/p2

# Plot drought severity by year
# 2016
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2016,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2016 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2016,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2017
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2017,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2017 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2017,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2018
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2018,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2018 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2018,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2019
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2019,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2019 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2019,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2020
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2020,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2020 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2020,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2021
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2021,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2021 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2021,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2022
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2022,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2022 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2022,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2023
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2023,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2023 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2023,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# 2024
p1 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2024,]), aes(color=Weighted_Ave_USDM)) +
  labs(color = 'Weighted Average DM') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 5, by = 1)) +
  ggtitle("2024 Drought at Phenocam Locations")
p2 <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(site_year_vect[site_year_vect$Year==2024,]), aes(color=Drought_Weeks)) +
  labs(color = 'Number of Drought Weeks') +
  scale_color_viridis_c(option = 'viridis') +
  expand_limits(colour = seq(0, 52, by = 13))
p1 / p2

# Plot primary veg type
primary_rois_df <- left_join(primary_rois,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))
primary_rois_df <- filter(primary_rois_df, Primary_Veg_Type=="AG" | Primary_Veg_Type=="EN" | Primary_Veg_Type=="DB" | Primary_Veg_Type=="GR" | Primary_Veg_Type=="SH")

# Convert primary ROIs into vector
primary_rois_vect <- vect(primary_rois_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)

# Only look at ROIs in the contiguous United States
primary_rois_vect <- terra::intersect(primary_rois_vect,project(usa_states_vect,crs(site_year_vect)))

# Plot (all years)
p <- ggplot() +
  geom_sf(data = sf::st_as_sf(usa_states_vect)) +
  geom_sf(data = sf::st_as_sf(primary_rois_vect), aes(color=Primary_Veg_Type)) +
  labs(color = 'Primary Veg Type') +
  ggtitle("PhenoCam locations")
p

###################################
# Other plots (testing with Darren)
###################################

# facet_peak <- ggplot(data=filter(temp_data, Primary_Veg_Type=="AG" | Primary_Veg_Type=="EN" | Primary_Veg_Type=="DB" | Primary_Veg_Type=="GR" | Primary_Veg_Type=="SH"), aes(x=Cumulative_DM,y=Peak_GCC, col = Year)) +
#   geom_point() +
#   geom_smooth(method = loess) +
#   facet_wrap(~Primary_Veg_Type,ncol=2) +
#   ggtitle("Peak GCC vs. Cumulative DM (not including D0)")
# facet_peak
# 
# facet_peak <- ggplot(data=filter(temp_data, Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH"), aes(x=Cumulative_DM_with_D0,y=Peak_GCC, col = as.factor(Year))) +
#   geom_point() +
#   geom_smooth(method = loess) +
#   facet_wrap(~Veg_Type,ncol=2,scales = "free_y") +
#   ggtitle("Peak GCC vs. Weighted Average DM (including D0)")
# facet_peak
# 
# # Plot residuals for the above by veg type
# linear_model_1 <- temp_data %>%
#   filter(Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH") %>%
#   group_by(Veg_Type) %>%
#   do(model = lm(Peak_GCC ~ Cumulative_DM, data = .))
# 
# linear_model_2 <- temp_data %>%
#   filter(Veg_Type=="AG" | Veg_Type=="EN" | Veg_Type=="DB" | Veg_Type=="GR" | Veg_Type=="SH") %>%
#   group_by(Veg_Type) %>%
#   do(model = lm(Peak_GCC ~ Cumulative_DM_with_D0, data = .))
# 
# residuals_df_1 <- data.frame()
# residuals_df_2 <- data.frame()
# 
# for (i in 1:5){
#   veg_type <- linear_model_1$Veg_Type[i]
#   cum_dm_1 <- linear_model_1[[2]][[i]][["model"]][["Cumulative_DM"]]
#   peak_gcc <- linear_model_1[[2]][[i]][["model"]][["Peak_GCC"]]
#   resid_1 <- linear_model_1[[2]][[i]][["residuals"]]
#   temp_df <- data.frame("Veg_Type" = veg_type,
#                         "Cumulative_DM" = cum_dm_1,
#                         "Peak_GCC" = peak_gcc,
#                         "Residuals_no_D0" = resid_1)
#   residuals_df_1 <- rbind(residuals_df_1,temp_df)
# }
# 
# for (i in 1:5){
#   veg_type <- linear_model_2$Veg_Type[i]
#   peak_gcc <- linear_model_2[[2]][[i]][["model"]][["Peak_GCC"]]
#   cum_dm_2 <- linear_model_2[[2]][[i]][["model"]][["Cumulative_DM_with_D0"]]
#   resid_2 <- linear_model_2[[2]][[i]][["residuals"]]
#   temp_df <- data.frame("Veg_Type" = veg_type,
#                         "Cumulative_DM_with_D0" = cum_dm_2,
#                         "Peak_GCC" = peak_gcc,
#                         "Residuals_with_D0" = resid_2)
#   residuals_df_2 <- rbind(residuals_df_2,temp_df)
# }
# 
# rm(cum_dm_1,cum_dm_2,i,peak_gcc,resid_1,resid_2,veg_type)
# 
# facet_resid <- ggplot(data=residuals_df_1, aes(x=Cumulative_DM,y=Residuals_no_D0)) +
#   geom_point() +
#   facet_wrap(~Veg_Type,ncol=2) +
#   ggtitle("Residuals for Peak GCC vs. Cumulative DM (not including D0)")
# facet_resid
# 
# facet_resid <- ggplot(data=residuals_df_2, aes(x=Cumulative_DM_with_D0,y=Residuals_with_D0)) +
#   geom_point() +
#   facet_wrap(~Veg_Type,ncol=2) +
#   ggtitle("Residuals for Peak GCC vs. Weighted Average DM (including D0)")
# facet_resid
# 
# for(veg_type in unique(temp_data$Veg_Type)){
#   p1 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=Cumulative_GCC_yearly)) +
#     geom_point() +
#     ggtitle(paste("Veg type =",veg_type))
#   
#   p2 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=Peak_GCC)) +
#     geom_point()
#   
#   p3 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=GSL_25)) +
#     geom_point()
#   
#   p4 <- ggplot(data=filter(temp_data,Veg_Type==veg_type), aes(x=Cumulative_DM,y=SOS_25_DOY)) +
#     geom_point()
#   
#   print((p1 + p2)/(p3 + p4))
# }
# 
# # Plot cumulative DM vs cumulative and peak GCC by ecoregion
# temp_data <- left_join(temp_data,distinct(select(phen_usdm_modified_filtered,Phenocam,NA_L1NAME)))
# 
# for(ecoregion in unique(temp_data$NA_L1NAME)){
#   p1 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=Cumulative_GCC_yearly)) +
#     geom_point() +
#     ggtitle(paste("Ecoregion =",ecoregion))
#   
#   p2 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=Peak_GCC)) +
#     geom_point()
#   
#   p3 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=GSL_25)) +
#     geom_point()
#   
#   p4 <- ggplot(data=filter(temp_data,NA_L1NAME==ecoregion), aes(x=Cumulative_DM,y=SOS_25_DOY)) +
#     geom_point()
#   
#   print((p1 + p2)/(p3 + p4))
# }

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

# ################################
# # Phenocam Drought Indicator #2:
# # Low Growing Season GCC
# ################################
# 
# # Add columns to growing_seasons_phen_number
# growing_seasons_phen_number <- cbind(growing_seasons_phen_number,
#                                      Peak_GCC_Percent_of_Ave=0, Peak_GCC_STD_from_Ave = 0)
# 
# # Calculate growing season peak GCC as a percentage of the average
# for(i in 1:nrow(growing_seasons_phen_number)){
#   phenocam <- growing_seasons_phen_number$Phenocam[i]
#   veg_type <- growing_seasons_phen_number$Veg_Type[i]
#   roi <- growing_seasons_phen_number$ROI[i]
#   peak_gcc_mean <- filter(growing_seasons_phen_ave, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi)$Peak_GCC_mean
#   peak_gcc_std <- filter(growing_seasons_phen_ave, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi)$Peak_GCC_std
#   
#   growing_seasons_phen_number$Peak_GCC_Percent_of_Ave[i] <- growing_seasons_phen_number$Peak_GCC[i]/peak_gcc_mean
#   growing_seasons_phen_number$Peak_GCC_STD_from_Ave[i] <- (growing_seasons_phen_number$Peak_GCC[i]-peak_gcc_mean)/peak_gcc_std
# }
# 
# rm(i, phenocam, veg_type, roi, peak_gcc_mean, peak_gcc_std)
# 
# # Low peak GCC site-years
# low_gcc <- filter(growing_seasons_phen_number, Peak_GCC_STD_from_Ave <= -1)
# 
# # Add columns for USDM statistics
# low_gcc <- cbind(low_gcc, DM_total_yearly=0, Drought_Weeks_yearly=0,
#                  DM_total_gs=NA, Drought_Weeks_gs=NA, Weeks_gs=NA)
# 
# # Calculate USDM statistics for calendar year
# for(i in 1:nrow(low_gcc)){
#   phenocam <- low_gcc$Phenocam[i]
#   veg_type <- low_gcc$Veg_Type[i]
#   roi <- low_gcc$ROI[i]
#   year <- low_gcc$Year[i]
#   temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
#   low_gcc$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
#   low_gcc$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
#   if(low_gcc$Number_Growing_Seasons[i]==1){
#     sos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$SOS_25
#     eos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$EOS_25
#     temp_gs_df <- filter(temp_df, Date >= sos_25 & Date <= eos_25)
#     
#     low_gcc$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
#     low_gcc$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
#     low_gcc$Weeks_gs[i] <- nrow(temp_gs_df)
#   }
# }
# 
# rm(i, phenocam, roi, veg_type, year, temp_df, temp_gs_df, sos_25, eos_25)
# 
# # Calculate USDM statistics for same sites as above, other years
# normal_gcc <- filter(growing_seasons_phen_number, Peak_GCC_STD_from_Ave > -1)
# normal_gcc <- right_join(normal_gcc, distinct(select(low_gcc,Phenocam,Veg_Type,ROI)), join_by(Phenocam,Veg_Type,ROI))
# 
# # Add columns for USDM statistics
# normal_gcc <- cbind(normal_gcc, DM_total_yearly=0, Drought_Weeks_yearly=0,
#                     DM_total_gs=NA, Drought_Weeks_gs=NA, Weeks_gs=NA)
# 
# # Calculate USDM statistics for calendar year
# for(i in 1:nrow(normal_gcc)){
#   phenocam <- normal_gcc$Phenocam[i]
#   veg_type <- normal_gcc$Veg_Type[i]
#   roi <- normal_gcc$ROI[i]
#   year <- normal_gcc$Year[i]
#   temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year==year)
#   normal_gcc$DM_total_yearly[i] <- sum(filter(temp_df, DM>=1)$DM)
#   normal_gcc$Drought_Weeks_yearly[i] <- nrow(filter(temp_df, DM>=1))
#   
#   if(normal_gcc$Number_Growing_Seasons[i]==1){
#     sos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$SOS_25
#     eos_25 <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$EOS_25
#     temp_gs_df <- filter(temp_df, Date >= sos_25 & Date <= eos_25)
#     
#     normal_gcc$DM_total_gs[i] <- sum(filter(temp_gs_df, DM>=1)$DM)
#     normal_gcc$Drought_Weeks_gs[i] <- nrow(filter(temp_gs_df, DM>=1))
#     normal_gcc$Weeks_gs[i] <- nrow(temp_gs_df)
#   }
# }
# 
# rm(i, phenocam, veg_type, roi, year, temp_df, temp_gs_df, sos_25, eos_25)

#########################
# Plot PhenoCam Locations
#########################

# # Load shapefiles
# type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
# usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
# usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska" &
#                                      usa_states_vect$NAME != "Hawaii" &
#                                      usa_states_vect$NAME != "Puerto Rico",]
# level_1_ecoregions_vect <- vect("data/ecoregions/NA_CEC_Eco_Level1.shp")
# 
# # Convert phenocam shapefile into dataframe with coordinates
# type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")
# 
# # Add XY coordinates to primary ROIs
# primary_rois_df <- left_join(primary_rois,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))
# primary_rois_df <- filter(primary_rois_df, Primary_Veg_Type=="AG" | Primary_Veg_Type=="EN" | Primary_Veg_Type=="DB" | Primary_Veg_Type=="GR" | Primary_Veg_Type=="SH")
# 
# # Convert primary ROIs into vector
# primary_rois_vect <- vect(primary_rois_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# 
# # Convert primary ROIs into sf
# primary_rois_sf <- sf::st_as_sf(primary_rois_vect)
# 
# # Plot primary ROIs (terra)
# plot(usa_states_vect)
# plot(primary_rois_vect,2,add=TRUE)
# 
# # Plot primary ROIs (ggplot)
# ggplot() +
#   #geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = primary_rois_sf, aes(color=Primary_Veg_Type)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   #scale_color_stepsn(colors=c("purple","blue","green","yellow","orange","red","brown"),
#   #                   breaks=c(0,10,20,30,40,50,100,200),
#   #                   guide = "colorsteps") +
#   labs(color = 'Primary Veg Type') +
#   ggtitle("PhenoCam locations")

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
# 
# # Add XY coordinates to drought statistics
# drought_statistics_df <- left_join(drought_statistics_filtered,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))
# 
# # Select top 10, 20, 50, and 100 droughts by cumulative severity
# top_10 <- drought_statistics_df %>% slice_max(Sum_USDM, n=10)
# top_20 <- drought_statistics_df %>% slice_max(Sum_USDM, n=20)
# top_50 <- drought_statistics_df %>% slice_max(Sum_USDM, n=50)
# top_100 <- drought_statistics_df %>% slice_max(Sum_USDM, n=100)
# top_100 <- left_join(top_100, unique(select(temp_data,Phenocam,Primary_Veg_Type,NA_L1NAME)))
# 
# # Convert drought statistics into vector
# drought_statistics_vect <- vect(drought_statistics_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# top_20_vect <- vect(top_20, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# top_50_vect <- vect(top_50, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# top_100_vect <- vect(top_100, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# 
# # Convert drought statistics into sf
# drought_statistics_sf <- sf::st_as_sf(drought_statistics_vect)
# top_20_sf <- sf::st_as_sf(top_20_vect)
# top_50_sf <- sf::st_as_sf(top_50_vect)
# top_100_sf <- sf::st_as_sf(top_100_vect)
# 
# # Plot drought statistics (terra)
# plot(usa_states_vect)
# plot(drought_statistics_vect,6,add=TRUE) # Number of weeks of drought
# plot(usa_states_vect)
# plot(drought_statistics_vect,7,add=TRUE) # Sum of drought monitor index
# 
# # Plot drought statistics (ggplot)
# ggplot(data = drought_statistics_sf, aes(color=Number_Weeks)) +
#   #scale_color_distiller(palette = 'YlOrRd', direction = 1) +
#   scale_color_viridis_c(option = 'viridis') +
#   geom_sf() +
#   labs(color = 'Drought length (weeks)') +
#   ggtitle("Number of Weeks of Drought")
# 
# ggplot(data = drought_statistics_sf, aes(color=Sum_USDM)) +
#   #scale_color_distiller(palette = 'YlOrRd', direction = 1) +
#   scale_color_viridis_c(option = 'viridis') +
#   geom_sf() +
#   labs(color = 'Sum of drought index') +
#   ggtitle("Cumulative drought severity")
# 
# ggplot(data = top_20_sf, aes(color=Sum_USDM)) +
#   #scale_color_distiller(palette = 'YlOrRd', direction = 1) +
#   scale_color_viridis_c(option = 'viridis') +
#   geom_sf() +
#   labs(color = 'Sum of drought index') +
#   ggtitle("Cumulative drought severity")
# 
# # Plot drought statistics (ggplot, contiguous USA)
# drought_statistics_contig_vect <- terra::intersect(drought_statistics_vect,project(usa_states_vect,crs(drought_statistics_vect)))
# drought_statistics_contig_sf <- sf::st_as_sf(drought_statistics_contig_vect)
# drought_statistics_contig_sf <- mutate(drought_statistics_contig_sf, Number_Weeks_Binned = cut(Number_Weeks, breaks=c(0,10,20,40,60,100,200)))
# drought_statistics_contig_sf <- mutate(drought_statistics_contig_sf, Sum_USDM_Binned = cut(Sum_USDM, breaks=c(0,25,50,100,200,400,600)))
# 
# ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = drought_statistics_contig_sf, aes(color=Number_Weeks_Binned)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   #scale_color_stepsn(colors=c("purple","blue","green","yellow","orange","red","brown"),
#   #                   breaks=c(0,10,20,30,40,50,100,200),
#   #                   guide = "colorsteps") +
#   labs(color = 'Drought length (weeks)') +
#   ggtitle("Number of Weeks of Continuous Drought")
# 
# ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = drought_statistics_contig_sf, aes(color=Sum_USDM_Binned)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   labs(color = 'Sum of drought index') +
#   ggtitle("Cumulative Severity (DM) of Continuous Drought")
# 
# ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = top_50_sf, aes(color=Sum_USDM)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   labs(color = 'Sum of drought index') +
#   ggtitle("Cumulative Severity (DM) of Continuous Drought")
# 
# ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = top_100_sf, aes(color=NA_L1NAME)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   labs(color = 'Ecoregion name') +
#   ggtitle("Top 100 droughts by cumulative severity, 2016-2024")
# 
# ggplot() +
#   geom_sf(data = sf::st_as_sf(usa_states_vect)) +
#   geom_sf(data = top_100_sf, aes(shape=Primary_Veg_Type)) +
#   #scale_color_viridis_c(option = 'viridis') +
#   labs(shape = 'Primary Veg Type') +
#   ggtitle("Top 100 USA droughts, measured by cumulative severity")

# #######################
# # Plot no-gs site-years
# #######################
# 
# # Load shapefiles
# # type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
# # usa_states_vect <- vect("data/shapefiles/cb_2018_us_state_20m.shp")
# # usa_states_vect <- usa_states_vect[usa_states_vect$NAME != "Alaska"]
# 
# # Convert phenocam shapefile into dataframe with coordinates
# # type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect,geom="XY")
# 
# # Add XY coordinates to no_gs
# no_gs_df <- left_join(no_gs,select(type1_sites_usa_df,site,x,y),join_by(Phenocam==site))
# 
# # Convert drought statistics into vector
# no_gs_vect <- vect(no_gs_df, geom=c("x", "y"), crs=crs(type1_sites_usa_vect), keepgeom=FALSE)
# 
# # Plot no gs site-years
# plot(usa_states_vect)
# plot(no_gs_vect,2,add=TRUE) # Veg Type
# plot(no_gs_vect,2) # Veg Type

# ##########################################
# # Plot SOS vs. EOS with standard deviation
# ##########################################
# 
# # plot of SOS vs EOS, colored by veg type
# ggplot(data = growing_seasons_phen_ave_filtered, aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   ggtitle("All Primary ROIs")
# 
# # plot without error bars, filtered so that max stdev is no more than 20 days
# ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=20 & EOS_25_std<=20), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   ggtitle("Phenocam ROIs with SOS and EOS stdev <20 days")
# 
# # plot of SOS vs EOS, with error bars
# ggplot(data = growing_seasons_phen_ave_filtered, aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
#   geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
#   ggtitle("All Phenocam ROIs")
# 
# # plot with error bars, but filtered so that max stdev is no more than 30 days
# ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=30 & EOS_25_std<=30), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
#   geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
#   ggtitle("SOS and EOS stdev <30 days")
# 
# # plot with error bars, but filtered so that max stdev is no more than 20 days
# ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=20 & EOS_25_std<=20), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
#   geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
#   ggtitle("SOS and EOS stdev <20 days")
# 
# # plot with error bars, but filtered so that max stdev is no more than 10 days
# ggplot(data = filter(growing_seasons_phen_ave_filtered, SOS_25_std<=10 & EOS_25_std<=10), aes(x=SOS_25_mean, y=EOS_25_mean, color=Veg_Type)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=EOS_25_mean-EOS_25_std, ymax=EOS_25_mean+EOS_25_std)) +
#   geom_errorbarh(aes(xmin=SOS_25_mean-SOS_25_std, xmax=SOS_25_mean+SOS_25_std)) +
#   ggtitle("SOS and EOS stdev <10 days")

# ####################
# # Preliminary plots
# ####################
# 
# drought_indicator_1 <- right_join(select(no_gs, Phenocam, Veg_Type, ROI, Year, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
#                                   select(gs_df, Phenocam, Veg_Type, ROI, Year, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
#                                   join_by(Phenocam,Veg_Type,ROI))
# 
# p1 <- ggplot(data = drought_indicator_1, aes(x = DM_total_yearly.x, y = DM_total_yearly.y)) +
#   geom_point() +
#   xlab("Total DM for no growing season year") +
#   ylab("Total DM for growing season years, same site")
# p2 <- ggplot(data = drought_indicator_1, aes(x = Drought_Weeks_yearly.x, y = Drought_Weeks_yearly.y)) +
#   geom_point() +
#   xlab("# Weeks in Drought for no growing season year") +
#   ylab("# Weeks in Drought for growing season years, same site")
# p1
# p2
# 
# drought_indicator_2 <- right_join(select(low_gcc, Phenocam, Veg_Type, ROI, Year, Peak_GCC_Percent_of_Ave, Peak_GCC_STD_from_Ave, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
#                                   select(normal_gcc, Phenocam, Veg_Type, ROI, Year, Peak_GCC_Percent_of_Ave, Peak_GCC_STD_from_Ave, DM_total_yearly, Drought_Weeks_yearly, DM_total_gs, Drought_Weeks_gs, Weeks_gs),
#                                   join_by(Phenocam,Veg_Type,ROI), relationship = "many-to-many")
# 
# p3 <- ggplot(data = drought_indicator_2, aes(x = DM_total_yearly.x, y = DM_total_yearly.y)) +
#   geom_point() +
#   xlab("Total DM for low GCC year") +
#   ylab("Total DM for normal GCC years, same site")
# p3
# 
# p4 <- ggplot() +
#   geom_point(data = low_gcc, aes(x = Peak_GCC_STD_from_Ave, y = DM_total_yearly), color = "red") +
#   geom_point(data = normal_gcc, aes(x = Peak_GCC_STD_from_Ave, y = DM_total_yearly), color = "black") +
#   xlim(-10,10) +
#   xlab("Standard deviations of peak GCC away from mean") +
#   ylab("Total DM for site-year")
# p4
# 
# p5 <- ggplot() +
#   geom_point(data = low_gcc, aes(x = Peak_GCC_STD_from_Ave, y = Drought_Weeks_yearly), color = "red") +
#   geom_point(data = normal_gcc, aes(x = Peak_GCC_STD_from_Ave, y = Drought_Weeks_yearly), color = "black") +
#   xlim(-10,10) +
#   xlab("Standard deviations of peak GCC away from mean") +
#   ylab("# Weeks in drought for site-year")
# p5

