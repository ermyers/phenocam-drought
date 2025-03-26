# Load in SPEI from CSV, calculate SPEI anomalies

# Note 3-25-2025: Modified from original.
# Added: Standard deviation (per site-year), number of drought weeks (per site-year)
# Removed: Maximum (per site-year), average (per site), anomaly (per site-year,
# compared to per-site average)

library(dplyr)
library(lubridate)

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

# Calculate mean, peak, and cumulative SPEI by site-year
spei_statistics <- data.frame()

unique_siteyears <- spei_df %>% distinct(Phenocam,Year)
for (i in 1:nrow(unique_siteyears)){
  site <- unique_siteyears$Phenocam[i]
  year <- unique_siteyears$Year[i]
  print(paste("i = ",i,", site = ",site,", year = ",year,sep=""))
  temp_df <- filter(spei_df, Phenocam==site & Year==year)
  new_row <- data.frame("Phenocam" = site,
                        "Year" = year,
                        "Mean_SPEI_1month" = mean(temp_df$SPEI_1month, na.rm=TRUE),
                        "Mean_SPEI_3month" = mean(temp_df$SPEI_3month, na.rm=TRUE),
                        "Mean_SPEI_6month" = mean(temp_df$SPEI_6month, na.rm=TRUE),
                        "Mean_SPEI_12month" = mean(temp_df$SPEI_12month, na.rm=TRUE),
                        "Stdev_SPEI_1month" = sd(temp_df$SPEI_1month, na.rm=TRUE),
                        "Stdev_SPEI_3month" = sd(temp_df$SPEI_3month, na.rm=TRUE),
                        "Stdev_SPEI_6month" = sd(temp_df$SPEI_6month, na.rm=TRUE),
                        "Stdev_SPEI_12month" = sd(temp_df$SPEI_12month, na.rm=TRUE),
                        "Cumulative_SPEI_1month" = sum(temp_df$SPEI_1month, na.rm=TRUE),
                        "Cumulative_SPEI_3month" = sum(temp_df$SPEI_3month, na.rm=TRUE),
                        "Cumulative_SPEI_6month" = sum(temp_df$SPEI_6month, na.rm=TRUE),
                        "Cumulative_SPEI_12month" = sum(temp_df$SPEI_12month, na.rm=TRUE),
                        "Number_of_SPEI_Weeks" = nrow(temp_df),
                        "SPEI_Drought_Weeks_1month" = nrow(filter(temp_df, SPEI_1month<=-1)),
                        "SPEI_Drought_Weeks_3month" = nrow(filter(temp_df, SPEI_3month<=-1)),
                        "SPEI_Drought_Weeks_6month" = nrow(filter(temp_df, SPEI_6month<=-1)),
                        "SPEI_Drought_Weeks_12month" = nrow(filter(temp_df, SPEI_12month<=-1)))
  spei_statistics <- rbind(spei_statistics, new_row)
}

rm(temp_df, unique_siteyears, i, site, year)

# Remove years with <52 weeks of data from cumulative SPEI calculation
spei_statistics <- spei_statistics %>%
  mutate(Cumulative_SPEI_1month=replace(Cumulative_SPEI_1month, Number_of_SPEI_Weeks<52, NA)) %>%
  mutate(Cumulative_SPEI_3month=replace(Cumulative_SPEI_3month, Number_of_SPEI_Weeks<52, NA)) %>%
  mutate(Cumulative_SPEI_6month=replace(Cumulative_SPEI_6month, Number_of_SPEI_Weeks<52, NA)) %>%
  mutate(Cumulative_SPEI_12month=replace(Cumulative_SPEI_12month, Number_of_SPEI_Weeks<52, NA))

# # Mean statistics
# spei_summary <- data.frame()
# 
# for (i in 1:nrow(distinct(spei_statistics, Phenocam))){
#   site <- spei_statistics$Phenocam[i]
#   temp_df <- filter(spei_statistics, Phenocam==site)
#   new_row <- data.frame("Phenocam" = site,
#                         "Average_Max_SPEI_1month" = mean(temp_df$Max_SPEI_1month, na.rm=TRUE),
#                         "Average_Mean_SPEI_1month" = mean(temp_df$Mean_SPEI_1month, na.rm=TRUE),
#                         "Average_Cumulative_SPEI_1month" = mean(temp_df$Cumulative_SPEI_1month, na.rm=TRUE),
#                         "Average_Max_SPEI_3month" = mean(temp_df$Max_SPEI_3month, na.rm=TRUE),
#                         "Average_Mean_SPEI_3month" = mean(temp_df$Mean_SPEI_3month, na.rm=TRUE),
#                         "Average_Cumulative_SPEI_3month" = mean(temp_df$Cumulative_SPEI_3month, na.rm=TRUE),
#                         "Average_Max_SPEI_6month" = mean(temp_df$Max_SPEI_6month, na.rm=TRUE),
#                         "Average_Mean_SPEI_6month" = mean(temp_df$Mean_SPEI_6month, na.rm=TRUE),
#                         "Average_Cumulative_SPEI_6month" = mean(temp_df$Cumulative_SPEI_6month, na.rm=TRUE),
#                         "Average_Max_SPEI_12month" = mean(temp_df$Max_SPEI_12month, na.rm=TRUE),
#                         "Average_Mean_SPEI_12month" = mean(temp_df$Mean_SPEI_12month, na.rm=TRUE),
#                         "Average_Cumulative_SPEI_12month" = mean(temp_df$Cumulative_SPEI_12month, na.rm=TRUE))
#   spei_summary <- rbind(spei_summary, new_row)
# }
# 
# rm(i, site, temp_df, new_row)
# 
# # Calculate SPEI anomalies
# spei_statistics <- left_join(spei_statistics,spei_summary)
# 
# spei_statistics <- mutate(spei_statistics,
#                          Max_SPEI_1month_Anomaly = Max_SPEI_1month-Average_Max_SPEI_1month,
#                          Mean_SPEI_1month_Anomaly = Mean_SPEI_1month-Average_Mean_SPEI_1month,
#                          Cumulative_SPEI_1month_Anomaly = Cumulative_SPEI_1month-Average_Cumulative_SPEI_1month,
#                          Max_SPEI_3month_Anomaly = Max_SPEI_3month-Average_Max_SPEI_3month,
#                          Mean_SPEI_3month_Anomaly = Mean_SPEI_3month-Average_Mean_SPEI_3month,
#                          Cumulative_SPEI_3month_Anomaly = Cumulative_SPEI_3month-Average_Cumulative_SPEI_3month,
#                          Max_SPEI_6month_Anomaly = Max_SPEI_6month-Average_Max_SPEI_6month,
#                          Mean_SPEI_6month_Anomaly = Mean_SPEI_6month-Average_Mean_SPEI_6month,
#                          Cumulative_SPEI_6month_Anomaly = Cumulative_SPEI_6month-Average_Cumulative_SPEI_6month,
#                          Max_SPEI_12month_Anomaly = Max_SPEI_12month-Average_Max_SPEI_12month,
#                          Mean_SPEI_12month_Anomaly = Mean_SPEI_12month-Average_Mean_SPEI_12month,
#                          Cumulative_SPEI_12month_Anomaly = Cumulative_SPEI_12month-Average_Cumulative_SPEI_12month)

# Save results
save(spei_statistics, file="outputs/spei_statistics.RData")
