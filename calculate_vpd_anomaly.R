# Load in VPD from CSV, calculate VPD anomaly

library(dplyr)
library(lubridate)

vpd_df <- data.frame()

# Read in VPD data
for (filename in list.files("data/vpd", full.names=TRUE)){
  temp_df <- read.csv(filename)
  vpd_df <- rbind(vpd_df, temp_df)
}

rm(temp_df, filename)

# Rename variables and remove invalid values
vpd_df <- vpd_df %>%
  rename(VPD = values, Phenocam = site, Date = date) %>%
  mutate(Date = as.Date(Date), Year=year(Date), Month=month(Date), DOY=yday(Date)) %>%
  filter(VPD >= 0)

# Calculate mean, peak, and cumulative VPD by site-year
vpd_statistics <- data.frame()

unique_siteyears <- vpd_df %>% distinct(Phenocam,Year)
for (i in 1:nrow(unique_siteyears)){
  site <- unique_siteyears$Phenocam[i]
  year <- unique_siteyears$Year[i]
  print(paste("i = ",i,", site = ",site,", year = ",year,sep=""))
  temp_df <- filter(vpd_df, Phenocam==site & Year==year)
  temp_gs_df <- filter(temp_df, Month>=5, Month<=9)
  new_row <- data.frame("Phenocam" = site,
                        "Year" = year,
                        "Max_VPD" = max(temp_df$VPD, na.rm=TRUE),
                        "Mean_VPD" = mean(temp_df$VPD, na.rm=TRUE),
                        "Mean_May_Sept_VPD" = mean(temp_gs_df$VPD, na.rm=TRUE),
                        "Cumulative_VPD" = sum(temp_df$VPD, na.rm=TRUE),
                        "Number_of_Days" = nrow(temp_df))
  vpd_statistics <- rbind(vpd_statistics, new_row)
}

rm(new_row, temp_df, temp_gs_df, unique_siteyears, i, site, year)

# Remove years with <365 days of data from cumulative VPD calculation
vpd_statistics <- vpd_statistics %>% mutate(Cumulative_VPD=replace(Cumulative_VPD, Number_of_Days<365, NA))

# Mean statistics
vpd_summary <- data.frame()

for (i in 1:nrow(distinct(vpd_statistics, Phenocam))){
  site <- vpd_statistics$Phenocam[i]
  temp_df <- filter(vpd_statistics, Phenocam==site)
  new_row <- data.frame("Phenocam" = site,
                        "Average_Max_VPD" = mean(temp_df$Max_VPD, na.rm=TRUE),
                        "Average_Mean_VPD" = mean(temp_df$Mean_VPD, na.rm=TRUE),
                        "Average_Mean_May_Sept_VPD" = mean(temp_df$Mean_May_Sept_VPD, na.rm=TRUE),
                        "Average_Cumulative_VPD" = mean(temp_df$Cumulative_VPD, na.rm=TRUE))
  vpd_summary <- rbind(vpd_summary, new_row)
}

rm(i, site, temp_df, new_row)

# Calculate VPD anomalies
vpd_statistics <- left_join(vpd_statistics,vpd_summary)

vpd_statistics <- mutate(vpd_statistics,
                         Max_VPD_Anomaly = Max_VPD-Average_Max_VPD,
                         Mean_VPD_Anomaly = Mean_VPD-Average_Mean_VPD,
                         Mean_May_Sept_VPD_Anomaly = Mean_May_Sept_VPD-Average_Mean_May_Sept_VPD,
                         Cumulative_VPD_Anomaly = Cumulative_VPD-Average_Cumulative_VPD)

# Save results
save(vpd_statistics, file="outputs/vpd_statistics.RData")