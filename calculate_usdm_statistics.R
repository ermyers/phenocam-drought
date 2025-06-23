# Load in USDM data, calculate drought statistics

library(dplyr)
library(lubridate)

# Load USDM data
load("outputs/phen_usdm_2010.RData")
load("outputs/phen_usdm_2011.RData")
load("outputs/phen_usdm_2012.RData")
load("outputs/phen_usdm_2013.RData")
load("outputs/phen_usdm_2014.RData")
load("outputs/phen_usdm_2015.RData")
load("outputs/phen_usdm_2016.RData")
load("outputs/phen_usdm_2017.RData")
load("outputs/phen_usdm_2018.RData")
load("outputs/phen_usdm_2019.RData")
load("outputs/phen_usdm_2020.RData")
load("outputs/phen_usdm_2021.RData")
load("outputs/phen_usdm_2022.RData")
load("outputs/phen_usdm_2023.RData")
load("outputs/phen_usdm_2024.RData")

usdm_df <- rbind(phen_usdm_2010,phen_usdm_2011,phen_usdm_2012,phen_usdm_2013,
                 phen_usdm_2014,phen_usdm_2015,phen_usdm_2016,phen_usdm_2017,
                 phen_usdm_2018,phen_usdm_2019,phen_usdm_2020,phen_usdm_2021,
                 phen_usdm_2022,phen_usdm_2023,phen_usdm_2024)

usdm_df[is.na(usdm_df$DM),]$DM <- -1

usdm_df <- usdm_df %>%
  mutate(Month=month(Date)) %>%
  select(Phenocam,Date,Year,Month,DOY,DM)

# Optional - remove yearly USDM
rm(phen_usdm_2010,phen_usdm_2011,phen_usdm_2012,phen_usdm_2013,phen_usdm_2014,
   phen_usdm_2015,phen_usdm_2016,phen_usdm_2017,phen_usdm_2018,phen_usdm_2019,
   phen_usdm_2020,phen_usdm_2021,phen_usdm_2022,phen_usdm_2023,phen_usdm_2024)

# Calculate drought statistics by site-year

usdm_statistics <- data.frame()

unique_siteyears <- usdm_df %>% distinct(Phenocam,Year)
for (i in 1:nrow(unique_siteyears)){
  site <- unique_siteyears$Phenocam[i]
  year <- unique_siteyears$Year[i]
  print(paste("i = ",i,", site = ",site,", year = ",year,sep=""))
  temp_df <- filter(usdm_df, Phenocam==site & Year==year)
  # whole calendar year
  weighted_ave_usdm <- (nrow(filter(temp_df,DM==0)) + 2*nrow(filter(temp_df,DM==1)) +
                          3*nrow(filter(temp_df,DM==2)) + 4*nrow(filter(temp_df,DM==3)) +
                          5*nrow(filter(temp_df,DM==4)))/nrow(temp_df)
  temp_gs_df <- filter(temp_df, Month>=5, Month<=9)
  # growing season (May-September)
  weighted_ave_gs_usdm <- (nrow(filter(temp_gs_df,DM==0)) + 2*nrow(filter(temp_gs_df,DM==1)) +
                             3*nrow(filter(temp_gs_df,DM==2)) + 4*nrow(filter(temp_gs_df,DM==3)) +
                             5*nrow(filter(temp_gs_df,DM==4)))/nrow(temp_gs_df)
  temp_cs_df <- rbind(filter(usdm_df, Phenocam==site, Year==(year-1), Month>=10),
                      filter(usdm_df, Phenocam==site, Year==year, Month<=2))
  # antecedent cooler season (October-February)
  weighted_ave_cs_usdm <- (nrow(filter(temp_cs_df,DM==0)) + 2*nrow(filter(temp_cs_df,DM==1)) +
                             3*nrow(filter(temp_cs_df,DM==2)) + 4*nrow(filter(temp_cs_df,DM==3)) +
                             5*nrow(filter(temp_cs_df,DM==4)))/nrow(temp_cs_df)
  new_row <- data.frame("Phenocam" = site,
                        "Year" = year,
                        "Weighted_Ave_USDM" = weighted_ave_usdm,
                        "Weighted_Ave_May_Sept_USDM" = weighted_ave_gs_usdm,
                        "Weighted_Ave_Oct_Feb_USDM" = weighted_ave_cs_usdm,
                        "Drought_Weeks" = nrow(filter(temp_df, DM>=1)),
                        "Drought_Weeks_May_Sept" = nrow(filter(temp_gs_df, DM>=1)),
                        "Drought_Weeks_Oct_Feb" = nrow(filter(temp_cs_df, DM>=1)))
  usdm_statistics <- rbind(usdm_statistics, new_row)
}

rm(new_row, temp_df, temp_cs_df, temp_gs_df, unique_siteyears, i, site, 
   weighted_ave_cs_usdm, weighted_ave_gs_usdm, weighted_ave_usdm, year)

# Calculate drought statistics for previous year

usdm_statistics <- cbind(usdm_statistics,
                         "Prev_Year_Weighted_Ave_USDM" = NA,
                         "Prev_Year_Weighted_Ave_May_Sept_USDM" = NA)

for (i in 1:nrow(usdm_statistics)){
  site <- usdm_statistics$Phenocam[i]
  year <- usdm_statistics$Year[i]
  prev_year_df <- filter(usdm_statistics, Phenocam==site, Year==(year-1))
  if(nrow(prev_year_df)==1){
    usdm_statistics$Prev_Year_Weighted_Ave_USDM[i] <- prev_year_df$Weighted_Ave_USDM
    usdm_statistics$Prev_Year_Weighted_Ave_May_Sept_USDM[i] <- prev_year_df$Weighted_Ave_May_Sept_USDM
  }
}

rm(prev_year_df, i, site, year)

# Save results
save(usdm_statistics, file="outputs/usdm_statistics.RData")
