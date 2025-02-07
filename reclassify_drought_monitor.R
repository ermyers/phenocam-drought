# Reclassify drought monitor data

library(tidyverse)
library(patchwork)

# Load data
load("outputs/phen_gcc.RData")
#load("outputs/vci_weekly.RData")
load("outputs/phen_usdm.RData")
load("outputs/growing_seasons_phen.RData")

##########
# For VCI
##########

# # Sort DM into binary categories:
# # -1-0: no drought -> 0
# # 1-4: drought -> 1
# 
# vci_weekly <- cbind(vci_weekly,USDM_reclassified=NA)
# vci_weekly[vci_weekly$USDM_weekly<=0,]$USDM_reclassified <- 0 # not drought
# vci_weekly[vci_weekly$USDM_weekly>=1,]$USDM_reclassified <- 1 # drought
# 
# # Record whether DM is increasing, decreasing, or same
# 
# vci_weekly_modified <- data.frame()
# unique_phenocams <- vci_weekly %>% distinct(Phenocam,Veg_Type,ROI)
# 
# for (i in 1:nrow(unique_phenocams)){
#   phenocam <- unique_phenocams$Phenocam[i]
#   veg_type <- unique_phenocams$Veg_Type[i]
#   roi <- unique_phenocams$ROI[i]
#   temp_df <- vci_weekly %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
#   temp_df <- cbind(temp_df, full_USDM_change=NA, reclass_USDM_change=NA)
#   for (j in 2:nrow(temp_df)){
#     temp_df$full_USDM_change[j] <- temp_df$USDM_weekly[j] - temp_df$USDM_weekly[j-1]
#     temp_df$reclass_USDM_change[j] <- temp_df$USDM_reclassified[j] - temp_df$USDM_reclassified[j-1]
#   }
#   vci_weekly_modified <- rbind(vci_weekly_modified,temp_df)
# }
# 
# # Optional
# rm(temp_df,phenocam,veg_type,roi,i,j)

###########
# For USDM
###########

# Sort DM into binary categories:
# -1-0: no drought -> 0
# 1-4: drought -> 1

phen_usdm <- cbind(phen_usdm,USDM_reclassified=NA)
phen_usdm[phen_usdm$DM<=0,]$USDM_reclassified <- 0 # not drought
phen_usdm[phen_usdm$DM>=1,]$USDM_reclassified <- 1 # drought

# Record whether DM is increasing, decreasing, or same

phen_usdm_modified <- data.frame()
unique_phenocams <- phen_usdm %>% distinct(Phenocam)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  temp_df <- phen_usdm %>% filter(Phenocam == phenocam)
  temp_df <- cbind(temp_df, full_USDM_change=NA, reclass_USDM_change=NA)
  for (j in 2:nrow(temp_df)){
    temp_df$full_USDM_change[j] <- temp_df$DM[j] - temp_df$DM[j-1]
    temp_df$reclass_USDM_change[j] <- temp_df$USDM_reclassified[j] - temp_df$USDM_reclassified[j-1]
  }
  phen_usdm_modified <- rbind(phen_usdm_modified,temp_df)
}

# Optional
rm(temp_df,phenocam,i,j)

# Remove drought periods with no corresponding phenocam data
phen_usdm_modified_filtered <- data.frame()
growing_seasons <- distinct(growing_seasons_phen_number, Phenocam, Year)
unique_phenocams <- growing_seasons %>% distinct(Phenocam)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  year_min <- min(filter(growing_seasons, Phenocam==phenocam)$Year, na.rm=TRUE)
  year_max <- max(filter(growing_seasons, Phenocam==phenocam)$Year, na.rm=TRUE)
  temp_df <- filter(phen_usdm_modified, Phenocam==phenocam & Year>=year_min & Year<=year_max)
  phen_usdm_modified_filtered <- rbind(phen_usdm_modified_filtered, temp_df)
}

#phen_usdm_modified <- right_join(phen_usdm_modified,growing_seasons,by=c("Phenocam","Year"))
phen_usdm_modified <- phen_usdm_modified_filtered

# Optional
rm(i, phenocam, year_min, year_max, temp_df, phen_usdm_modified_filtered)

# Record drought period statistics
# (site, number of consecutive weeks, start date and year, end date and year, sum)
drought_statistics <- data.frame(matrix(ncol=7,nrow=0))
colnames(drought_statistics) <- c("Phenocam","Start_Date","Start_Year","End_Date","End_Year","Number_Weeks","Sum_USDM")

unique_phenocams <- phen_usdm_modified %>% distinct(Phenocam)
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  temp_df <- phen_usdm_modified %>% filter(Phenocam == phenocam)
  temp_drought_change <- temp_df %>% filter(reclass_USDM_change != 0)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,sep=""))
  
  # Check whether there are entries in dataframe
  if (nrow(temp_drought_change)>0){
    # Check whether first date is a drought end
    if (temp_drought_change$reclass_USDM_change[1]==-1){
      # If first row is drought end, add a new row with drought end date and no other statistics
      new_row <- data.frame("Phenocam" = phenocam, "Start_Date" = NA, "Start_Year" = NA,
                            "End_Date" = temp_drought_change$Date[1], "End_Year" = temp_drought_change$Year[1],
                            "Number_Weeks" = NA, "Sum_USDM" = NA)
      drought_statistics <- rbind(drought_statistics,new_row)
      
      # Remove first date from temp_drought_change
      temp_drought_change <- temp_drought_change[-1,]
    }
    # Check whether last date is a drought onset
    if (nrow(temp_drought_change)>0){
      if (tail(temp_drought_change,n=1)$reclass_USDM_change==1){
        # If last row is drought onset, add a new row with drought start date and no other statistics
        new_row <- data.frame("Phenocam" = phenocam, "Start_Date" = tail(temp_drought_change,n=1)$Date,
                              "Start_Year" = tail(temp_drought_change,n=1)$Year,
                              "End_Date" = NA, "End_Year" = NA, "Number_Weeks" = NA, "Sum_USDM" = NA)
        drought_statistics <- rbind(drought_statistics,new_row)
        
        # Remove last date from temp_drought_change
        temp_drought_change <- temp_drought_change[1:nrow(temp_drought_change)-1,]
      }
    }
    # Check for at least one complete drought period
    if (nrow(temp_drought_change)>=2){
      temp_start_dates <- filter(temp_drought_change, reclass_USDM_change == 1)
      temp_end_dates <- filter(temp_drought_change, reclass_USDM_change == -1)
      # Iterate through all drought periods
      for (k in 1:nrow(temp_start_dates)){
        # Calculate drought period statistics
        start_date <- temp_start_dates$Date[k]
        start_year <- temp_start_dates$Year[k]
        end_date <- temp_end_dates$Date[k]
        end_year <- temp_end_dates$Year[k]
        temp_drought_period <- filter(temp_df, Date >= start_date & Date < end_date)
        number_weeks <- nrow(temp_drought_period)
        sum_usdm <- sum(temp_drought_period$DM)
        
        new_row <- data.frame("Phenocam" = phenocam, "Start_Date" = start_date, "Start_Year" = start_year,
                              "End_Date" = end_date, "End_Year" = end_year,
                              "Number_Weeks" = number_weeks, "Sum_USDM" = sum_usdm)
        drought_statistics <- rbind(drought_statistics,new_row)
      }
    }
  }
}

drought_statistics$End_Date <- as.Date(drought_statistics$End_Date)

# Optional
rm(unique_phenocams, i, phenocam, temp_df, temp_drought_change, new_row,
   temp_start_dates, temp_end_dates, k, start_date, start_year, end_date, end_year,
   temp_drought_period, number_weeks, sum_usdm)

# USDM histogram
ggplot(data = filter(drought_statistics, Sum_USDM>0), aes(x=Number_Weeks)) + geom_histogram(binwidth=1)
ggplot(data = filter(drought_statistics, Sum_USDM>0), aes(x=Sum_USDM)) + geom_histogram(binwidth=1)

# Optional - save outputs
save(drought_statistics,phen_usdm_modified, file="outputs/drought_statistics.RData")

# # Sample plots
# p1 <- ggplot(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="XX")) +
#   geom_line(aes(x=Date, y=GCC_weekly), color="black") +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change<0)$Date, color="green") +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change>0)$Date, color="red") +
# p2 <- ggplot(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="XX")) +
#   geom_point(aes(x=Date, y=USDM_weekly))
# p1 / p2
# 
# p1 <- ggplot(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="GR")) +
#   geom_line(aes(x=Date, y=GCC_weekly), color="black") +
#   # drought -> more severe drought (USDM change > 0, USDM_weekly - full_USDM_change >= 1)
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change>0, USDM_weekly-full_USDM_change>=1)$Date, color="brown") +
#   # not drought -> drought (USDM change > 0, USDM_weekly >= 1, USDM_weekly - full_USDM_change < 1)
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change>0, USDM_weekly>=1, USDM_weekly-full_USDM_change<1)$Date, color="orange") +
#   # drought -> less severe drought (USDM change < 0, USDM_weekly > 0)
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change<0, USDM_weekly>=1)$Date, color="yellow") +
#   # drought -> not drought (USDM change < 0, USDM_weekly < 1, USDM_weekly - full_USDM_change > 0)
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change<0, USDM_weekly<1, USDM_weekly-full_USDM_change>=1)$Date, color="green")
# p2 <- ggplot(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="XX")) +
#   geom_point(aes(x=Date, y=USDM_weekly))
# p1 / p2
# 
# p3 <- ggplot() +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change<0, USDM_weekly>-1)$Date, color="green") +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="ibp", full_USDM_change>0, USDM_weekly>0)$Date, color="red") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="GR" & USDM_weekly>=1), aes(x=Date, y=GCC_weekly), color="brown") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="GR" & USDM_weekly<1), aes(x=Date, y=GCC_weekly), color="forestgreen") +
#   ggtitle("IBP GR weekly GCC")
# p4 <- ggplot() +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="GR" & USDM_weekly>=1), aes(x=Date, y=USDM_weekly), color="brown") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="ibp" & Veg_Type=="GR" & USDM_weekly<1), aes(x=Date, y=USDM_weekly), color="forestgreen")
# p3 / p4
# 
# p5 <- ggplot() +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="sevilletagrass", full_USDM_change<0, USDM_weekly>-1)$Date, color="green") +
#   geom_vline(xintercept = filter(vci_weekly_modified, Phenocam=="sevilletagrass", full_USDM_change>0, USDM_weekly>0)$Date, color="red") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="sevilletagrass" & Veg_Type=="GR" & USDM_weekly>=1), aes(x=Date, y=GCC_weekly), color="brown") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="sevilletagrass" & Veg_Type=="GR" & USDM_weekly<1), aes(x=Date, y=GCC_weekly), color="forestgreen") +
#   ggtitle("Sevilleta GR weekly GCC, vertical lines show increase (red) or decrease (green) in drought category")
# p6 <- ggplot() +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="sevilletagrass" & Veg_Type=="GR" & USDM_weekly>=1), aes(x=Date, y=USDM_weekly), color="brown") +
#   geom_point(data = filter(vci_weekly_modified, Phenocam=="sevilletagrass" & Veg_Type=="GR" & USDM_weekly<1), aes(x=Date, y=USDM_weekly), color="forestgreen") +
#   ggtitle("Sevilleta weekly USDM, color shows drought (dark red) or not drought (dark green)")
# p5 / p6