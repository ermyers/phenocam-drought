# Code to find and subtract GCC baseline from PhenoCam data

# Note 04/23/2025: Redid yearly maxima to use primary peak values (calculated
# from original phenometrics) when they exist, instead of maximum yearly GCC

# Note 04/18/2025: Baseline values not between two minima (i.e. before the first
# minimum or after the last minimum for a given site-year) are currently set to
# the value of the nearest minimum, as per na.approx(args, rule=2).

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(dplyr)
library(ggplot2)
library(zoo)

# Load GCC data
load('outputs/phen_gcc.RData')
load('outputs/growing_seasons_phen.RData')

##########################
# Calculate yearly maxima
##########################

# Create a list of unique phenocam ROIs
unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)

# Create a dataframe to hold information about yearly maxima
yearly_peaks <- data.frame()

# Find and record max in each site-year
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type,", roi = ",roi,sep=""))
  
  # Loop over each year
  for (year in unique(phen_subset$Year)){
    year_subset <- filter(phen_subset, Year==year)
    growing_seasons_subset <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year & Primary_Growing_Season=="yes")
    
    if(nrow(growing_seasons_subset)==1){
      peak_gcc <- growing_seasons_subset$Peak_GCC[1]
      peak_date <- growing_seasons_subset$Peak[1]
    } else{
      peak_gcc <- max(year_subset$smooth_gcc_90, na.rm=TRUE)
      peak_date <- max(filter(year_subset, smooth_gcc_90==peak_gcc)$Date)
    }
    #peak_gcc <- max(year_subset$smooth_gcc_90, na.rm=TRUE)
    
    # only record peak for year with at least 200 days of non-interpolated data
    if(nrow(filter(year_subset, is.na(int_flag)==TRUE))>=200){
      yearly_peaks <- rbind(yearly_peaks,
                            data.frame(Phenocam = phenocam,
                                       Veg_Type = veg_type,
                                       ROI = roi,
                                       Year = year,
                                       Peak_GCC = peak_gcc,
                                       Peak_Date = peak_date))
    }
    
  }
}

rm(growing_seasons_subset, phen_subset, year_subset, unique_phenocams, i, 
   peak_date, peak_gcc, phenocam, roi, veg_type, year)

#################################
# Calculate minima between peaks
#################################

# Create a list of unique phenocam ROIs within yearly peaks dataframe
unique_phenocams <- yearly_peaks %>% distinct(Phenocam,Veg_Type,ROI)

# Create a dataframe to hold information about minima
minima <- data.frame()

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- yearly_peaks %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type,", roi = ",roi,sep=""))
  
  # Look for minimum before first peak
  between_peak_subset <- filter(phen_gcc, Phenocam==phenocam & Veg_Type == veg_type & ROI == roi & Date<phen_subset$Peak_Date[1])
  min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
  minima <- rbind(minima,
                  data.frame(Phenocam = phenocam,
                             Veg_Type = veg_type,
                             ROI = roi,
                             Min_GCC = min_gcc,
                             Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
  
  # Look for minima between peaks
  if (nrow(phen_subset)>1){
    for (j in 2:nrow(phen_subset)){
      between_peak_subset <- filter(phen_gcc, Phenocam==phenocam & Veg_Type == veg_type & ROI == roi & Date>phen_subset$Peak_Date[j-1] & Date<phen_subset$Peak_Date[j])
      min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
      min_date <- filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date
      minima <- rbind(minima,
                      data.frame(Phenocam = phenocam,
                                 Veg_Type = veg_type,
                                 ROI = roi,
                                 Min_GCC = min_gcc,
                                 Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
    }
  }
  
  # Look for minimum after last peak
  between_peak_subset <- filter(phen_gcc, Phenocam==phenocam & Veg_Type == veg_type & ROI == roi & Date > phen_subset$Peak_Date[nrow(phen_subset)])
  min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
  minima <- rbind(minima,
                  data.frame(Phenocam = phenocam,
                             Veg_Type = veg_type,
                             ROI = roi,
                             Min_GCC = min_gcc,
                             Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
}

rm(between_peak_subset, phen_subset, unique_phenocams, i, j, min_date, min_gcc, phenocam, roi, veg_type)

###################################################
# Interpolate between minima to get a GCC baseline
###################################################

# Create a list of unique phenocam ROIs
unique_phenocams <- yearly_peaks %>% distinct(Phenocam,Veg_Type,ROI)

# Create a dataframe to hold phenocam GCC and GCC baseline
phen_gcc_with_baseline <- data.frame()

# Add linearly interpolated GCC baseline to phenocam GCC
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type,", roi = ",roi,sep=""))
  
  phen_subset <- left_join(phen_subset, minima, by=join_by(Phenocam, Veg_Type, ROI, Date==Min_Date))
  phen_subset <- mutate(phen_subset, Baseline_GCC = na.approx(Min_GCC, rule=2))
  phen_gcc_with_baseline <- rbind(phen_gcc_with_baseline, phen_subset)
}

rm(phen_subset, unique_phenocams, i, phenocam, roi, veg_type)

#############################
# Subtract baseline from GCC
#############################

phen_gcc_with_baseline <- mutate(phen_gcc_with_baseline, smooth_gcc_baseline_subtracted = smooth_gcc_90-Baseline_GCC)

###############
# Plot results
###############

ggplot() +
  geom_point(data=filter(phen_gcc_with_baseline, Phenocam=="ibp" & Veg_Type=="XX" & ROI=="1000"), aes(x=Date,y=smooth_gcc_90)) +
  geom_line(data=filter(phen_gcc_with_baseline, Phenocam=="ibp" & Veg_Type=="XX" & ROI=="1000"), aes(x=Date,y=Baseline_GCC), color='red') +
  ggtitle("ibp")

ggplot() +
  geom_point(data=filter(phen_gcc_with_baseline, Phenocam=="ibp" & Veg_Type=="XX" & ROI=="1000"), aes(x=Date,y=smooth_gcc_baseline_subtracted)) +
  ggtitle("ibp")

ggplot() +
  geom_point(data=filter(phen_gcc_with_baseline, Phenocam=="vaira" & Veg_Type=="GR" & ROI=="1000"), aes(x=Date,y=smooth_gcc_90)) +
  geom_line(data=filter(phen_gcc_with_baseline, Phenocam=="vaira" & Veg_Type=="GR" & ROI=="1000"), aes(x=Date,y=Baseline_GCC), color='red') +
  ggtitle("vaira")

ggplot() +
  geom_point(data=filter(phen_gcc_with_baseline, Phenocam=="vaira" & Veg_Type=="GR" & ROI=="1000"), aes(x=Date,y=smooth_gcc_baseline_subtracted)) +
  ggtitle("vaira")

###############
# Save results
###############

save(phen_gcc_with_baseline, file="outputs/phen_gcc_with_baseline.RData")

