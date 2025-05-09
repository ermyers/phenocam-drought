# Code to find and subtract GCC baseline from PhenoCam data

# Note 04/29/2025: Redid yearly maxima to use highest yearly peak instead of
# maximum GCC.

# Note 04/18/2025: Baseline values not between two minima (i.e. before the first
# minimum or after the last minimum for a given site-year) are currently set to
# the value of the nearest minimum, as per na.approx(args, rule=2).

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(dplyr)
library(ggplot2)
library(zoo)
library(pracma)

# Load GCC data
load('outputs/phen_gcc.RData')
#load('outputs/growing_seasons_phen.RData')

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
  
  # Loop over each year, with 40 days padding at the beginning and end
  for (year in unique(phen_subset$Year)){
    # Find the start of the year, minus 40 days
    date_min <- head(filter(phen_subset, Year==year),n=1)$Date - 40
    # Find the end of the year, plus 40 days
    date_max <- tail(filter(phen_subset, Year==year),n=1)$Date + 40
    # Subset to the year with 40-day padding
    year_subset <- filter(phen_subset, Date>=date_min & Date<= date_max)
    # Find candidate peaks in the year subset
    candidate_peaks <- findpeaks(year_subset$smooth_gcc_90, peakpat = "[+]{1,}[0]*[-]{1,}")
    # Extract rows containing candidate peaks
    candidate_peaks_chrono <- year_subset[candidate_peaks[,2],]
    all_candidate_peaks_chrono <- candidate_peaks_chrono # unchanging record of candidate peaks
    # Sort candidate peaks in order of ascending value
    candidate_peaks_sorted <- dplyr::arrange(candidate_peaks_chrono,smooth_gcc_90)
    all_candidate_peaks_sorted <- candidate_peaks_sorted # unchanging record of candidate peaks
    
    # Check for a non-zero number of candidate peaks
    if(nrow(all_candidate_peaks_sorted)<1){
      # do nothing
    } else{
      # Loop over all candidate peaks in order of smallest to largest,
      # in order to eliminate smaller peaks near larger peaks
      for (j in 1:nrow(all_candidate_peaks_sorted)){
        date <- all_candidate_peaks_sorted$Date[j]
        peak_index <- which(candidate_peaks_chrono$Date==date)
        days_to_prev <- 365 # Placeholder value for number of days to previous chronological peak
        days_to_next <- 365 # Placeholder value for number of days to next chronological peak
        # If this is not the first peak in the year subset, find number of days to previous peak
        if(peak_index>1){
          days_to_prev <- min(365,date - candidate_peaks_chrono[peak_index-1,]$Date)
        }
        # If this is not the first peak in the year subset, find number of days to next peak
        if(peak_index<nrow(candidate_peaks_chrono)){
          days_to_next <- min(365,candidate_peaks_chrono[peak_index+1,]$Date - date)
        }
        # Eliminate candidate peaks with a larger peak within 30 days
        if(days_to_prev<30 | days_to_next<30){
          candidate_peaks_sorted <- candidate_peaks_sorted[candidate_peaks_sorted$Date != date,]
          candidate_peaks_chrono <- candidate_peaks_chrono[candidate_peaks_chrono$Date != date,]
        }
      }
      
      # Eliminate peaks occurring outside of the calendar year
      candidate_peaks_chrono <- filter(candidate_peaks_chrono, Year==year)
      candidate_peaks_sorted <- filter(candidate_peaks_sorted, Year==year)
    }
    
    # Select the greatest remaining peak as the yearly peak
    yearly_peak <- tail(candidate_peaks_sorted,1)
    
    # Record peak for year if:
    # 1) There is a yearly peak
    # 2) There are at least 200 days of non-interpolated data
    if(nrow(filter(phen_subset, Year==year & is.na(int_flag)==TRUE))>=200){
      yearly_peaks <- rbind(yearly_peaks,
                            data.frame(Phenocam = phenocam,
                                       Veg_Type = veg_type,
                                       ROI = roi,
                                       Year = year,
                                       Peak_GCC = yearly_peak$smooth_gcc_90,
                                       Peak_Date = yearly_peak$Date))
    }
    
    #year_subset <- filter(phen_subset, Year==year)
    #growing_seasons_subset <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year & Primary_Growing_Season=="yes")
    
    # if(nrow(growing_seasons_subset)==1){
    #   peak_gcc <- growing_seasons_subset$Peak_GCC[1]
    #   peak_date <- growing_seasons_subset$Peak[1]
    # } else{
    #   peak_gcc <- max(year_subset$smooth_gcc_90, na.rm=TRUE)
    #   peak_date <- max(filter(year_subset, smooth_gcc_90==peak_gcc)$Date)
    # }
    # #peak_gcc <- max(year_subset$smooth_gcc_90, na.rm=TRUE)
    # 
    # # only record peak for year with at least 200 days of non-interpolated data
    # if(nrow(filter(year_subset, is.na(int_flag)==TRUE))>=200){
    #   yearly_peaks <- rbind(yearly_peaks,
    #                         data.frame(Phenocam = phenocam,
    #                                    Veg_Type = veg_type,
    #                                    ROI = roi,
    #                                    Year = year,
    #                                    Peak_GCC = peak_gcc,
    #                                    Peak_Date = peak_date))
    # }
    
  }
}

# rm(growing_seasons_subset, phen_subset, year_subset, unique_phenocams, i, 
#    peak_date, peak_gcc, phenocam, roi, veg_type, year)

rm(all_candidate_peaks_chrono, all_candidate_peaks_sorted, candidate_peaks,
   candidate_peaks_chrono, candidate_peaks_sorted, phen_subset, unique_phenocams,
   year_subset, yearly_peak, date, date_max, date_min, days_to_next, days_to_prev,
   i, j, peak_index, phenocam, roi, veg_type, year)

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

