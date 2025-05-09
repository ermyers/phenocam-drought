# Calculate phenometrics using baseline-subtracted GCC

# Note 04/18/2025: Mostly copied code over from calculate_phenometrics.R,
# changing variable names as appropriate. Got code to run, but haven't checked
# extensively to see if it's doing what I want.

# Note 3/19/25: Added code under count number of site-years section to determine
# "primary" growing season in case of multiple growing seasons per year.
# Have not yet decided whether to calculate average phenometrics from "primary"
# growing season or from earliest start and latest end.

# Note 10/29/24: Code to eliminate smaller peaks with larger peaks within 30 days
# throws a warning for peaks at very beginning and end of time series (but should
# still return good results - these peaks will be eliminated at the check for 35% variation)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(dplyr)
library(pracma)
library(ggplot2)
library(lubridate)

# Load GCC data
load('outputs/phen_gcc_with_baseline.RData')

###############################################
# Calculate growing season dates and statistics
###############################################

# Create dataframe to hold growing season information

# Peak GCC: GCC value at peak
# Peak_GCC_baseline_subtracted: Baseline-subtracted GCC value at peak
# Amplitude_GCC_yearly: Yearly GCC amplitude (max - min of yearly baseline-subtracted GCC)
# Amplitude_GCC_peak: Peak GCC amplitude (peak - the lower of pre-/post-peak min, from baseline-subtracted GCC)
# Cumulative_GCC_yearly: Sum of baseline-subtracted GCC over the whole year
# Cumulative_GCC_to_EOS: Sum of baseline-subtracted GCC from start of year to EOS_25

growing_seasons_phen <- data.frame(matrix(ncol = 24, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Year",
       "Peak_GCC","Peak_GCC_baseline_subtracted",
       "Amplitude_GCC_yearly","Amplitude_GCC_peak",
       "Cumulative_GCC_yearly","Cumulative_GCC_to_EOS",
       "SOS_15","SOS_25","SOS_50","Peak","EOS_50","EOS_25","EOS_15",
       "SOS_15_flag","SOS_25_flag","SOS_50_flag","Peak_flag","EOS_50_flag","EOS_25_flag","EOS_15_flag")
colnames(growing_seasons_phen) <- x

# Loop over each phenocam

# Create a list of unique phenocam ROIs
unique_phenocams <- phen_gcc_with_baseline %>% distinct(Phenocam,Veg_Type,ROI)
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc_with_baseline %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Loop over each year, plus or minus half-year buffer
  for (year in unique(phen_subset$Year)){
    # Find the start of the year, minus 185 days
    date_min <- head(filter(phen_subset, Year==year),n=1)$Date - 185
    # Find the end of the year, plus 185 days
    date_max <- tail(filter(phen_subset, Year==year),n=1)$Date + 185
    # Subset to the year with 185-day padding
    date_subset <- filter(phen_subset, Date>=date_min & Date<= date_max)
    # Find the minimum value of GCC within this date subset (should be near 0)
    min_GCC <- min(date_subset$smooth_gcc_baseline_subtracted,na.rm=TRUE)
    # Find the maximum value of GCC within this date subset
    max_GCC <- max(date_subset$smooth_gcc_baseline_subtracted,na.rm=TRUE)
    # Find candidate peaks in the date subset
    candidate_peaks <- findpeaks(date_subset$smooth_gcc_baseline_subtracted, peakpat = "[+]{1,}[0]*[-]{1,}")
    # Extract rows containing candidate peaks
    candidate_peaks_chrono <- date_subset[candidate_peaks[,2],]
    all_candidate_peaks_chrono <- candidate_peaks_chrono # unchanging record of candidate peaks
    # Sort candidate peaks in order of ascending value
    candidate_peaks_sorted <- dplyr::arrange(candidate_peaks_chrono,smooth_gcc_baseline_subtracted) # candidate peaks sorted in ascending order of value
    all_candidate_peaks_sorted <- candidate_peaks_sorted # unchanging record of candidate peaks
    
    # Check for a nonzero number of candidate peaks
    if(nrow(all_candidate_peaks_sorted)<1){
      # do nothing
    } else{
      
      # Loop over all candidate peaks in order of smallest to largest,
      # eliminating peaks that don't meet the criteria for inclusion
      for (j in 1:nrow(all_candidate_peaks_sorted)){
        date <- all_candidate_peaks_sorted$Date[j]
        peak_index <- which(candidate_peaks_chrono$Date==date)
        
        # Define the search window before and after the peak
        if(peak_index>1){
          dist_to_prev <- min(365,date - candidate_peaks_chrono[peak_index-1,]$Date)
        } else{
          dist_to_prev <- 365
        }
        if(peak_index<nrow(candidate_peaks_chrono)){
          dist_to_next <- min(365,candidate_peaks_chrono[peak_index+1,]$Date - date)
        } else{
          dist_to_next <- 365
        }
        
        # Eliminate candidate peaks with larger peaks within 30 days
        if(dist_to_prev<30 | dist_to_next<30){
          candidate_peaks_sorted <- candidate_peaks_sorted[candidate_peaks_sorted$Date != date,]
          candidate_peaks_chrono <- candidate_peaks_chrono[candidate_peaks_chrono$Date != date,]
        } else{
          # Subset to the time period around this candidate peak
          peak_subset <- filter(phen_subset, Date>=date-dist_to_prev & Date<=date+dist_to_next)
          # Calculate the pre-peak minimum value and date
          pre_peak_min <- min(filter(peak_subset,Date<=date-30)$smooth_gcc_baseline_subtracted,na.rm=TRUE)
          pre_peak_min_date <- min(filter(peak_subset,smooth_gcc_baseline_subtracted==pre_peak_min)$Date,na.rm=TRUE)
          # Calculate the post-peak minimum value and date
          post_peak_min <- min(filter(peak_subset, Date>=date+30)$smooth_gcc_baseline_subtracted,na.rm=TRUE)
          post_peak_min_date <- max(filter(peak_subset,smooth_gcc_baseline_subtracted==post_peak_min)$Date,na.rm=TRUE)
          # Calculate the peak value (baseline subtracted GCC) and peak GCC value (non-baseline-subtracted GCC)
          peak_value <- filter(peak_subset, Date==date)$smooth_gcc_baseline_subtracted
          peak_gcc_value <- filter(peak_subset, Date==date)$smooth_gcc_90
          
          # Eliminate peaks with difference between minimum and maximum <35% of total variation within date subset
          if((peak_value-pre_peak_min)/(max_GCC-min_GCC)<0.35 | (peak_value-post_peak_min)/(max_GCC-min_GCC)<0.35){
            candidate_peaks_sorted <- candidate_peaks_sorted[candidate_peaks_sorted$Date != date,]
            candidate_peaks_chrono <- candidate_peaks_chrono[candidate_peaks_chrono$Date != date,]
          } else{
            
            # Calculate peak statistics for any remaining candidate peaks
            # that occur within the calendar year
            if(peak_subset[peak_subset$Date==date,]$Year == year){
              
              # Calculate peak timing and season start and end dates, based on % thresholds
              peak_amplitude <- peak_value-min(c(pre_peak_min,post_peak_min))
              SOS_15_thresh <- 0.15*(peak_value-pre_peak_min) + pre_peak_min
              SOS_25_thresh <- 0.25*(peak_value-pre_peak_min) + pre_peak_min
              SOS_50_thresh <- 0.50*(peak_value-pre_peak_min) + pre_peak_min
              EOS_50_thresh <- 0.50*(peak_value-post_peak_min) + post_peak_min
              EOS_25_thresh <- 0.25*(peak_value-post_peak_min) + post_peak_min
              EOS_15_thresh <- 0.15*(peak_value-post_peak_min) + post_peak_min
              SOS_15 <- min(filter(peak_subset, smooth_gcc_baseline_subtracted>SOS_15_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              SOS_25 <- min(filter(peak_subset, smooth_gcc_baseline_subtracted>SOS_25_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              SOS_50 <- min(filter(peak_subset, smooth_gcc_baseline_subtracted>SOS_50_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              EOS_50 <- max(filter(peak_subset, smooth_gcc_baseline_subtracted>EOS_50_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              EOS_25 <- max(filter(peak_subset, smooth_gcc_baseline_subtracted>EOS_25_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              EOS_15 <- max(filter(peak_subset, smooth_gcc_baseline_subtracted>EOS_15_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              date_peak <- peak_subset[peak_subset$Date==date,]$Date
              
              # Calculate cumulative GCC (start of year to EOS 25) and cumulative yearly GCC (start to end of year)
              year_subset <- filter(phen_subset, Year==year)
              if(min(year_subset$DOY, na.rm=TRUE)!=1){
                cumulative_gcc <- NA
                cumulative_gcc_yearly <- NA
              } else{
                cumulative_gcc <- sum(filter(year_subset, Date <= EOS_25)$smooth_gcc_baseline_subtracted, na.rm=TRUE)
                if(max(year_subset$DOY, na.rm=TRUE)<365){
                  cumulative_gcc_yearly <- NA
                } else{
                  cumulative_gcc_yearly <- sum(year_subset$smooth_gcc_baseline_subtracted, na.rm=TRUE)
                }
              }
              
              # Calculate yearly GCC amplitude
              yearly_amplitude <- max(year_subset$smooth_gcc_baseline_subtracted) - min(year_subset$smooth_gcc_baseline_subtracted)
              
              # Check for int_flags within two weeks of key points
              SOS_15_flag <- sum(filter(peak_subset, Date >= SOS_15-7 & Date <= SOS_15+7)$int_flag, na.rm=TRUE)
              SOS_25_flag <- sum(filter(peak_subset, Date >= SOS_25-7 & Date <= SOS_25+7)$int_flag, na.rm=TRUE)
              SOS_50_flag <- sum(filter(peak_subset, Date >= SOS_50-7 & Date <= SOS_50+7)$int_flag, na.rm=TRUE)
              EOS_50_flag <- sum(filter(peak_subset, Date >= EOS_50-7 & Date <= EOS_50+7)$int_flag, na.rm=TRUE)
              EOS_25_flag <- sum(filter(peak_subset, Date >= EOS_25-7 & Date <= EOS_25+7)$int_flag, na.rm=TRUE)
              EOS_15_flag <- sum(filter(peak_subset, Date >= EOS_15-7 & Date <= EOS_15+7)$int_flag, na.rm=TRUE)
              peak_flag <- sum(filter(peak_subset, Date >= date_peak-7 & Date <= date_peak+7)$int_flag, na.rm=TRUE)
              
              # Record peak values
              new_row <- data.frame(Phenocam = phenocam,
                                    Veg_Type = veg_type,
                                    ROI = roi,
                                    Year = year,
                                    Peak_GCC = peak_gcc_value,
                                    Peak_GCC_baseline_subtracted = peak_value,
                                    Amplitude_GCC_yearly = yearly_amplitude,
                                    Amplitude_GCC_peak = peak_amplitude,
                                    Cumulative_GCC_yearly = cumulative_gcc_yearly,
                                    Cumulative_GCC_to_EOS = cumulative_gcc,
                                    SOS_15 = SOS_15,
                                    SOS_25 = SOS_25,
                                    SOS_50 = SOS_50,
                                    Peak = date_peak,
                                    EOS_50 = EOS_50,
                                    EOS_25 = EOS_25,
                                    EOS_15 = EOS_15,
                                    SOS_15_flag = SOS_15_flag,
                                    SOS_25_flag = SOS_25_flag,
                                    SOS_50_flag = SOS_50_flag,
                                    Peak_flag = peak_flag,
                                    EOS_50_flag = EOS_50_flag,
                                    EOS_25_flag = EOS_25_flag,
                                    EOS_15_flag = EOS_15_flag)
              growing_seasons_phen <- rbind(growing_seasons_phen,new_row)
              
            }
          }
        }
      }
    }
  }
}

# Add a new field (Primary_Growing_Season) to growing_seasons_phen
# Primary growing seasons will be designated in the next section of code
growing_seasons_phen <- cbind(growing_seasons_phen, "Primary_Growing_Season" = "no")

# Optional - remove variables
rm(all_candidate_peaks_chrono,all_candidate_peaks_sorted,candidate_peaks,candidate_peaks_chrono,
   candidate_peaks_sorted,date_subset,new_row,peak_subset,phen_subset,unique_phenocams,year_subset)
rm(cumulative_gcc,cumulative_gcc_yearly,date,date_max,date_min,date_peak,dist_to_next,dist_to_prev,EOS_15,EOS_15_flag,
   EOS_15_thresh,EOS_25,EOS_25_flag,EOS_25_thresh,EOS_50,EOS_50_flag,EOS_50_thresh,
   i,j,max_GCC,min_GCC,peak_amplitude,peak_flag,peak_gcc_value,peak_index,peak_value,phenocam,post_peak_min,
   post_peak_min_date,pre_peak_min,pre_peak_min_date,roi,SOS_15,SOS_15_flag,SOS_15_thresh,
   SOS_25,SOS_25_flag,SOS_25_thresh,SOS_50,SOS_50_flag,SOS_50_thresh,veg_type,x,year,yearly_amplitude)

# Optional - plot selected results
ggplot() +
  geom_point(data=filter(phen_gcc_with_baseline, Phenocam=="ibp" & Veg_Type=="XX"), aes(x=Date,y=smooth_gcc_baseline_subtracted)) +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp" & Veg_Type=="XX")$SOS_15, color="yellow") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp"& Veg_Type=="XX")$EOS_15, color="brown") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp" & Veg_Type=="XX")$SOS_50, color="green") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp"& Veg_Type=="XX")$EOS_50, color="red")

############################################################################
# Count number of growing seasons per site-year, including 0 or more than 1
############################################################################

# Create dataframe to hold site-year information

# Peak_GCC_yearly: Max GCC over entire site-year
# Peak_GCC_baseline_subtracted_yearly: Max baseline-subtracted GCC over entire site-year
# Amplitude_GCC_yearly: Yearly GCC amplitude (max - min of yearly baseline-subtracted GCC)
# Cumulative_GCC_yearly: Sum of baseline-subtracted GCC over the whole year
# Number_Growing_Seasons: Number of valid growing seasons
# Cols. 10-15: Growing season statistics for the primary growing season for a given site-year.
#              Given as NA when there is no growing season.

site_year_statistics <- data.frame(matrix(ncol = 15, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Year",
       "Peak_GCC_yearly","Peak_GCC_baseline_subtracted_yearly",
       "Amplitude_GCC_yearly", "Cumulative_GCC_yearly",
       "Number_Growing_Seasons","SOS_25","Peak","EOS_25","GSL_25",
       "Peak_GCC","Peak_GCC_baseline_subtracted")
colnames(site_year_statistics) <- x

# Loop over each phenocam

# Create a list of unique phenocam ROIs
unique_phenocams <- phen_gcc_with_baseline %>% distinct(Phenocam,Veg_Type,ROI)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc_with_baseline %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Loop over each year
  for (year in unique(phen_subset$Year)){
    # check that we have non-interpolated data for most of the year
    if(nrow(filter(phen_subset,Year==year & is.na(int_flag)==TRUE))>=200){
      # Calculate yearly statistics
      year_gcc_subset <- filter(phen_subset, Year==year)
      peak_gcc_yearly <- max(year_gcc_subset$smooth_gcc_90, na.rm=TRUE)
      peak_gcc_baseline_subtracted_yearly <- max(year_gcc_subset$smooth_gcc_baseline_subtracted, na.rm=TRUE)
      amplitude_gcc_yearly <- peak_gcc_baseline_subtracted_yearly - min(year_gcc_subset$smooth_gcc_baseline_subtracted, na.rm=TRUE)
      if(nrow(year_gcc_subset)>=365){
        cumulative_gcc_yearly <- sum(year_gcc_subset$smooth_gcc_baseline_subtracted, na.rm=TRUE)
      } else{
        cumulative_gcc_yearly <- NA
      }
      # Calculate growing season statistics
      growing_season_subset <- filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)
      # Count number of growing seasons
      number_growing_seasons <- nrow(growing_season_subset)
      # If no valid growing seasons, set growing season statistics to NA
      if(number_growing_seasons==0){
        sos_25 <- NA
        peak <- NA
        eos_25 <- NA
        gsl_25 <- NA
        peak_gcc <- NA
        peak_gcc_baseline_subtracted <- NA
      } else{
        # If at least one valid growing season, record growing season statistics
        # for growing season with the greatest peak
        growing_season_subset <- growing_season_subset %>% slice(which.max(Peak_GCC_baseline_subtracted))
        sos_25 <- growing_season_subset$SOS_25
        peak <- growing_season_subset$Peak
        eos_25 <- growing_season_subset$EOS_25
        gsl_25 <- as.numeric(eos_25-sos_25)
        peak_gcc <- growing_season_subset$Peak_GCC
        peak_gcc_baseline_subtracted <- growing_season_subset$Peak_GCC_baseline_subtracted
      }
      
      new_row <- data.frame(Phenocam = phenocam,
                            Veg_Type = veg_type,
                            ROI = roi,
                            Year = year,
                            Peak_GCC_yearly = peak_gcc_yearly,
                            Peak_GCC_baseline_subtracted_yearly = peak_gcc_baseline_subtracted_yearly,
                            Amplitude_GCC_yearly = amplitude_gcc_yearly,
                            Cumulative_GCC_yearly = cumulative_gcc_yearly,
                            Number_Growing_Seasons = number_growing_seasons,
                            SOS_25 = sos_25,
                            Peak = peak,
                            EOS_25 = eos_25,
                            GSL_25 = gsl_25,
                            Peak_GCC = peak_gcc,
                            Peak_GCC_baseline_subtracted = peak_gcc_baseline_subtracted)
      site_year_statistics <- rbind(site_year_statistics,new_row)
      
      # Denote primary growing season in growing_seasons_phen
      # Primary growing season is currently determined by the growing season with the greatest peak GCC
      if(number_growing_seasons>0){
        growing_seasons_phen[growing_seasons_phen$Phenocam==phenocam & growing_seasons_phen$Veg_Type==veg_type & growing_seasons_phen$ROI==roi & growing_seasons_phen$Year==year & growing_seasons_phen$Peak_GCC==peak_gcc,]$Primary_Growing_Season <- "yes"
      }
    }
  }
}

# Optional - remove variables
rm(growing_season_subset, new_row, phen_subset, unique_phenocams, year_gcc_subset,
   amplitude_gcc_yearly, cumulative_gcc_yearly, eos_25, gsl_25, i, number_growing_seasons,
   peak, peak_gcc, peak_gcc_baseline_subtracted, peak_gcc_baseline_subtracted_yearly,
   peak_gcc_yearly, phenocam, roi, sos_25, veg_type, x, year)

###############################################################
# Identify typical growing season values for each phenocam ROI
###############################################################

growing_seasons_phen_ave <- data.frame(matrix(ncol = 25, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Peak_GCC_mean","Peak_GCC_baseline_subtracted_mean",
       "Amplitude_GCC_yearly_mean","Amplitude_GCC_peak_mean",
       "Cumulative_GCC_yearly_mean","Cumulative_GCC_mean","SOS_15_mean","SOS_25_mean",
       "SOS_50_mean","Peak_mean","EOS_50_mean","EOS_25_mean","EOS_15_mean",
       "Peak_GCC_std","SOS_15_std","SOS_25_std","SOS_50_std","Peak_std",
       "EOS_50_std","EOS_25_std","EOS_15_std","Number_Years")
colnames(growing_seasons_phen_ave) <- x

unique_phenocams <- phen_gcc_with_baseline %>% distinct(Phenocam,Veg_Type,ROI)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- growing_seasons_phen %>%
    filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi) %>%
    mutate(SOS_15 = as.numeric(SOS_15 - as.Date(paste(Year-1,"-12-31",sep=""))),
           SOS_25 = as.numeric(SOS_25 - as.Date(paste(Year-1,"-12-31",sep=""))),
           SOS_50 = as.numeric(SOS_50 - as.Date(paste(Year-1,"-12-31",sep=""))),
           Peak = yday(Peak),
           EOS_50 = as.numeric(EOS_50 - as.Date(paste(Year-1,"-12-31",sep=""))),
           EOS_25 = as.numeric(EOS_25 - as.Date(paste(Year-1,"-12-31",sep=""))),
           EOS_15 = as.numeric(EOS_15 - as.Date(paste(Year-1,"-12-31",sep=""))))
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Combine multiple growing seasons if necessary
  for(year in unique(phen_subset$Year)){
    if (nrow(filter(phen_subset, Year==year)) > 1){
      peak_gcc <- max(filter(phen_subset, Year==year)$Peak_GCC, na.rm = TRUE)
      peak_gcc_baseline_subtracted <- max(filter(phen_subset, Year==year)$Peak_GCC_baseline_subtracted, na.rm = TRUE)
      yearly_amplitude <- max(filter(phen_subset, Year==year)$Amplitude_GCC_yearly, na.rm = TRUE)
      peak_amplitude <- max(filter(phen_subset, Year==year)$Amplitude_GCC_peak, na.rm = TRUE)
      cumulative_gcc_yearly <- max(filter(phen_subset, Year==year)$Cumulative_GCC_yearly, na.rm = TRUE)
      cumulative_gcc <- max(filter(phen_subset, Year==year)$Cumulative_GCC_to_EOS, na.rm = TRUE)
      if(cumulative_gcc==-Inf){
        cumulative_gcc <- NA
      }
      sos_15 <- min(filter(phen_subset, Year==year)$SOS_15, na.rm = TRUE)
      sos_25 <- min(filter(phen_subset, Year==year)$SOS_25, na.rm = TRUE)
      sos_50 <- min(filter(phen_subset, Year==year)$SOS_50, na.rm = TRUE)
      peak <- filter(phen_subset, Year==year & Peak_GCC==peak_gcc)$Peak
      eos_50 <- max(filter(phen_subset, Year==year)$EOS_50, na.rm = TRUE)
      eos_25 <- max(filter(phen_subset, Year==year)$EOS_25, na.rm = TRUE)
      eos_15 <- max(filter(phen_subset, Year==year)$EOS_15, na.rm = TRUE)
      
      new_row <- data.frame(Phenocam = phenocam,
                            Veg_Type = veg_type,
                            ROI = roi,
                            Year = year,
                            Peak_GCC = peak_gcc,
                            Peak_GCC_baseline_subtracted = peak_gcc_baseline_subtracted,
                            Amplitude_GCC_yearly = yearly_amplitude,
                            Amplitude_GCC_peak = peak_amplitude,
                            Cumulative_GCC_yearly = cumulative_gcc_yearly,
                            Cumulative_GCC_to_EOS = cumulative_gcc,
                            SOS_15 = sos_15,
                            SOS_25 = sos_25,
                            SOS_50 = sos_50,
                            Peak = peak,
                            EOS_50 = eos_50,
                            EOS_25 = eos_25,
                            EOS_15 = eos_15,
                            SOS_15_flag = filter(phen_subset, Year==year & SOS_15==sos_15)$SOS_15_flag,
                            SOS_25_flag = filter(phen_subset, Year==year & SOS_25==sos_25)$SOS_25_flag,
                            SOS_50_flag = filter(phen_subset, Year==year & SOS_50==sos_50)$SOS_50_flag,
                            Peak_flag = filter(phen_subset, Year==year & Peak==peak)$Peak_flag,
                            EOS_50_flag = filter(phen_subset, Year==year & EOS_50==eos_50)$EOS_50_flag,
                            EOS_25_flag = filter(phen_subset, Year==year & EOS_25==eos_25)$EOS_25_flag,
                            EOS_15_flag = filter(phen_subset, Year==year & EOS_15==eos_15)$EOS_15_flag,
                            Primary_Growing_Season = "yes")
      phen_subset <- filter(phen_subset, Year != year)
      phen_subset <- rbind(phen_subset, new_row)
      rm(new_row)
    }
  }
  
  # Calculate average values and add to typical growing seasons
  new_row <- data.frame(Phenocam = phenocam,
                        Veg_Type = veg_type,
                        ROI = roi,
                        Peak_GCC_mean = mean(phen_subset$Peak_GCC, na.rm=TRUE),
                        Peak_GCC_baseline_subtracted_mean = mean(phen_subset$Peak_GCC_baseline_subtracted, na.rm=TRUE),
                        Amplitude_GCC_yearly_mean = mean(phen_subset$Amplitude_GCC_yearly, na.rm=TRUE),
                        Amplitude_GCC_peak_mean = mean(phen_subset$Amplitude_GCC_peak, na.rm=TRUE),
                        Cumulative_GCC_yearly_mean = mean(phen_subset$Cumulative_GCC_yearly, na.rm=TRUE),
                        Cumulative_GCC_to_EOS_mean = mean(phen_subset$Cumulative_GCC_to_EOS, na.rm=TRUE),
                        SOS_15_mean = mean(phen_subset$SOS_15, na.rm=TRUE),
                        SOS_25_mean = mean(phen_subset$SOS_25, na.rm=TRUE),
                        SOS_50_mean = mean(phen_subset$SOS_50, na.rm=TRUE),
                        Peak_mean = mean(phen_subset$Peak, na.rm=TRUE),
                        EOS_50_mean = mean(phen_subset$EOS_50, na.rm=TRUE),
                        EOS_25_mean = mean(phen_subset$EOS_25, na.rm=TRUE),
                        EOS_15_mean = mean(phen_subset$EOS_15, na.rm=TRUE),
                        Peak_GCC_std = sd(phen_subset$Peak_GCC, na.rm=TRUE),
                        SOS_15_std = sd(phen_subset$SOS_15, na.rm=TRUE),
                        SOS_25_std = sd(phen_subset$SOS_25, na.rm=TRUE),
                        SOS_50_std = sd(phen_subset$SOS_50, na.rm=TRUE),
                        Peak_std = sd(phen_subset$Peak, na.rm=TRUE),
                        EOS_50_std = sd(phen_subset$EOS_50, na.rm=TRUE),
                        EOS_25_std = sd(phen_subset$EOS_25, na.rm=TRUE),
                        EOS_15_std = sd(phen_subset$EOS_15, na.rm=TRUE),
                        Number_Years = nrow(phen_subset))
  growing_seasons_phen_ave <- rbind(growing_seasons_phen_ave,new_row)
}

# Optional - remove variables
rm(cumulative_gcc_yearly, cumulative_gcc, new_row, phen_subset, unique_phenocams, eos_15, eos_25, eos_50, 
   i, peak, peak_amplitude, peak_gcc, peak_gcc_baseline_subtracted, phenocam, roi, sos_15, sos_25, sos_50, veg_type, x, year, yearly_amplitude)

# Optional - save results
save(growing_seasons_phen, site_year_statistics, growing_seasons_phen_ave, file = "outputs/growing_seasons_phen_baseline_subtracted.RData")

##########################
# Limit data to 2016-2024
##########################

# Season start and end
growing_seasons_phen <- filter(growing_seasons_phen, Year >= 2016)

# No or multiple growing seasons
site_year_statistics <- filter(site_year_statistics, Year >= 2016)

# Typical growing season values (2016-2024 only)
growing_seasons_phen_ave <- data.frame(matrix(ncol = 25, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Peak_GCC_mean","Peak_GCC_baseline_subtracted_mean",
       "Amplitude_GCC_yearly_mean","Amplitude_GCC_peak_mean",
       "Cumulative_GCC_yearly_mean","Cumulative_GCC_mean","SOS_15_mean","SOS_25_mean",
       "SOS_50_mean","Peak_mean","EOS_50_mean","EOS_25_mean","EOS_15_mean",
       "Peak_GCC_std","SOS_15_std","SOS_25_std","SOS_50_std","Peak_std",
       "EOS_50_std","EOS_25_std","EOS_15_std","Number_Years")
colnames(growing_seasons_phen_ave) <- x

unique_phenocams <- growing_seasons_phen %>% distinct(Phenocam,Veg_Type,ROI)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- growing_seasons_phen %>%
    filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi) %>%
    mutate(SOS_15 = as.numeric(SOS_15 - as.Date(paste(Year-1,"-12-31",sep=""))),
           SOS_25 = as.numeric(SOS_25 - as.Date(paste(Year-1,"-12-31",sep=""))),
           SOS_50 = as.numeric(SOS_50 - as.Date(paste(Year-1,"-12-31",sep=""))),
           Peak = yday(Peak),
           EOS_50 = as.numeric(EOS_50 - as.Date(paste(Year-1,"-12-31",sep=""))),
           EOS_25 = as.numeric(EOS_25 - as.Date(paste(Year-1,"-12-31",sep=""))),
           EOS_15 = as.numeric(EOS_15 - as.Date(paste(Year-1,"-12-31",sep=""))))
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Combine multiple growing seasons if necessary
  for(year in unique(phen_subset$Year)){
    if (nrow(filter(phen_subset, Year==year)) > 1){
      peak_gcc <- max(filter(phen_subset, Year==year)$Peak_GCC, na.rm = TRUE)
      peak_gcc_baseline_subtracted <- max(filter(phen_subset, Year==year)$Peak_GCC_baseline_subtracted, na.rm = TRUE)
      yearly_amplitude <- max(filter(phen_subset, Year==year)$Amplitude_GCC_yearly, na.rm = TRUE)
      peak_amplitude <- max(filter(phen_subset, Year==year)$Amplitude_GCC_peak, na.rm = TRUE)
      cumulative_gcc_yearly <- max(filter(phen_subset, Year==year)$Cumulative_GCC_yearly, na.rm = TRUE)
      cumulative_gcc <- max(filter(phen_subset, Year==year)$Cumulative_GCC_to_EOS, na.rm = TRUE)
      if(cumulative_gcc==-Inf){
        cumulative_gcc <- NA
      }
      sos_15 <- min(filter(phen_subset, Year==year)$SOS_15, na.rm = TRUE)
      sos_25 <- min(filter(phen_subset, Year==year)$SOS_25, na.rm = TRUE)
      sos_50 <- min(filter(phen_subset, Year==year)$SOS_50, na.rm = TRUE)
      peak <- filter(phen_subset, Year==year & Peak_GCC==peak_gcc)$Peak
      eos_50 <- max(filter(phen_subset, Year==year)$EOS_50, na.rm = TRUE)
      eos_25 <- max(filter(phen_subset, Year==year)$EOS_25, na.rm = TRUE)
      eos_15 <- max(filter(phen_subset, Year==year)$EOS_15, na.rm = TRUE)
      
      new_row <- data.frame(Phenocam = phenocam,
                            Veg_Type = veg_type,
                            ROI = roi,
                            Year = year,
                            Peak_GCC = peak_gcc,
                            Peak_GCC_baseline_subtracted = peak_gcc_baseline_subtracted,
                            Amplitude_GCC_yearly = yearly_amplitude,
                            Amplitude_GCC_peak = peak_amplitude,
                            Cumulative_GCC_yearly = cumulative_gcc_yearly,
                            Cumulative_GCC_to_EOS = cumulative_gcc,
                            SOS_15 = sos_15,
                            SOS_25 = sos_25,
                            SOS_50 = sos_50,
                            Peak = peak,
                            EOS_50 = eos_50,
                            EOS_25 = eos_25,
                            EOS_15 = eos_15,
                            SOS_15_flag = filter(phen_subset, Year==year & SOS_15==sos_15)$SOS_15_flag,
                            SOS_25_flag = filter(phen_subset, Year==year & SOS_25==sos_25)$SOS_25_flag,
                            SOS_50_flag = filter(phen_subset, Year==year & SOS_50==sos_50)$SOS_50_flag,
                            Peak_flag = filter(phen_subset, Year==year & Peak==peak)$Peak_flag,
                            EOS_50_flag = filter(phen_subset, Year==year & EOS_50==eos_50)$EOS_50_flag,
                            EOS_25_flag = filter(phen_subset, Year==year & EOS_25==eos_25)$EOS_25_flag,
                            EOS_15_flag = filter(phen_subset, Year==year & EOS_15==eos_15)$EOS_15_flag,
                            Primary_Growing_Season = "yes")
      phen_subset <- filter(phen_subset, Year != year)
      phen_subset <- rbind(phen_subset, new_row)
      rm(new_row)
    }
  }
  
  # Calculate average values and add to typical growing seasons
  new_row <- data.frame(Phenocam = phenocam,
                        Veg_Type = veg_type,
                        ROI = roi,
                        Peak_GCC_mean = mean(phen_subset$Peak_GCC, na.rm=TRUE),
                        Peak_GCC_baseline_subtracted_mean = mean(phen_subset$Peak_GCC_baseline_subtracted, na.rm=TRUE),
                        Amplitude_GCC_yearly_mean = mean(phen_subset$Amplitude_GCC_yearly, na.rm=TRUE),
                        Amplitude_GCC_peak_mean = mean(phen_subset$Amplitude_GCC_peak, na.rm=TRUE),
                        Cumulative_GCC_yearly_mean = mean(phen_subset$Cumulative_GCC_yearly, na.rm=TRUE),
                        Cumulative_GCC_to_EOS_mean = mean(phen_subset$Cumulative_GCC_to_EOS, na.rm=TRUE),
                        SOS_15_mean = mean(phen_subset$SOS_15, na.rm=TRUE),
                        SOS_25_mean = mean(phen_subset$SOS_25, na.rm=TRUE),
                        SOS_50_mean = mean(phen_subset$SOS_50, na.rm=TRUE),
                        Peak_mean = mean(phen_subset$Peak, na.rm=TRUE),
                        EOS_50_mean = mean(phen_subset$EOS_50, na.rm=TRUE),
                        EOS_25_mean = mean(phen_subset$EOS_25, na.rm=TRUE),
                        EOS_15_mean = mean(phen_subset$EOS_15, na.rm=TRUE),
                        Peak_GCC_std = sd(phen_subset$Peak_GCC, na.rm=TRUE),
                        SOS_15_std = sd(phen_subset$SOS_15, na.rm=TRUE),
                        SOS_25_std = sd(phen_subset$SOS_25, na.rm=TRUE),
                        SOS_50_std = sd(phen_subset$SOS_50, na.rm=TRUE),
                        Peak_std = sd(phen_subset$Peak, na.rm=TRUE),
                        EOS_50_std = sd(phen_subset$EOS_50, na.rm=TRUE),
                        EOS_25_std = sd(phen_subset$EOS_25, na.rm=TRUE),
                        EOS_15_std = sd(phen_subset$EOS_15, na.rm=TRUE),
                        Number_Years = nrow(phen_subset))
  growing_seasons_phen_ave <- rbind(growing_seasons_phen_ave,new_row)
}

rm(cumulative_gcc_yearly, cumulative_gcc, new_row, phen_subset, unique_phenocams, eos_15, eos_25, eos_50, 
   i, peak, peak_amplitude, peak_gcc, peak_gcc_baseline_subtracted, phenocam, roi, sos_15, sos_25, sos_50, veg_type, x, year, yearly_amplitude)

# Save results
save(growing_seasons_phen, site_year_statistics, growing_seasons_phen_ave, file = "outputs/growing_seasons_phen_2016_to_2024_baseline_subtracted.RData")
