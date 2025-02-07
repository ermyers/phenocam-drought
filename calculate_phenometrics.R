# Calculate phenometrics

# Note 10/29/24: Code to eliminate smaller peaks with larger peaks within 30 days
# throws a warning for peaks at very beginning and end of time series (but should
# still return good results - these peaks will be eliminated at the check for 35% variation)

# Load packages
library(dplyr)
library(pracma)
library(ggplot2)
library(lubridate)

# Load GCC data
load('outputs/phen_gcc.RData')

#######################################
# Calculate season start and end dates
#######################################

growing_seasons_phen <- data.frame(matrix(ncol = 21, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Year","Peak_GCC","Cumulative_GCC_yearly","Cumulative_GCC","SOS_15","SOS_25","SOS_50","Peak","EOS_50","EOS_25","EOS_15",
       "SOS_15_flag","SOS_25_flag","SOS_50_flag","Peak_flag","EOS_50_flag","EOS_25_flag","EOS_15_flag")
colnames(growing_seasons_phen) <- x

# Loop over each phenocam
unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)
# unique_phenocams <- filter(phen_gcc, Phenocam=="ibp") %>% distinct(Phenocam,Veg_Type,ROI) # test with just one
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Loop over each year, plus or minus some buffer
  for (year in unique(phen_subset$Year)){
    date_min <- head(filter(phen_subset, Year==year),n=1)$Date - 185
    date_max <- tail(filter(phen_subset, Year==year),n=1)$Date + 185
    date_subset <- filter(phen_subset, Date>=date_min & Date<= date_max)
    
    min_GCC <- min(date_subset$smooth_gcc_90,na.rm=TRUE)
    max_GCC <- max(date_subset$smooth_gcc_90,na.rm=TRUE)
    candidate_peaks <- findpeaks(date_subset$smooth_gcc_90, peakpat = "[+]{1,}[0]*[-]{1,}")
    candidate_peaks_chrono <- date_subset[candidate_peaks[,2],] # extract rows containing candidate peaks
    all_candidate_peaks_chrono <- candidate_peaks_chrono
    candidate_peaks_sorted <- dplyr::arrange(candidate_peaks_chrono,smooth_gcc_90) # candidate peaks sorted in ascending order of value
    all_candidate_peaks_sorted <- candidate_peaks_sorted
    
    # print statements
    #print(" ")
    #print(paste("year = ",year,", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
    #print(paste("candidate peaks chrono"))
    #print(paste(candidate_peaks_chrono$Date))
    #print(paste("candidate peaks sorted"))
    #print(paste(candidate_peaks_sorted$Date))
    
    # Check for a nonzero number of candidate peaks
    if(nrow(all_candidate_peaks_sorted)<1){
    } else{
      
      # Loop over each candidate peak and eliminate peaks that don't meet the criteria,
      # moving from smallest to largest
      for (j in 1:nrow(all_candidate_peaks_sorted)){
        date <- all_candidate_peaks_sorted$Date[j]
        peak_index <- which(candidate_peaks_chrono$Date==date)
        
        # print statements
        #print(paste("j = ",j,", date = ",date,", peak index = ",peak_index,sep=""))
        #print(paste(candidate_peaks_chrono$Date))
        
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
        
        # print statements
        # print(paste("dist to prev = ",dist_to_prev,", dist to next = ",dist_to_next,sep=""))
        
        # Eliminate small peaks with larger peaks within 30 days
        if(dist_to_prev<30 | dist_to_next<30){
          candidate_peaks_sorted <- candidate_peaks_sorted[candidate_peaks_sorted$Date != date,]
          candidate_peaks_chrono <- candidate_peaks_chrono[candidate_peaks_chrono$Date != date,]
        } else{
          peak_subset <- filter(phen_subset, Date>=date-dist_to_prev & Date<=date+dist_to_next)
          pre_peak_min <- min(filter(peak_subset,Date<=date-30)$smooth_gcc_90,na.rm=TRUE)
          pre_peak_min_date <- min(filter(peak_subset,smooth_gcc_90==pre_peak_min)$Date,na.rm=TRUE)
          post_peak_min <- min(filter(peak_subset, Date>=date+30)$smooth_gcc_90,na.rm=TRUE)
          post_peak_min_date <- max(filter(peak_subset,smooth_gcc_90==post_peak_min)$Date,na.rm=TRUE)
          peak_value <- filter(peak_subset, Date==date)$smooth_gcc_90
          
          # print statements
          # if(pre_peak_min == Inf | post_peak_min == Inf){
          #   print(paste("Date =",date,"Peak =",peak_value))
          #   print(paste("pre peak min =",pre_peak_min,"post peak min =",post_peak_min))
          #   print(paste("pre peak min date =",pre_peak_min_date,"post peak min date =",post_peak_min_date))
          # }
          # print(paste("pre peak min =",pre_peak_min,"post peak min =",post_peak_min))
          # print(paste("pre peak min date =",pre_peak_min_date,"post peak min date =",post_peak_min_date))
          # print(paste("peak value =",peak_value))
          
          # Eliminate peaks with difference between minimum and maximum <35% of total variation
          if((peak_value-pre_peak_min)/(max_GCC-min_GCC)<0.35 | (peak_value-post_peak_min)/(max_GCC-min_GCC)<0.35){
            candidate_peaks_sorted <- candidate_peaks_sorted[candidate_peaks_sorted$Date != date,]
            candidate_peaks_chrono <- candidate_peaks_chrono[candidate_peaks_chrono$Date != date,]
            
            # print statements
            # print(paste("% of total variation =",(peak_value-pre_peak_min)/(max_GCC-min_GCC),(peak_value-post_peak_min)/(max_GCC-min_GCC)))

          } else{
            
            # Only consider peak values that are in the calendar year
            if(peak_subset[peak_subset$Date==date,]$Year == year){
              
              # Calculate peak timing and season start and end dates, based on % thresholds
              SOS_15_thresh <- 0.15*(peak_value-pre_peak_min) + pre_peak_min
              SOS_25_thresh <- 0.25*(peak_value-pre_peak_min) + pre_peak_min
              SOS_50_thresh <- 0.50*(peak_value-pre_peak_min) + pre_peak_min
              EOS_50_thresh <- 0.50*(peak_value-post_peak_min) + post_peak_min
              EOS_25_thresh <- 0.25*(peak_value-post_peak_min) + post_peak_min
              EOS_15_thresh <- 0.15*(peak_value-post_peak_min) + post_peak_min
              SOS_15 <- min(filter(peak_subset, smooth_gcc_90>SOS_15_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              SOS_25 <- min(filter(peak_subset, smooth_gcc_90>SOS_25_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              SOS_50 <- min(filter(peak_subset, smooth_gcc_90>SOS_50_thresh & Date >= pre_peak_min_date)$Date, na.rm=TRUE)
              EOS_50 <- max(filter(peak_subset, smooth_gcc_90>EOS_50_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              EOS_25 <- max(filter(peak_subset, smooth_gcc_90>EOS_25_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              EOS_15 <- max(filter(peak_subset, smooth_gcc_90>EOS_15_thresh & Date <= post_peak_min_date)$Date, na.rm=TRUE)
              date_peak <- peak_subset[peak_subset$Date==date,]$Date
              
              # Calculate cumulative GCC (start of year to EOS 25) and cumulative yearly GCC (start to end of year)
              year_subset <- filter(phen_subset, Year==year)
              if(min(year_subset$DOY, na.rm=TRUE)!=1){
                cumulative_gcc <- NA
                cumulative_gcc_yearly <- NA
              } else{
                cumulative_gcc <- sum(filter(year_subset, Date <= EOS_25)$smooth_gcc_90, na.rm=TRUE)
                if(max(year_subset$DOY, na.rm=TRUE)<365){
                  cumulative_gcc_yearly <- NA
                  } else{
                  cumulative_gcc_yearly <- sum(year_subset$smooth_gcc_90, na.rm=TRUE)
                  }
                }
              
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
                                    Peak_GCC = peak_value,
                                    Cumulative_GCC_yearly = cumulative_gcc_yearly,
                                    Cumulative_GCC = cumulative_gcc,
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
              
              # print statements
              # print("new row")
              # print(paste(new_row))
            } else{}
          }
        }
      }
    }
  }
}

# Optional - remove variables
rm(all_candidate_peaks_chrono,all_candidate_peaks_sorted,candidate_peaks,candidate_peaks_chrono,
   candidate_peaks_sorted,date_subset,new_row,peak_subset,phen_subset,unique_phenocams,year_subset)
rm(cumulative_gcc,cumulative_gcc_yearly,date,date_max,date_min,date_peak,dist_to_next,dist_to_prev,EOS_15,EOS_15_flag,
   EOS_15_thresh,EOS_25,EOS_25_flag,EOS_25_thresh,EOS_50,EOS_50_flag,EOS_50_thresh,
   i,j,max_GCC,min_GCC,peak_flag,peak_index,peak_value,phenocam,post_peak_min,
   post_peak_min_date,pre_peak_min,pre_peak_min_date,roi,SOS_15,SOS_15_flag,SOS_15_thresh,
   SOS_25,SOS_25_flag,SOS_25_thresh,SOS_50,SOS_50_flag,SOS_50_thresh,veg_type,x,year)

# Optional - plot selected results
ggplot() +
  geom_point(data=filter(phen_gcc, Phenocam=="ibp" & Veg_Type=="XX"), aes(x=Date,y=smooth_gcc_90)) +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp" & Veg_Type=="XX")$SOS_15, color="yellow") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp"& Veg_Type=="XX")$EOS_15, color="brown") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp" & Veg_Type=="XX")$SOS_50, color="green") +
  geom_vline(xintercept = filter(growing_seasons_phen, Phenocam=="ibp"& Veg_Type=="XX")$EOS_50, color="red")

##########################################################
# Identify site-years with no or multiple growing seasons
##########################################################

growing_seasons_phen_number <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Year","Number_Growing_Seasons","Peak_GCC")
colnames(growing_seasons_phen_number) <- x

unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)
for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- phen_gcc %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  # print statements
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  for (year in unique(phen_subset$Year)){
    # check that we have non-interpolated data for most of the year
    if(nrow(filter(phen_subset,Year==year & is.na(int_flag)==TRUE))>=200){
      number_growing_seasons <- nrow(filter(growing_seasons_phen, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year))
      peak_gcc <- max(filter(phen_gcc, Phenocam==phenocam & Veg_Type==veg_type & ROI==roi & Year==year)$smooth_gcc_90, na.rm=TRUE)
      new_row <- data.frame(Phenocam = phenocam,
                            Veg_Type = veg_type,
                            ROI = roi,
                            Year = year,
                            Number_Growing_Seasons = number_growing_seasons,
                            Peak_GCC = peak_gcc)
      growing_seasons_phen_number <- rbind(growing_seasons_phen_number,new_row)
    }
  }
}

# Optional - remove variables
rm(new_row, phen_subset,unique_phenocams,i,number_growing_seasons,peak_gcc,phenocam,roi,veg_type,x,year)

###############################################################
# Identify typical growing season values for each phenocam ROI
###############################################################

growing_seasons_phen_ave <- data.frame(matrix(ncol = 22, nrow = 0))
x <- c("Phenocam","Veg_Type","ROI","Peak_GCC_mean","Cumulative_GCC_yearly_mean","Cumulative_GCC_mean","SOS_15_mean","SOS_25_mean",
       "SOS_50_mean","Peak_mean","EOS_50_mean","EOS_25_mean","EOS_15_mean",
       "Peak_GCC_std","SOS_15_std","SOS_25_std","SOS_50_std","Peak_std",
       "EOS_50_std","EOS_25_std","EOS_15_std","Number_Years")
colnames(growing_seasons_phen_ave) <- x

unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)

for (i in 1:nrow(unique_phenocams)){
  phenocam <- unique_phenocams$Phenocam[i]
  veg_type <- unique_phenocams$Veg_Type[i]
  roi <- unique_phenocams$ROI[i]
  phen_subset <- growing_seasons_phen %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi)
  phen_subset <- mutate(phen_subset, SOS_15 = yday(SOS_15), SOS_25 = yday(SOS_25), SOS_50 = yday(SOS_50), Peak = yday(Peak), EOS_50 = yday(EOS_50), EOS_25 = yday(EOS_25), EOS_15 = yday(EOS_15))
  print(paste("i = ",i," out of ",nrow(unique_phenocams),", phenocam = ",phenocam,", veg_type = ",veg_type," roi = ",roi,sep=""))
  
  # Combine multiple growing seasons if necessary
  for(year in unique(phen_subset$Year)){
    if (nrow(filter(phen_subset, Year==year)) > 1){
      peak_gcc <- max(filter(phen_subset, Year==year)$Peak_GCC, na.rm = TRUE)
      cumulative_gcc_yearly <- max(filter(phen_subset, Year==year)$Cumulative_GCC_yearly, na.rm = TRUE)
      cumulative_gcc <- max(filter(phen_subset, Year==year)$Cumulative_GCC, na.rm = TRUE)
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
                            Cumulative_GCC_yearly = cumulative_gcc_yearly,
                            Cumulative_GCC = cumulative_gcc,
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
                            EOS_15_flag = filter(phen_subset, Year==year & EOS_15==eos_15)$EOS_15_flag)
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
                        Cumulative_GCC_yearly_mean = mean(phen_subset$Cumulative_GCC_yearly, na.rm=TRUE),
                        Cumulative_GCC_mean = mean(phen_subset$Cumulative_GCC, na.rm=TRUE),
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
   i, peak, peak_gcc, phenocam, roi, sos_15, sos_25, sos_50, veg_type, x, year)

# Optional - save results
save(growing_seasons_phen, growing_seasons_phen_number, growing_seasons_phen_ave, file = "outputs/growing_seasons_phen.RData")