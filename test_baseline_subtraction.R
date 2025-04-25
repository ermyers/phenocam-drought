# Test code for baseline normalization

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(dplyr)
#library(pracma)
library(ggplot2)
#library(lubridate)
library(zoo)

# Load GCC data
load('outputs/phen_gcc.RData')

# Load list of manually filtered PhenoCam ROIs
unique_phenocam_rois_annotated <- read.csv("data/unique_phenocam_rois_annotated_03-09-2025.csv")

# Filter list to include only primary ROIs
primary_rois <- filter(unique_phenocam_rois_annotated, action=="keep" & rank==1)

# Select ibp and 9 other random ROIs from primary ROI list
sample_rois <- slice_sample(filter(primary_rois, Phenocam!="ibp" & Phenocam!="vaira"), n=8)
sample_rois <- rbind(sample_rois, filter(primary_rois, Phenocam=="ibp"))
sample_rois <- rbind(sample_rois, filter(primary_rois, Phenocam=="vaira"))
sample_rois$ROI <- as.character(sample_rois$ROI)

# Filter GCC data to only include the 10 selected ROIs
sample_gcc <- right_join(phen_gcc, select(sample_rois, Phenocam, Veg_Type, ROI), by=join_by(Phenocam, Veg_Type, ROI))

# Baseline normalization

# Find max in each year
sample_yearly_peaks <- data.frame()

for (phenocam in unique(sample_gcc$Phenocam)){
  phen_subset <- filter(sample_gcc, Phenocam==phenocam)
  
  for (year in unique(phen_subset$Year)){
    year_subset <- filter(phen_subset, Year==year)
    peak_gcc <- max(year_subset$smooth_gcc_90, na.rm=TRUE)
    
    # only record peak for year with at least 200 days of non-interpolated data
    if(nrow(filter(year_subset, is.na(int_flag)==TRUE))>=200){
      sample_yearly_peaks <- rbind(sample_yearly_peaks,
                                   data.frame(Phenocam = phenocam,
                                              Year = year,
                                              Peak_GCC = peak_gcc,
                                              Peak_Date = max(filter(year_subset, smooth_gcc_90==peak_gcc)$Date)))
    }
    
  }
}

rm(phen_subset, year_subset, peak_gcc, phenocam, year)

# Find minima between peaks
sample_minima <- data.frame()

for (phenocam in unique(sample_yearly_peaks$Phenocam)){
  phen_subset <- filter(sample_yearly_peaks, Phenocam==phenocam)
  
  # Look for minimum before first peak
  between_peak_subset <- filter(sample_gcc, Phenocam==phenocam & Date<phen_subset$Peak_Date[1])
  min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
  sample_minima <- rbind(sample_minima,
                          data.frame(Phenocam = phenocam,
                                     Min_GCC = min_gcc,
                                     Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
  # Look for minima between peaks
  for (i in 2:nrow(phen_subset)){
    between_peak_subset <- filter(sample_gcc, Phenocam==phenocam & Date>phen_subset$Peak_Date[i-1] & Date<phen_subset$Peak_Date[i])
    min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
    min_date <- filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date
    if (length(min_date)>1){
      # not sure
    }
    sample_minima <- rbind(sample_minima,
                           data.frame(Phenocam = phenocam,
                                      Min_GCC = min_gcc,
                                      Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
  }

  # Look for minimum after last peak
  between_peak_subset <- filter(sample_gcc, Phenocam==phenocam & Date > phen_subset$Peak_Date[nrow(phen_subset)])
  min_gcc <- min(between_peak_subset$smooth_gcc_90, na.rm=TRUE)
  sample_minima <- rbind(sample_minima,
                         data.frame(Phenocam = phenocam,
                                    Min_GCC = min_gcc,
                                    Min_Date = filter(between_peak_subset, smooth_gcc_90==min_gcc)$Date))
}

rm(between_peak_subset, phen_subset, phenocam, i, min_gcc, min_date)

# Interpolate between minima
sample_gcc_with_baseline <- data.frame()

for (phenocam in unique(sample_gcc$Phenocam)){
  phen_subset <- filter(sample_gcc, Phenocam==phenocam)
  phen_subset <- left_join(phen_subset, sample_minima, by=join_by(Phenocam, Date==Min_Date))
  phen_subset <- mutate(phen_subset, Baseline_GCC = na.approx(Min_GCC, rule=2))
  sample_gcc_with_baseline <- rbind(sample_gcc_with_baseline, phen_subset)
}

rm(phenocam)

# Subtract minima from GCC
sample_gcc_with_baseline <- mutate(sample_gcc_with_baseline, smooth_gcc_baseline_subtracted = smooth_gcc_90-Baseline_GCC)

# Plot results

ggplot() +
  geom_point(data=filter(sample_gcc, Phenocam=="vaira"), aes(x=Date,y=smooth_gcc_90)) +
  geom_line(data=filter(sample_gcc_with_baseline, Phenocam=="vaira"), aes(x=Date,y=Baseline_GCC), color='red') +
  ggtitle("vaira")

ggplot() +
  geom_point(data=filter(sample_gcc_with_baseline, Phenocam=="vaira"), aes(x=Date,y=smooth_gcc_baseline_subtracted)) +
  ggtitle("vaira")
