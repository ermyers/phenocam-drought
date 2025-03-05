# Create dataframes with PhenoCam GCC data and USDM data
#
# Phenocam GCC (phen_gcc): Includes all GCC values for downloaded phenocam data
# (which is limited to Type I USA sites with at least 4 years of data), along with
# ecoregion data for each phenocam.
# Phenocam USDM (phen_usdm): Includes USDM values for all Type I USA Phenocam sites
# between the years of 2010-2024, along with ecoregion data for each phenocam.

library(stringr)
library(dplyr)

######
# GCC
######

gcc_list <- list.files("data/phenocam-gcc")
load("outputs/phenocam_lists.RData")

# Create dataframe to hold phenocam GCC
phen_gcc <- data.frame(matrix(ncol=9,nrow=0))
colnames(phen_gcc) <- c("Phenocam","Veg_Type","ROI","Date","Year","DOY","gcc_90","smooth_gcc_90","int_flag")

# For now, filter GCC list to only include Type I sites in USA
for (i in 1:length(gcc_list)){
  phenocam <- str_split(gcc_list[i],pattern = "_",simplify=TRUE)[1]
  veg_type <- str_split(gcc_list[i],pattern = "_",simplify=TRUE)[2]
  roi <- str_split(gcc_list[i],pattern = "_",simplify=TRUE)[3]
  
  if (phenocam %in% type1_sitenames_usa){
    gcc_df <- read.csv(paste("data/phenocam-gcc/",gcc_list[i],sep=""),skip=24)
    gcc_df <- gcc_df %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
    gcc_df <- data.frame("Phenocam"=phenocam, "Veg_Type"=veg_type, "ROI"=roi,
                         "Date"=gcc_df$date, "Year"=gcc_df$year, "DOY"=gcc_df$doy,
                         "gcc_90"=gcc_df$gcc_90, "smooth_gcc_90"=gcc_df$smooth_gcc_90,
                         "int_flag"=gcc_df$int_flag)
    phen_gcc <- rbind(phen_gcc,gcc_df)
  }
}

# Optional - remove unneeded variables
rm(gcc_df, gcc_list, i, phenocam, veg_type, roi)

# Add level I ecoregion data to GCC
load("outputs/phen_with_ecoregion.RData")
phen_with_ecoregion_df <- select(phen_with_ecoregion_df, "site", "ecoregion", "NA_L1CODE", "NA_L1NAME")
phen_gcc <- left_join(phen_gcc,phen_with_ecoregion_df,by=join_by("Phenocam"=="site"))

# Save GCC
save(phen_gcc, file = "outputs/phen_gcc.RData")

#######
# USDM
#######

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

phen_usdm <- rbind(phen_usdm_2010,phen_usdm_2011,phen_usdm_2012,phen_usdm_2013,
                   phen_usdm_2014,phen_usdm_2015,phen_usdm_2016,phen_usdm_2017,
                   phen_usdm_2018,phen_usdm_2019,phen_usdm_2020,phen_usdm_2021,
                   phen_usdm_2022,phen_usdm_2023,phen_usdm_2024)

phen_usdm[is.na(phen_usdm$DM),]$DM <- -1

# Optional - remove yearly USDM
rm(phen_usdm_2010,phen_usdm_2011,phen_usdm_2012,phen_usdm_2013,phen_usdm_2014,
   phen_usdm_2015,phen_usdm_2016,phen_usdm_2017,phen_usdm_2018,phen_usdm_2019,
   phen_usdm_2020,phen_usdm_2021,phen_usdm_2022,phen_usdm_2023,phen_usdm_2024)

# Add level I ecoregion data to USDM
load("outputs/phen_with_ecoregion.RData")
phen_with_ecoregion_df <- select(phen_with_ecoregion_df, "site", "ecoregion", "NA_L1CODE", "NA_L1NAME")
phen_usdm <- left_join(phen_usdm,phen_with_ecoregion_df,by=join_by("Phenocam"=="site"))

# Save USDM
save(phen_usdm, file = "outputs/phen_usdm.RData")