# Load in GCC and USDM and compute modified VCI based on GCC

library(stringr)
library(dplyr)

gcc_list <- list.files("data/phenocam-gcc")
load("outputs/phenocam_lists.RData")

# Create dataframe to hold phenocam GCC
phen_gcc <- data.frame(matrix(ncol=9,nrow=0))
colnames(phen_gcc) <- c("Phenocam","Veg_Type","ROI","Date","Year","DOY","gcc_90","smooth_gcc_90","int_flag")

# Had to delete flagstaff2_GR_2000, spruceT6P16E_DN_0001

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

# Calculate modified VCI (typically NDVI, but we will use GCC)
# Modified VCI = 100*(GCC - GCC_min)/(GCC_max - GCC_min)
# GCC = current GCC (or average over current time period)
# GCC_min = minimum historical GCC over the same time period
# GCC_max = maximum historical GCC over the same time period

# vci_weekly <- data.frame(matrix(ncol=9,nrow=0))
# colnames(vci_weekly) <- c("Phenocam","Veg_Type","ROI","Date","Year","DOY","GCC_weekly","VCI_weekly","USDM_weekly")
# 
# unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI)
# #unique_phenocams <- phen_gcc %>% distinct(Phenocam,Veg_Type,ROI) %>% filter(Phenocam=="harvard")
# for (i in 1:nrow(unique_phenocams)){
#   phenocam <- unique_phenocams$Phenocam[i]
#   veg_type <- unique_phenocams$Veg_Type[i]
#   roi <- unique_phenocams$ROI[i]
#   temp_df <- phen_gcc %>% filter(Phenocam == phenocam, Veg_Type == veg_type, ROI == roi, is.na(int_flag) == TRUE)
#   
#   # Check for at least 5 years of data
#   if (max(temp_df$Date) - min(temp_df$Date) >= 1826){
#     
#     # Filter drought monitor data to only look at phenocam of interest
#     temp_dm <- filter(phen_usdm, Phenocam==phenocam, Date>=min(temp_df$Date), Date<=max(temp_df$Date))
#     
#     # Iterate over each row of filtered drought monitor data
#     for (j in 1:nrow(temp_dm)){
#       date <- temp_dm$Date[j]
#       dm <- temp_dm$DM[j]
#       year <- temp_dm$Year[j]
#       doy <- temp_dm$DOY[j]
#       datelist <- filter(temp_df, DOY==doy)$Date # list all years with available data on current date
#       temp_gcc <- filter(temp_df, Date>(date-7), Date<=date) # only the GCC values in the week leading up to current date
#       
#       # Check that there are at least 4 days of GCC values
#       # and that there are at least 5 years with available data on current date
#       # Assuming both conditions are met, calculate mean GCC over that week of data
#       if(nrow(temp_gcc)>3 & length(datelist)>=5){
#         mean_gcc <- mean(temp_gcc$smooth_gcc_90) # mean weekly GCC
#         mean_gcc_list <- 0 # initialize mean GCC list - the 0 is a placeholder that will be removed later
#         
#         # Calculate mean GCC for each year with available data on current date
#         for (temp_date in datelist){
#           temp_gcc <- filter(temp_df, Date>(temp_date-7), Date<=temp_date)
#           # Check that there are at least 4 days of GCC values
#           # If yes, add the mean GCC to the mean GCC list
#           if(nrow(temp_gcc)>3){
#             mean_gcc_list <- c(mean_gcc_list,mean(temp_gcc$smooth_gcc_90))
#           }
#         }
#         mean_gcc_list <- mean_gcc_list[-1] # remove placeholder element
#         
#         # Check if there are at least 5 years of weekly GCC data
#         # If yes, calculate min and max GCC and use for VCI
#         if (length(mean_gcc_list) >= 5){
#           max_gcc <- max(mean_gcc_list)
#           min_gcc <- min(mean_gcc_list)
#           vci <- 100*(mean_gcc-min_gcc)/(max_gcc-min_gcc)
#           vci_weekly <- rbind(vci_weekly,
#                               data.frame("Phenocam"=phenocam,
#                                          "Veg_Type"=veg_type,
#                                          "ROI"=roi,
#                                          "Date"=date,
#                                          "Year"=year,
#                                          "DOY"=doy,
#                                          "GCC_weekly"=mean_gcc,
#                                          "VCI_weekly"=vci,
#                                          "USDM_weekly"=dm))
#         }
#       }
#     }
#   }
# }
# 
# # Optional - remove unneeded variables
# rm(temp_df,temp_dm,temp_gcc,unique_phenocams,date,datelist,dm,doy,i,j,
#    max_gcc,mean_gcc,mean_gcc_list,min_gcc,phenocam,roi,temp_date,vci,veg_type,year)
# 
# # Add level I ecoregion data to VCI
# load("outputs/phen_with_ecoregion.RData")
# phen_with_ecoregion_df <- select(phen_with_ecoregion_df, "site", "ecoregion", "NA_L1CODE", "NA_L1NAME")
# vci_weekly <- left_join(vci_weekly,phen_with_ecoregion_df,by=join_by("Phenocam"=="site"))
# 
# # Save VCI
# save(vci_weekly, file = "outputs/vci_weekly.RData")
