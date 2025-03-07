# Extract USDM data at PhenoCam point locations

# USDM shapefiles are downloaded from the US Drought Monitor Webpage
# For example, the download link for all USDM shapefiles in 2024 is:
# https://droughtmonitor.unl.edu/data/shapefiles_m//2024_USDM_M.zip
# A sample command to download USDM files into the folder of interest would be:
# download.file("https://droughtmonitor.unl.edu/data/shapefiles_m//2024_USDM_M.zip", "data/usdm/2024/2024_USDM_M.zip")
# Once downloaded, the yearly file unzips into individual zipped shapefiles, which also need to be unzipped

# Code is parallelized by running it on one calendar year at a time (specified by year_list parameter)
# Output dataframes are later combined together into one large dataframe in load_phenocam_usdm_data.R

library(terra)
library(lubridate)

# Define a function to extract the date from the USDM filename
date_from_shp <- function(filename){
  date_text <- substring(sub(".*USDM_", "", filename),1,8)
  date <- as.Date(date_text, format = "%Y%m%d")
  date
}

# Load in the phenocam locations
type1_sites_usa_vect <- vect("outputs/shapefiles/type1_sites_usa.shp")
type1_sites_usa_df <- as.data.frame(type1_sites_usa_vect)

# Extract USDM values at different points
#year_list <- c("2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023","2024")
year_list <- c("2010")
phen_usdm <- data.frame(matrix(ncol=5,nrow=0))
colnames(phen_usdm) <- c("Phenocam","Date","Year","DOY","DM")

for (year in year_list){
  file_list <- list.files(paste("data/usdm/",year,sep=""), pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  for (file in file_list){
    date <- date_from_shp(file)
    print(date)
    year <- as.numeric(format(date, "%Y"))
    doy <- yday(date)
    usdm <- vect(file)
    usdm_at_phen <- extract(usdm,type1_sites_usa_vect)
    phen_usdm <- rbind(phen_usdm,
                       data.frame("Phenocam"=type1_sites_usa_df$site,
                                  "Date"=date,
                                  "Year"=year,
                                  "DOY"=doy,
                                  "DM"=usdm_at_phen$DM))
  }
}

phen_usdm_2010 <- phen_usdm
save(phen_usdm_2010, file="outputs/phen_usdm_2010.RData")
