# Extract USDM data at PhenoCam point locations

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
year_list <- c("2024")
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

phen_usdm_2024 <- phen_usdm
save(phen_usdm_2024, file="outputs/phen_usdm_2024.RData")