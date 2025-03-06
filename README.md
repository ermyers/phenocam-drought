# phenocam-drought

## Overview

This project explores the relationship between PhenoCam GCC metrics and drought.

## Datasets

| Name | Description | Website | Download Information |
| ---- | ----------- | ------- | -------------------- |
| PhenoCam Green Chromatic Coordinate (GCC) | Smoothed daily GCC values from all PhenoCams meeting a specific set of requirements (located within the contiguous USA, Type I cameras, at least 5 years of data between 2016 and 2024, and several other requirements detailed below) | https://phenocam.nau.edu/webcam/ | [phenocamr package for downloading PhenoCam data using R](https://cran.r-project.org/web/packages/phenocamr/vignettes/phenocamr-vignette.html) |
| U.S. Drought Monitor (USDM) | Weekly map of drought conditions in the United States, with 6 possible categories ranging from no drought to severe drought. Data are available in multiple formats, but this study extracts USDM information from weekly vector shapefiles. | https://droughtmonitor.unl.edu/ | [USDM GIS data](https://droughtmonitor.unl.edu/DmData/GISData.aspx) |
| Level I Ecoregions of North America | 15 broad ecological regions in North America, which group areas of similar ecosystem type and natural resource availability. Data are available in vector shapefile format. | https://www.epa.gov/eco-research/ecoregions-north-america | [Direct download link](https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip) |
| U.S. Counties | 2023 U.S. county boundaries, provided in shapefile format by the U.S. Census Bureau. | https://www.census.gov/geographies/mapping-files.html | [Direct download link](https://www2.census.gov/geo/tiger/TIGER2023/COUNTY/tl_2023_us_county.zip) |
| Vapor Pressure Deficit (VPD) | Daily VPD values at phenocam locations, provided by the National Drought Mitigation Center (NDMC) and derived from 5 KM raster data. | https://drought.unl.edu/ | N/A |

## Getting started

### Download shapefiles
U.S. county shapefile is too large for github and needs to be downloaded locally. USDM and Level I Ecoregions shapefiles should already exist in the correct folders, but can be manually re-downloaded if needed.
- U.S. county shapefile should be downloaded and unzipped into data/shapefiles. (File path: data/shapefiles/tl_2023_us_county.shp)
- USDM shapefiles should be located in data/usdm/YEAR, with years between 2010 and 2024. (Example file path: data/usdm/2024/USDM_20241231.shp)
- Level I ecoregions shapefile should be downloaded and unzipped into data/ecoregions. (File path: data/ecoregions/NA_CEC_Eco_Level1.shp)

### Run data download and preparation scripts
**create_phenocam_list.R** generates sitename lists and point location shapefiles of Type I PhenoCams located within the United States with at least 3, 4, and 5 years of GCC data, and saves them in outputs/phenocam_lists.RData and outputs/shapefiles. Certain unsuitable sites (e.g. NEON understory, spruce temperature- or rainfall-manipulated plots) are removed from the analysis at this stage.

**download_phenocam_gcc.R** downloads all phenocam GCC from the 4 year list into data/phenocam-gcc.

**extract_ecoregions_at_phen_locations.R** extracts level I ecoregion values for all PhenoCam locations, and saves this information in outputs/phen_with_ecoregion.RData.

**extract_usdm_at_phen_locations.R** extracts USDM values for a given year (given by year_list) for all PhenoCam locations, and saves these values in outputs/phen_usdm_YEAR.RData. This code should be run for each year of interest.

**load_phen_usdm_data.R** reads downloaded phenocam GCC into a dataframe, combines yearly USDM values into a single dataframe, and adds level I ecoregion information to both of these dataframes. Outputs are saved as outputs/phen_gcc.RData and outputs/phen_usdm.RData.
