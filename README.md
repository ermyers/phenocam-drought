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
| Standardized Precipitation and Evapotranspiration Index (SPEI) | Weekly SPEI values at phenocam locations, provided by the National Drought Mitigation Center and derived from 5 KM raster data. | https://drought.unl.edu/ | N/A |
| Aridity Index | Average of annual aridity index for the years 1970-2000, provided at 30 arc second resolution. | [Future Global Aridity Index and PET Database (CMIP_6)](https://doi.org/10.57760/sciencedb.nbsdc.00086) | [1970-2000 aridity and PET](https://cstr.cn/16666.47.sciencedb.nbsdc.00086.0021605F.V2.V3.V4.V5.V6.V7)|

## Getting started

### Download shapefiles
The U.S. county shapefile is too large for github and needs to be downloaded locally. USDM and Level I Ecoregions shapefiles should already exist in the correct folders, but can be manually re-downloaded if needed.
- U.S. county shapefile should be downloaded and unzipped into data/shapefiles. (File path: data/shapefiles/tl_2023_us_county.shp)
- USDM shapefiles should be located in data/usdm/YEAR, with years between 2010 and 2024. (Example file path: data/usdm/2024/USDM_20241231.shp)
- Level I ecoregions shapefile should be downloaded and unzipped into data/ecoregions. (File path: data/ecoregions/NA_CEC_Eco_Level1.shp)

### Download aridity data
The global aridity index tif is too large for github and needs to be downloaded locally.
- From the download link for Aridity Index, WorldClim_2_1_1970-2000.zip should be downloaded and unzipped into data/WorldClim_2_1_1970-2000. (File path: data/shapefiles/WorldClim_2_1_1970-2000/aridity_index.tif)

### Run data download and preparation scripts
**create_phenocam_list.R** generates sitename lists and point location shapefiles of Type I PhenoCams located within the United States with at least 3, 4, and 5 years of GCC data, and saves them in outputs/phenocam_lists.RData and outputs/shapefiles. Certain unsuitable sites (e.g. NEON understory, spruce temperature- or rainfall-manipulated plots) are removed from the analysis at this stage.

**download_phenocam_gcc.R** downloads all phenocam GCC from the 4 year list into data/phenocam-gcc.

**extract_ecoregions_at_phen_locations.R** extracts level I ecoregion values and aridity index values for all PhenoCam locations, and saves this information in outputs/phen_with_ecoregion.RData.

**extract_usdm_at_phen_locations.R** extracts USDM values for a given year (given by year_list) for all PhenoCam locations, and saves these values in outputs/phen_usdm_YEAR.RData. This code should be run for each year of interest.

**load_phen_usdm_data.R** reads downloaded phenocam GCC into a dataframe, combines yearly USDM values into a single dataframe, and adds level I ecoregion information to both of these dataframes. Outputs are saved as outputs/phen_gcc.RData and outputs/phen_usdm.RData.

**subtract_phenocam_baseline.R** calculates baseline-subtracted GCC values for each phenocam site. Outputs are saved as outputs/phen_gcc_with_baseline.RData. The baseline is calculated by identifying the highest peak in each year, identifying the minimum GCC values between peaks, and linearly interpolating between minima.

## Calculate relevant metrics

### Calculate phenological metrics
**calculate_phenometrics_baseline_subtracted.R** calculates start of season (SOS), timing and value of peak GCC (actual peak and baseline-subtracted peak values are both recorded), and end of season (EOS) for each valid growing season for each site and year of PhenoCam data. SOS and EOS are calculated at multiple thresholds (15%, 25%, and 50% of the peak amplitude), and it is possible for sites to have no valid growing seasons or more than 1 valid growing season in a given year. In addition, yearly metrics like maximum GCC amplitude and cumulative GCC are calculated for all site-years, regardless of number of growing seasons. Outputs are saved as outputs/growing_seasons_phen.RData (all years of available PhenoCam data) and outputs/growing_seasons_phen_2016_to_2024.RData (data collected between 2016 and 2024).

### Calculate drought metrics
**calculate_vpd_anomaly.R** calculates the mean and cumulative VPD values for each site-year of PhenoCam data, and calculates the VPD anomalies per site-year as the difference between the mean for that year and the mean across all years for that site. Outputs are saved as outputs/vpd_statistics.RData.

**calculate_spei_anomaly.R** calculates the mean SPEI values for each site-year of PhenoCam data, and saves them in outputs/spei_statistics.RData.

**calculate_drought_statistics.R** calculates the start, end, number of weeks, and cumulative severity of continuous droughts (USDM >= 1) at PhenoCam sites. Outputs are saved in outputs/drought_statistics.RData.

## Analysis

### Manual Filtering and Site Selection
Some steps of the site selection are not easily automated and need to be performed manually. These include:
- In cases where a single PhenoCam had multiple valid ROIs (e.g. covering different vegetation types or different portions of the image), only one ROI was manually selected for analysis, with priority given to ROIs covering larger spatial or temporal extents.
- In cases where multiple cameras were positioned very close together (define how close?), only one ROI from the entire group of cameras was used for analysis.
- Landscape (XX) and some agricultural (AG) ROIs were re-labeled to match the primary vegetation type captured in the ROI.

To perform manual filtering, **list_phenocam_rois.R** was used to generate a CSV of all PhenoCam ROIs included in the growing season metrics calculation. Each ROI in the CSV was then manually annotated (based on visual inspection of ROIs) with a designation of "keep" or "remove", and a "rank" (1 for primary ROI, 2 or 3 for secondary or tertiary ROIs from the same camera or location). The annotated CSV was uploaded as data/unique_phenocam_rois_annotated_04-30-2025.csv and used in the analysis to consider only primary ROIs.

### Perform statistical analysis
**drought_analysis.R** calculates weighted average drought metric by site-year and performs linear regression between phenological metrics and drought metrics from USDM, VPD, and SPEI.
