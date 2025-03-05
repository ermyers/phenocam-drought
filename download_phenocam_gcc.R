# Download phenocam GCC from previously created list of phenocams

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(phenocamr)

# Load previous lists of phenocam data
load("outputs/phenocam_lists.RData")

# list <- type1_sitenames_5years_usa[-c(125, 127, 128, 137, 139, 142, 149, 152,
#                                       159, 177, 180, 183, 185, 187, 188, 190,
#                                       194, 269)]
# remove sites with no files available for download:
# NEON.D01.HOPB.DP1.20002, NEON.D02.LEWI.DP1.20002, NEON.D02.POSE.DP1.20002,
# NEON.D04.CUPE.DP1.20002, NEON.D05.LIRO.DP1.20002, NEON.D06.MCDI.DP1.20002,
# NEON.D07.LECO.DP1.20002, NEON.D08.MAYF.DP1.20002, NEON.D13.COMO.DP1.20002,
# NEON.D13.WLOU.DP1.20002, NEON.D14.SYCA.DP1.20002, NEON.D15.REDB.DP1.20002,
# NEON.D16.MART.DP1.20002, NEON.D16.MCRA.DP1.20002, NEON.D17.BIGC.DP1.20002,
# NEON.D17.TECR.DP1.20002, spruceA0P21

list <- type1_sitenames_4years_usa[-c(68, 146, 148, 149, 158, 160, 163, 170, 173,
                                      180, 198, 201, 204, 206, 208, 209, 211,
                                      215, 297)]

# remove sites with no files available for download:
# flagstaff2, NEON.D01.HOPB.DP1.20002, NEON.D02.LEWI.DP1.20002, NEON.D02.POSE.DP1.20002,
# NEON.D04.CUPE.DP1.20002, NEON.D04.GUIL.DP1.20002, NEON.D05.LIRO.DP1.20002,
# NEON.D06.MCDI.DP1.20002, NEON.D07.LECO.DP1.20002, NEON.D08.MAYF.DP1.20002,
# NEON.D13.COMO.DP1.20002, NEON.D13.WLOU.DP1.20002, NEON.D14.SYCA.DP1.20002,
# NEON.D15.REDB.DP1.20002, NEON.D16.MART.DP1.20002, NEON.D16.MCRA.DP1.20002,
# NEON.D17.BIGC.DP1.20002, NEON.D17.TECR.DP1.20002, spruceA0P21

for(sitename in list){
  download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
}
