# Download phenocam GCC from previously created lists of phenocams

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(phenocamr)

# Load previous lists of phenocam data
load("outputs/phenocam_lists.RData")

# for(sitename in list3){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }
# 
# list3b <- list3[10:110]
# for(sitename in list3b){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }
# 
# list3c <- list3b[15:101]
# for(sitename in list3c){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }
# 
# list3d <- list3c[6:87]
# for(sitename in list3d){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }
# 
# list3e <- list3d[6:82]
# for(sitename in list3e){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }
# 
# list3f <- list3e[4:77]
# for(sitename in list3f){
#   download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
# }

download_phenocam(site = "NEON", out_dir="data/phenocam-gcc")

list3b <- list3[66:102]
for(sitename in list3g){
  download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
}

download_phenocam(site = "segaarboretum", out_dir = "data/phenocam-gcc")
