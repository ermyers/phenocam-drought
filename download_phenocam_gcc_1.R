# Download phenocam GCC from previously created lists of phenocams

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(phenocamr)

# Load previous lists of phenocam data
load("outputs/phenocam_lists.RData")

for(sitename in list1){
  download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
}

list1b <- list1[74:110]
for(sitename in list1b){
  download_phenocam(site = paste(sitename,"$",sep=""), out_dir="data/phenocam-gcc")
}
