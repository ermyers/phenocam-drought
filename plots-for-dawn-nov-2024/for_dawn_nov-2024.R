# Plots for Dawn's 11/11 presentation

# Load packages
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data
load("outputs/phen_gcc.RData")
load("outputs/drought_statistics.RData")

phen_usdm_modified <- mutate(phen_usdm_modified, DM=DM+1)

# IBP
# Last common date for IBP GCC and IBP USDM is 10/22/2024
# First common date for IBP GCC and IBP USDM is 05/14/2013
ibp_gcc <- ggplot(data=filter(phen_gcc,Phenocam=="ibp" & Veg_Type=="XX" & Date>="2013-05-14" & Date<="2024-10-22"), aes(x=Date,y=gcc_90)) +
  geom_point() +
  theme_bw() +
  ylab("GCC") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  ggtitle("IBP")
ibp_dm <- ggplot(data=filter(phen_usdm_modified,Phenocam=="ibp" & Date>="2013-05-14" & Date<="2024-10-22"), aes(x=Date,y=DM,fill=as.factor(DM))) +
  geom_col() +
  scale_fill_manual(values=c("#CCCCCC","#FFFF00","#FFCC66","#FF9900","#FF3300","#990000")) +
  theme_bw() +
  theme(legend.position="none") +
  ylab("Drought Index") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), labels = c("None","D0","D1","D2","D3","D4")) +
  expand_limits(y=5)
ibp_gcc / ibp_dm

# Jernort
# First common date for nort GCC and nort USDM is 03/18/2014
# Last common date for nort GCC and nort USDM is 10/22/2024
nort_gcc <- ggplot(data=filter(phen_gcc,Phenocam=="jernort" & Veg_Type=="XX" & Date>="2014-03-18" & Date<="2024-10-22"), aes(x=Date,y=gcc_90)) +
  geom_point() +
  theme_bw() +
  ylab("GCC") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  ggtitle("Jernort")
nort_dm <- ggplot(data=filter(phen_usdm_modified,Phenocam=="jernort" & Date>="2014-03-18" & Date<="2024-10-22"), aes(x=Date,y=DM,fill=as.factor(DM))) +
  geom_col() +
  scale_fill_manual(values=c("#CCCCCC","#FFFF00","#FFCC66","#FF9900","#FF3300","#990000")) +
  theme_bw() +
  theme(legend.position="none") +
  ylab("Drought Index") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), labels = c("None","D0","D1","D2","D3","D4")) +
  expand_limits(y=5)
nort_gcc / nort_dm

# Jerbajada
# First common date for bajada GCC and bajada USDM is 04/29/2014
# Last common date for bajada GCC and bajada USDM is 10/22/2024
bajada_gcc <- ggplot(data=filter(phen_gcc,Phenocam=="jerbajada" & Veg_Type=="SH" & Date>="2014-04-29" & Date<="2024-10-22"), aes(x=Date,y=gcc_90)) +
  geom_point() +
  theme_bw() +
  ylab("GCC") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  ggtitle("Jerbajada")
bajada_dm <- ggplot(data=filter(phen_usdm_modified,Phenocam=="jerbajada" & Date>="2014-04-29" & Date<="2024-10-22"), aes(x=Date,y=DM,fill=as.factor(DM))) +
  geom_col() +
  scale_fill_manual(values=c("#CCCCCC","#FFFF00","#FFCC66","#FF9900","#FF3300","#990000")) +
  theme_bw() +
  theme(legend.position="none") +
  ylab("Drought Index") +
  xlab("Year") +
  scale_x_date(breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), labels = c("None","D0","D1","D2","D3","D4")) +
  expand_limits(y=5)
bajada_gcc / bajada_dm