rm(list=ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)
library(mapdata)
library(lubridate)
library(maps)
library(maptools)
library(xts)
library(spacetime)
library(rgdal)
source('fahrenheit_to_celsius.R')

`%notin%` <- Negate(`%in%`)


library(dplyr)
#library(plyr)
library(lubridate)


# ------------------------------------------------------------------------------
# Load dataset with missing values #
# ------------------------------------------------------------------------------

daily_weather=read.table(
  '2_merged_dataset.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)


# ------------------------------------------------------------------------------
# Solar radiation data
# ------------------------------------------------------------------------------

#1) Reading the files
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/NASA")
NASA_FILES<-list()
for (i in list.files(pattern = ".txt")) {
  NASA_FILES[[i]]<-read.table(i,header=T)
}

NASA_FILES<-do.call(rbind,NASA_FILES)
NASA_FILES$day<-lubridate::yday(NASA_FILES$date)

cor<-list()
for (j in unique(NASA_FILES$Year_Exp)) {
  day_start<-unique(daily_weather[daily_weather$Year_Exp==j,'Date.Planted'])
  day_end<-unique(daily_weather[daily_weather$Year_Exp==j,'Date.Harvested'])
  ind<-day_start:day_end
  if (!all(is.na(daily_weather[daily_weather$Year_Exp==j,'incoming_radiation_MJm2']))){cor[[j]]<-cor(NASA_FILES[NASA_FILES$Year_Exp==j&NASA_FILES$day%in%ind,'ALLSKY_SFC_SW_DWN'],daily_weather[daily_weather$Year_Exp==j,'incoming_radiation_MJm2'],method = 'pearson',use = 'complete.obs')

  }}
