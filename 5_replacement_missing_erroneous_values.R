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

# 1) Reading the NASA files
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/NASA")
NASA_FILES<-list()
for (i in list.files(pattern = ".txt")) {
  NASA_FILES[[i]]<-read.table(i,header=T)
}

NASA_FILES<-do.call('rbind',NASA_FILES)
NASA_FILES$day<-lubridate::yday(NASA_FILES$date)

# 2) Correlation between field solar radiation values and NASA files
cor <- list()
for (j in unique(NASA_FILES$Year_Exp)) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == j, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == j, 'Date.Harvested'])
  ind <- day_start:day_end
  if (!all(is.na(daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2']))) {
    cor[[j]] <-
      cor(NASA_FILES[NASA_FILES$Year_Exp == j &
                       NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN'],
          daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2'],
          method = 'pearson',
          use = 'complete.obs')}
    else{
      cor[[j]] <- 'no field value'
    }
  
}

write.table(as.data.frame(cor),'cor.daily.solar.radiation.fieldstation.vs.NASA.txt',col.names=T,quote = F)

# 3) Replacement of all field values for which the correlation with NASA values is < 0.6

for (variable in c(names(which(cor<0.6)),names(which(cor=='no field value')))) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
  ind <- day_start:day_end
  daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2'] <- NASA_FILES[NASA_FILES$Year_Exp == variable &
                                                                                               NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN']
}

# 4) Replacement of missing values for other field trials - but only missing days.

for (variable in unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp)%notin%c(names(which(cor<0.6)),names(which(cor=='no field value')))]) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
  ind <- day_start:day_end
  ind_missing_days<-which(is.na(daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2']))
  daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2'][ind_missing_days] <- NASA_FILES[NASA_FILES$Year_Exp == variable &
                                                                                               NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN'][ind_missing_days]
}