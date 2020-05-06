rm(list=ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)
library(mapdata)
library(maps)
library(maptools)
library(xts)
library(spacetime)
library(rgdal)
source('fahrenheit_to_celsius.R')
source('safeguarding.R')
`%notin%` <- Negate(`%in%`)


library(dplyr)
#library(plyr)
library(lubridate)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")

stations <- ghcnd_stations() 
stations<-filter(stations,last_year==2020)
stations<-filter(stations,first_year<=2013)


####IMPUTE ALL DATASETS USED AFTER FOR ORDINARY KRIGING#######
source('download_GSOD.R')
source('download_GHCND.R')
source('download_solar.R')
##Data files are downloaded from 2 types, since some variables are available in a type of dataset and not in the other (e.g. average relative humidity)

##Data.frame containing names of experiments together with geographical coordinates
daily_weather=read.table(
  'daily_weather_temp_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

all_experiments=unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
                                                 c('2014_ONH1',
                                                   '2014_ONH2',
                                                   '2015_ONH1',
                                                   '2015_ONH2',
                                                   '2016_ONH1',
                                                   '2016_ONH2')]


lapply(all_experiments,
       function(x)
         safeguarding(download_GSOD(
           x,
           radius = 70,
           daily_weather = daily_weather
         )))
lapply(all_experiments,
       function(x)
         safeguarding(download_solar(
           x,
           daily_weather = daily_weather
         )))
lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='TMAX'
         )))

lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='TMIN'
         )))

lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='PRCP'
         )))

lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='PSUN'
         )))

lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='EVAP'
         )))

lapply(all_experiments,
       function(x)
         safeguarding(download_GHCND(
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='PSUN'
         )))
