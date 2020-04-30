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


#library(countyweather)


library(dplyr)
#library(plyr)
library(lubridate)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")


######IMPUTE temp MISSING VALUES #####
daily_weather=read.table(
  'daily_weather_temp_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)


daily_weather$stationID_NOAA=NA
daily_weather$dist=NA
cc=1
par(mfrow=c(2,2))
`%notin%` <- Negate(`%in%`)

source("impute_kriging_withGSOD_server.R")
source('fahrenheit_to_celsius.R')


all_experiments=unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
                                                 c('2014_ONH1',
                                                   '2014_ONH2',
                                                   '2015_ONH1',
                                                   '2015_ONH2',
                                                   '2016_ONH1',
                                                   '2016_ONH2')]




cores <- as.integer(Sys.getenv('SLURM_NTASKS'))

#library(doMPI)
library(doParallel)
#cl <- startMPIcluster(cores)
#registerDoMPI(cl)

results_tmin = mclapply(all_experiments,
                        function(x)
                          safe_impute_function(
                            x,
                            radius = 60,
                            meteo_variable_GSOD = 'MIN',
                            meteo_variable_in_table  = 'TMIN',
                            daily_weather = daily_weather
                          ),mc.cores=cores)

saveRDS(results_tmin,file = 'TMIN/results_tmin.RDS')

results_tmax = mclapply(all_experiments,
                        function(x)
                          safe_impute_function(
                            x,
                            radius = 60,
                            meteo_variable_GSOD = 'MAX',
                            meteo_variable_in_table  = 'TMAX',
                            daily_weather = daily_weather
                          ),mc.cores=cores)


saveRDS(results_tmax,file = 'TMAX/results_tmax.RDS')
#stopCluster(cl)












