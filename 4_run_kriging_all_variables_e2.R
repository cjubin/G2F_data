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

`%notin%` <- Negate(`%in%`)


library(dplyr)
#library(plyr)
library(lubridate)



######IMPUTE MISSING VALUES #####


daily_weather=read.table(
  '2_merged_dataset.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
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



cores <- as.integer(Sys.getenv('SLURM_NTASKS'))
library(doParallel)
source('impute_kriging_withGSOD.R')
source('safeguarding.R')


results_prcp = mclapply(all_experiments,
                        function(x)
                          safeguarding(impute_kriging_withGSOD(
                            x,
                            radius = 70,
                            meteo_variable_GSOD = 'PRCP',
                            meteo_variable_in_table  = 'PRCP',
                            variable_to_impute = 'PRCP',
                            daily_weather = daily_weather
                          )),mc.cores=cores)


results_HMEAN = mclapply(all_experiments,
                         function(x)
                           safeguarding(impute_kriging_withGSOD(
                             x,
                             radius = 70,
                             meteo_variable_GSOD = c('DEWP','TEMP'),
                             meteo_variable_in_table  = c('HMEAN'),
                             variable_to_impute = 'HMEAN',
                             daily_weather = daily_weather
                           )),mc.cores=cores)


#saveRDS(results_tmax,file = 'TMAX/results_tmax.RDS')

results_wind = mclapply(all_experiments,
                        function(x)
                          safeguarding(impute_kriging_withGSOD(
                            x,
                            radius = 70,
                            meteo_variable_GSOD = 'WDSP',
                            meteo_variable_in_table  = 'MEANWINDSPEED',
                            variable_to_impute = 'WDSP',
                            daily_weather = daily_weather
                          )),mc.cores=cores)




