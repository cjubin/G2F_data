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

results_tmin = mclapply(all_experiments,
                        function(x)
                          safeguarding(impute_kriging_withGSOD(
                            x,
                            radius = 70,
                            meteo_variable_GSOD = 'MIN',
                            meteo_variable_in_table  = 'TMIN',
                            variable_to_impute = 'MIN',
                            daily_weather = daily_weather
                          )),mc.cores=cores)

#saveRDS(results_tmin,file = 'TMIN/results_tmin.RDS')

results_tmax = mclapply(all_experiments,
                        function(x)
                          safeguarding(impute_kriging_withGSOD(
                            x,
                            radius = 70,
                            meteo_variable_GSOD = 'MAX',
                            meteo_variable_in_table  = 'TMAX',
                            variable_to_impute = 'MAX',
                            daily_weather = daily_weather
                          )),mc.cores=cores)


#saveRDS(results_tmax,file = 'TMAX/results_tmax.RDS')







