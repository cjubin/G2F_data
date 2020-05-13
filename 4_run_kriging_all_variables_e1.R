rm(list = ls())
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

args = commandArgs(trailingOnly = TRUE)
x1 = args[1]


######IMPUTE MISSING VALUES #####


daily_weather = read.table(
  '2_merged_dataset.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)



all_experiments = unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
                                                   c('2014_ONH1',
                                                     '2014_ONH2',
                                                     '2015_ONH1',
                                                     '2015_ONH2',
                                                     '2016_ONH1',
                                                     '2016_ONH2')]



cores <- as.integer(Sys.getenv('SLURM_NTASKS'))
library(doParallel)
source('impute_kriging_withGSOD.R')
source('impute_kriging_withGHCND.R')
source('safeguarding.R')


max <- 5
x <- seq_along(1:length(all_experiments))
d1 <- split(1:length(all_experiments), ceiling(x / max))



mclapply(all_experiments[d1[[x1]]],
         function(x)
           safeguarding(
             impute_kriging_withGHCND(
               x,
               radius = 70,
               meteo_variable_GHCND = 'tmin',
               variable_to_impute = 'TMIN',
               name_in_table = 'TMIN',
               daily_weather = daily_weather
             )
           ), mc.cores = cores)
