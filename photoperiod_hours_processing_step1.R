rm(list = ls())
# ------------------------------------------------------------------------------
#Load packages
# ------------------------------------------------------------------------------
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

`%notin%` <- Negate(`%in%`)


library(dplyr)
library(plyr)
library(lubridate)

# ------------------------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------------------------
weather = read.table(
  'weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

daily_weather = read.table(
  'daily_weather.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)


# ------------------------------------------------------------------------------
# Compute with formula daylight based on latitude and day of year
# ------------------------------------------------------------------------------

source('daylength.R')

daily_weather$daylength=daylength(lat=daily_weather$lat,day_of_year = daily_weather$Day.of.Year)
