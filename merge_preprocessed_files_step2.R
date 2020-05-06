# ------------------------------------------------------------------------------
# Merge all the pre-processed steps (step1) realized on the different meteorological variables
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Run and load all datasets
# ------------------------------------------------------------------------------


#Temperature pre-processed
source('photoperiod_hours_processing_step1.R')
source('incomingsolar_radiation_processing_step1.R')
source('wind_processing_step1.R')
source('rainfall_processing_step1.R')
source('humidity_processing_step1.R')
source('temperature_processing_step1.R')


rm(list = ls())
daylength = read.table(
  'daily_weather_daylength_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
prcp = read.table(
  'daily_weather_prcp_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
humidity = read.table(
  'daily_weather_humidity_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
wind = read.table(
  'daily_weather_wind_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
radiation = read.table(
  'daily_weather_solarrad_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
temp = 
  read.table(
  'daily_weather_temp_processed1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

list1=list(temp,radiation,daylength,humidity,prcp,wind)
multi_full <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  list1
)
