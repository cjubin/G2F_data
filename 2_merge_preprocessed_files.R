# ------------------------------------------------------------------------------
# Merge all the pre-processed steps (step1) realized on the different meteorological variables
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Run and load all datasets
# ------------------------------------------------------------------------------


#Temperature pre-processed
source('1_photoperiod_hours_processing.R')
source('1_incomingsolar_radiation_processing.R')
source('1_wind_processing.R')
source('1_rainfall_processing.R')
source('1_humidity_processing.R')
source('1_temperature_processing.R')


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
