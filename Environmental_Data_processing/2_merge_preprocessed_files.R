# ------------------------------------------------------------------------------
# Merge all the pre-processed steps (step1) realized on the different meteorological variables
# ------------------------------------------------------------------------------
rm(list = ls())
# ------------------------------------------------------------------------------
# Run and load all datasets
# ------------------------------------------------------------------------------
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing")

#Temperature pre-processed

# Load the processed datasets according to each variable (range, persistency tests differing according to the weather variable)
source('1_photoperiod_hours_processing.R')
source('1_incomingsolar_radiation_processing.R')
source('1_wind_processing.R')
source('1_rainfall_processing.R')
source('1_humidity_processing.R')
source('1_temperature_processing.R')


rm(list = ls())

## Load the datasets
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

list1=list(daylength,prcp,humidity,wind,radiation,temp)
multi_full <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  list1
)
tmp=grep('flagged',colnames(multi_full))
multi_full=multi_full[,-tmp]
dat<-multi_full

colnames(dat)[which(colnames(dat)=='sum_rainfall')]<-'PRCP'

# ------------------------------------------------------------------------------
# Add saturation vapor pressure (es, kPa) and actual vapor pressure (ea, kPa)
# saturation vapor pressure deficit es-ea, kPa
# ------------------------------------------------------------------------------
source('vapor_pressure.R')
dat$ea=get.ea(rhmin = dat$HMIN,rhmax = dat$HMAX,tmin=dat$TMIN,tmax=dat$TMAX)
dat$es=get.es(tmin = dat$TMIN,tmax = dat$TMAX)
dat$vpd=dat$es-dat$ea
dat=arrange(dat,Year_Exp,Day.of.Year)

write.table(dat,file='2_merged_dataset_before_interpolation.txt',col.names = T,row.names = F,sep = '\t',quote = F)




  

