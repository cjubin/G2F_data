rm(list = ls())


## Load or install packages
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)
source('fahrenheit_to_celsius.R')

packages = c(
  "xts",
  "data.table",
  "rgdal",
  "readxl",
  "opencage",
  "intrval",
  "dplyr",
  "USAboundaries",
  "ggplot2",
  "rnoaa",
  "stringr",
  "ggmap",
  "revgeo",
  "lubridate"
)


package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

`%notin%` <- Negate(`%in%`)



# ------------------------------------------------------------------------------
# Load daily dataset with missing values
# ------------------------------------------------------------------------------
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)
daily_weather = read.table('2_merged_dataset_before_interpolation.txt',
                           header = T,
                           sep = '\t')
daily_weather = plyr::arrange(daily_weather, Year_Exp, Day.of.Year)

row.names(daily_weather) = NULL
s = nrow(daily_weather)
set_locations = unique(daily_weather$Year_Exp)

# ------------------------------------------------------------------------------
# Replacement Solar Radiation values
# ------------------------------------------------------------------------------

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/SOLAR_RADIATION'
)
daily_weather$solar_radiation_NASA = NA

for (i in set_locations) {
  print(i)
  dates_to_replace <-
    daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
  dat <- read.table(paste0('NASA_', i, '.txt'),
                    header = T,
                    sep = '\t')
  dat$yday = lubridate::yday(dat$date)
  
  val <- dat[match(dates_to_replace, dat$yday), 'ALLSKY_SFC_SW_DWN']
  
  daily_weather[daily_weather$Year_Exp == i, 'solar_radiation_NASA'] <-
    val
  
  
}


##Some residual NA


ds <-
  split(daily_weather, daily_weather$Year_Exp)
for (i in 1:length(ds)) {
  name = names(ds[[i]])
  ind <- which(is.na(ds[[i]])[, "solar_radiation_NASA"])
  if (length(ind) != 0) {
    for (s in ind) {
      ds[[i]][s, "solar_radiation_NASA"] <-
        ds[[i]][s + 1, "solar_radiation_NASA"]
    }
  }
}

daily_weather<-bind_rows(ds, .id = "column_label")
  





print('Total replacement solar irradiance values done')
# ------------------------------------------------------------------------------
# Replacement Precipitation
# ------------------------------------------------------------------------------

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation'
)
prcp_summary = read.table('PRCP_summary.txt', header = T, sep = "\t")


for (i in set_locations) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL'
  )
  if (prcp_summary[prcp_summary$Year_Exp == i, 'Total_partial_replacement'] ==
      'Partial')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i &
                      is.na(daily_weather$PRCP), 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (prcp_summary[prcp_summary$Year_Exp == i, 'Winning_kriging_model'] !=
        'nearest_station') {
      dat$yday = lubridate::yday(dat$dates)
      val <- dat[match(dates_to_replace, dat$yday), 'Kriging.data']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'PRCP']
    }
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$PRCP), 'PRCP'] <- val
    
  }
  if (prcp_summary[prcp_summary$Year_Exp == i, 'Total_partial_replacement'] ==
      'Total')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (prcp_summary[prcp_summary$Year_Exp == i, 'Winning_kriging_model'] !=
        'nearest_station') {
      dat$yday = lubridate::yday(dat$dates)
      val <- dat[match(dates_to_replace, dat$yday), 'Kriging.data']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'PRCP']
    }
    
    daily_weather[daily_weather$Year_Exp == i, 'PRCP'] <- val
    
  }
}

daily_weather[which(daily_weather$PRCP < 0), 'PRCP'] <- 0

print('Replacement precipitation values done')
# ------------------------------------------------------------------------------
# Replacement Temperature
# ------------------------------------------------------------------------------

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation'
)
tmax_summary = read.table('TMAX_summary.txt', header = T, sep = "\t")
tmin_summary = read.table('TMIN_summary.txt', header = T, sep = "\t")

for (i in set_locations) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/TMIN'
  )
  if (tmin_summary[tmin_summary$Year_Exp == i, 'Replacement'] == 'Partial')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i &
                      is.na(daily_weather$TMIN), 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (tmin_summary[tmin_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      dat$yday = lubridate::yday(dat$date)
      val <- dat[match(dates_to_replace, dat$yday), 'TMIN']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'TMIN']
    }
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$TMIN), 'TMIN'] <- val
    
  }
  if (tmin_summary[tmin_summary$Year_Exp == i, 'Replacement'] == 'Total')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (tmin_summary[tmin_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      dat$yday = lubridate::yday(dat$date)
      val <- dat[match(dates_to_replace, dat$yday), 'TMIN']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'TMIN']
    }
    daily_weather[daily_weather$Year_Exp == i, 'TMIN'] <- val
    
  }
}


print('Replacement TMIN values done')


for (i in set_locations) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/TMAX'
  )
  if (tmax_summary[tmax_summary$Year_Exp == i, 'Replacement'] == 'Partial')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i &
                      is.na(daily_weather$TMAX), 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (tmax_summary[tmax_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      dat$yday = lubridate::yday(dat$date)
      val <- dat[match(dates_to_replace, dat$yday), 'TMAX']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'TMAX']
    }
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$TMAX), 'TMAX'] <- val
    
  }
  if (tmax_summary[tmax_summary$Year_Exp == i, 'Replacement'] == 'Total')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (tmax_summary[tmax_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      dat$yday = lubridate::yday(dat$date)
      val <- dat[match(dates_to_replace, dat$yday), 'TMAX']
    }
    
    else{
      val <- dat[match(dates_to_replace, dat$date), 'TMAX']
    }
    
    daily_weather[daily_weather$Year_Exp == i, 'TMAX'] <- val
    
  }
}

print('Replacement TMAX values done')

# ------------------------------------------------------------------------------
# Add TMEAN
# ------------------------------------------------------------------------------
daily_weather$TMEAN = (daily_weather$TMAX + daily_weather$TMIN) / 2

# ------------------------------------------------------------------------------
# Add the GDD and photothermal time (Photothermal time (PTT) is a product between growing degree-days (GDD) and day length (hours) for each day.)
# ------------------------------------------------------------------------------
daily_weather$GDD = daily_weather$TMEAN - 10
daily_weather$photothermal_time = daily_weather$daylength * daily_weather$GDD


# ------------------------------------------------------------------------------
# Replacement HMEAN
# ------------------------------------------------------------------------------

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation'
)
HMEAN_summary = read.table('HMEAN_summary.txt', header = T, sep = "\t")
HMIN_HMAX_summary = read.table('HMIN_HMAX_summary.txt', header = T, sep =
                                 "\t")

for (i in set_locations[set_locations %notin% c(
  '2014_ONH1',
  '2014_ONH2',
  '2015_ONH1',
  '2015_ONH2',
  '2016_ONH1',
  '2016_ONH2',
  '2017_ONH1',
  '2017_ONH2',
  '2018_ONH2'
)]) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN'
  )
  if (HMEAN_summary[HMEAN_summary$Year_Exp == i, 'Replacement'] == 'Partial')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i &
                      is.na(daily_weather$HMEAN), 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (HMEAN_summary[HMEAN_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      val <- dat[match(dates_to_replace, dat$dates), 'var1.pred']
    }
    else{
      val <- dat[match(dates_to_replace, dat$date), 'HMEAN']
    }
    
    
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$HMEAN), 'HMEAN'] <- val
    
  }
  if (HMEAN_summary[HMEAN_summary$Year_Exp == i, 'Replacement'] == 'Total')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (HMEAN_summary[HMEAN_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      val <- dat[match(dates_to_replace, dat$dates), 'var1.pred']
    }
    else{
      val <- dat[match(dates_to_replace, dat$date), 'HMEAN']
    }
    
    daily_weather[daily_weather$Year_Exp == i, 'HMEAN'] <- val
    
  }
}


print('Replacement HMEAN values done')

print('Replacement HMIN HMAX for canadian locations start')

for (i in set_locations[set_locations %in% c(
  '2014_ONH1',
  '2014_ONH2',
  '2015_ONH1',
  '2015_ONH2',
  '2016_ONH1',
  '2016_ONH2',
  '2017_ONH1',
  '2017_ONH2',
  '2018_ONH2'
)]) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN'
  )
  if (HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp == i, 'Replacement'] ==
      'Partial')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i &
                      is.na(daily_weather$HMIN), 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      print('interpolated')
      val1 <- dat[match(dates_to_replace, dat$dates), 'rhmin']
      val2 <- dat[match(dates_to_replace, dat$dates), 'rhmax']
    } else{
      val1 <-
        dat[match(dates_to_replace, dat$yday), 'nearest_station_values.rhmin']
      val2 <-
        dat[match(dates_to_replace, dat$yday), 'nearest_station_values.rhmax']
    }
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$HMIN), 'HMIN'] <- val1
    daily_weather[daily_weather$Year_Exp == i &
                    is.na(daily_weather$HMAX), 'HMAX'] <- val2
    
  }
  if (HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp == i, 'Replacement'] ==
      'Total')
  {
    dates_to_replace <-
      daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
    dat <-
      read.table(paste0('selected_model_', i, '.txt'),
                 header = T,
                 sep = '\t')
    if (HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
        'interpolated') {
      val1 <- dat[match(dates_to_replace, dat$dates), 'rhmin']
      val2 <- dat[match(dates_to_replace, dat$dates), 'rhmax']
    } else{
      val1 <-
        dat[match(dates_to_replace, dat$dates), 'nearest_station_values.rhmin']
      val2 <-
        dat[match(dates_to_replace, dat$dates), 'nearest_station_values.rhmax']
    }
    daily_weather[daily_weather$Year_Exp == i, 'HMIN'] <- val1
    daily_weather[daily_weather$Year_Exp == i, 'HMAX'] <- val2
    
  }
}

print('Replacement HMIN/HMAX values done')


# ------------------------------------------------------------------------------
# Add the derivation of ea based on HMAX and HMIN, if the variables could be obtained from local weather station or for canadian locations.
# Otherwise, approximation using HMEAN, if no value could be obtained with the formula using RHmax and RHmin
# es-ea: Saturation vapor pressure deficit, kPa
# ------------------------------------------------------------------------------
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)
source('vapor_pressure.R')
indx<-which(is.na(daily_weather$HMEAN))
print(unique(daily_weather[indx,'Year_Exp']))
daily_weather$HMEAN[indx]=(daily_weather$HMIN[indx]+daily_weather$HMAX[indx])/2
daily_weather$ea = get.ea(
  rhmin = daily_weather$HMIN,
  rhmax = daily_weather$HMAX,
  tmin = daily_weather$TMIN,
  tmax = daily_weather$TMAX
)
daily_weather$es = get.es(tmin = daily_weather$TMIN, tmax = daily_weather$TMAX)
index = which(is.na(daily_weather$ea))
daily_weather$ea[index] = get.ea.with.rhmean(tmin = daily_weather$TMIN[index],
                                             tmax = daily_weather$TMAX[index],
                                             rhmean = daily_weather$HMEAN[index])
daily_weather$vpd = daily_weather$es - daily_weather$ea
index = which(is.na(daily_weather$vpd))

# ------------------------------------------------------------------------------
# Replacement Wind
# ------------------------------------------------------------------------------

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation'
)
wind_summary = read.table('WDSP_summary.txt', header = T, sep = "\t")
daily_weather$wind.speed.m.s = NA

for (i in set_locations) {
  print(i)
  setwd(
    '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/AWND_GSOD_WeatherCAN'
  )
  dates_to_replace <-
    daily_weather[daily_weather$Year_Exp == i, 'Day.of.Year']
  dat <-
    read.table(paste0('selected_model_', i, '.txt'),
               header = T,
               sep = '\t')
  if (wind_summary[wind_summary$Year_Exp == i, 'near_or_interpolated_selected'] ==
      'interpolated') {
    val <- dat[match(dates_to_replace, dat$dates), 'transformed_u2']
  }
  
  else{
    val <- dat[match(dates_to_replace, dat$date), 'transformed_u2']
  }
  
  
  daily_weather[daily_weather$Year_Exp == i, 'wind.speed.m.s'] <- val
  
  
}

# ------------------------------------------------------------------------------
# Manually add environments with same weather data used (same location) but distinguished due to different Plating Date or irrigated/not irrigated.
# ------------------------------------------------------------------------------
daily_weather$Year_Exp=as.character(as.vector(daily_weather$Year_Exp))


add1<-daily_weather[which(daily_weather$Year_Exp=='2016_ILH1'),]
add1$Year_Exp='2016_ILH1.b'
add1$Date.Planted=117
daily_weather[daily_weather$Year_Exp=='2016_ILH1','Year_Exp']<-'2016_ILH1.a'
daily_weather[which(daily_weather$Year_Exp=='2016_ILH1.a'),'Date.Planted']<-127

add2<-daily_weather[which(daily_weather$Year_Exp=='2017_TXH1'),]
add2$Year_Exp='2017_TXH1-Dry'
add2$Date.Planted=62
add2$Date.Harvested=206
add3<-daily_weather[which(daily_weather$Year_Exp=='2017_TXH1'),]
add3$Year_Exp='2017_TXH1-Early'
add3$Date.Planted=62
add3$Date.Harvested=212

daily_weather[daily_weather$Year_Exp=='2017_TXH1','Year_Exp']<-'2017_TXH1-Late'
daily_weather[which(daily_weather$Year_Exp=='2016_ILH1.a'),'Date.Planted']<-96
daily_weather[which(daily_weather$Year_Exp=='2016_ILH1.a'),'Date.Harvested']<-222


add4<-daily_weather[which(daily_weather$Year_Exp=='2018_KSH1'),]
add4$Year_Exp='2018_KSH1.drought'
daily_weather[daily_weather$Year_Exp=='2018_KSH1','Year_Exp']<-'2018_KSH1.irrigated'

daily_weather=rbind(daily_weather,add1,add2,add3,add4)



daily_weather[daily_weather$Year_Exp=='2018_TXH1- Dry','Year_Exp']<-'2018_TXH1-Dry'
daily_weather[daily_weather$Year_Exp=='2018_TXH1- Early','Year_Exp']<-'2018_TXH1-Early'
daily_weather[daily_weather$Year_Exp=='2018_TXH1- Late','Year_Exp']<-'2018_TXH1-Late'

daily_weather = plyr::arrange(daily_weather, Year_Exp, Day.of.Year)
# ------------------------------------------------------------------------------
# Manually add irrigation values
# ------------------------------------------------------------------------------
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)

irrigation = read.table('IRRIGATION.txt', sep = '\t', header = T)
irrigation$Day.of.Year = yday(as.Date(irrigation$Date, format = '%m/%d/%y'))
merged_df = merge(
  daily_weather,
  unique(irrigation[, c('Year_Exp', 'Day.of.Year', 'Total_mm')]),
  by = c('Year_Exp', 'Day.of.Year'),
  all.x = T,
  all.y = F
)
daily_weather$irrigation_mm = merged_df$Total_mm
daily_weather$PRCP2 = rowSums(merged_df[, c("PRCP", "Total_mm")], na.rm =
                                TRUE)

# ------------------------------------------------------------------------------
# Add the reference evapotranspiration rate from the reference surface (ETo) [Penman potential evapotranspiration (mm/day)].
# ------------------------------------------------------------------------------
source('etp.R')

daily_weather$net.sum.solar.radiation <-
  net.sum.solar.radiation.MJ.m2(
    day_year = daily_weather$Day.of.Year,
    lat = daily_weather$lat,
    elevation = daily_weather$elev,
    Rs_measured = daily_weather$solar_radiation_NASA,
    Tmax = daily_weather$TMAX,
    Tmin = daily_weather$TMIN,
    vpd = daily_weather$vpd,
    with_measurement = T,
    only_calculated = F
  )
daily_weather$et0 = evapotranspiration(
  elevation = daily_weather$elev,
  mean.T = daily_weather$TMEAN,
  Rnet = daily_weather$net.sum.solar.radiation,
  daily.mean.wind.speed.ms.s = daily_weather$wind.speed.m.s,
  vpd = daily_weather$vpd
)

# ------------------------------------------------------------------------------
# The crop’s water use can be determined by multiplying the reference ETo by
# a crop coefficient (Kc): values reported by FAO are 0.3, 1.2, 0.3–0.6 for the initial, mid-season and late stage, respectively.
# ------------------------------------------------------------------------------
daily_weather$etp_initial = 0.3 * daily_weather$et0
daily_weather$etp_midseason = 1.2 * daily_weather$et0
daily_weather$etp_late = 0.5 * daily_weather$et0

daily_weather$daily_diffPrec_ETP.V=daily_weather$PRCP2-daily_weather$etp_initial
daily_weather$daily_diffPrec_ETP.F=daily_weather$PRCP2-daily_weather$etp_midseason
daily_weather$daily_diffPrec_ETP.G=daily_weather$PRCP2-daily_weather$etp_late
daily_weather$daily_diffPrec_ET0=daily_weather$PRCP2-daily_weather$et0

print(length(which(is.na(daily_weather$et0))))


write.table(
  daily_weather,
  'replaced_daily_weather.txt',
  col.names = T,
  row.names = F,
  sep = '\t',
  quote = F
)
