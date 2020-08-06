rm(list = ls())


## Load or install packages
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing")
source('fahrenheit_to_celsius.R')

packages = c("xts","data.table","rgdal","readxl","opencage","intrval","dplyr","USAboundaries","ggplot2","rnoaa","stringr","ggmap","revgeo","lubridate")


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
daily_weather = read.table(
  '2_merged_dataset_before_interpolation.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
daily_weather = plyr::arrange(daily_weather, Year_Exp, Day.of.Year)

row.names(daily_weather) = NULL
s = nrow(daily_weather)
set_locations=unique(daily_weather$Year_Exp)


# ------------------------------------------------------------------------------
# Replacement Precipitation
# ------------------------------------------------------------------------------

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation')
prcp_summary=read.table('PRCP_summary.txt',header=T,sep="\t")


for (i in set_locations){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL')
  if(prcp_summary[prcp_summary$Year_Exp==i,'Total_partial_replacement']=='Partial')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$PRCP),'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(prcp_summary[prcp_summary$Year_Exp==i,'Winning_kriging_model']!='nearest_station'){dat$yday=lubridate::yday(dat$dates)
    val<-dat[match(dates_to_replace,dat$yday),'Kriging.data']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'PRCP']}
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$PRCP),'PRCP']<-val
    
  }
  if(prcp_summary[prcp_summary$Year_Exp==i,'Total_partial_replacement']=='Total')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(prcp_summary[prcp_summary$Year_Exp==i,'Winning_kriging_model']!='nearest_station'){dat$yday=lubridate::yday(dat$dates)
    val<-dat[match(dates_to_replace,dat$yday),'Kriging.data']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'PRCP']}
    
    daily_weather[daily_weather$Year_Exp==i,'PRCP']<-val
    
  }
}


# ------------------------------------------------------------------------------
# Replacement Temperature
# ------------------------------------------------------------------------------

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation')
tmax_summary=read.table('TMAX_summary.txt',header=T,sep="\t")
tmin_summary=read.table('TMIN_summary.txt',header=T,sep="\t")

for (i in set_locations){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/TMIN')
  if(tmin_summary[tmin_summary$Year_Exp==i,'Replacement']=='Partial')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$TMIN),'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(tmin_summary[tmin_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){
      dat$yday=lubridate::yday(dat$date)
      val<-dat[match(dates_to_replace,dat$yday),'TMIN']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'TMIN']}
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$TMIN),'TMIN']<-val
    
  }
  if(tmin_summary[tmin_summary$Year_Exp==i,'Replacement']=='Total')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(tmin_summary[tmin_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date)
    val<-dat[match(dates_to_replace,dat$yday),'TMIN']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'TMIN']}
    daily_weather[daily_weather$Year_Exp==i,'TMIN']<-val
    
  }
}



for (i in set_locations){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/TMAX')
  if(tmax_summary[tmax_summary$Year_Exp==i,'Replacement']=='Partial')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$TMAX),'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(tmax_summary[tmax_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date)
    val<-dat[match(dates_to_replace,dat$yday),'TMAX']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'TMAX']}
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$TMAX),'TMAX']<-val
    
  }
  if(tmax_summary[tmax_summary$Year_Exp==i,'Replacement']=='Total')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(tmax_summary[tmax_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date)
    val<-dat[match(dates_to_replace,dat$yday),'TMAX']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'TMAX']}
   
    daily_weather[daily_weather$Year_Exp==i,'TMAX']<-val
    
  }
}

# ------------------------------------------------------------------------------
# Add TMEAN 
# ------------------------------------------------------------------------------
daily_weather$TMEAN=(daily_weather$TMAX+daily_weather$TMIN)/2

# ------------------------------------------------------------------------------
# Add the GDD and photothermal time (Photothermal time (PTT) is a product between growing degree-days (GDD) and day length (hours) for each day.)
# ------------------------------------------------------------------------------
daily_weather$GDD=daily_weather$TMEAN-10
daily_weather$photothermal_time=daily_weather$daylength*daily_weather$GDD


# ------------------------------------------------------------------------------
# Replacement HMEAN
# ------------------------------------------------------------------------------

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation')
HMEAN_summary=read.table('HMEAN_summary.txt',header=T,sep="\t")
HMIN_HMAX_summary=read.table('HMIN_HMAX_summary.txt',header=T,sep="\t")

for (i in set_locations[set_locations%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')]){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN')
  if(HMEAN_summary[HMEAN_summary$Year_Exp==i,'Replacement']=='Partial')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$HMEAN),'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(HMEAN_summary[HMEAN_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date)
    val<-dat[match(dates_to_replace,dat$yday),'HMEAN']}
    else{val<-dat[match(dates_to_replace,dat$date),'HMEAN']}
    
    
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$HMEAN),'HMEAN']<-val
    
  }
  if(HMEAN_summary[HMEAN_summary$Year_Exp==i,'Replacement']=='Total')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(HMEAN_summary[HMEAN_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date)
    val<-dat[match(dates_to_replace,dat$yday),'HMEAN']}
    else{val<-dat[match(dates_to_replace,dat$date),'HMEAN']}
    
    daily_weather[daily_weather$Year_Exp==i,'HMEAN']<-val
    
  }
}


for (i in set_locations[set_locations%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')]){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN')
  if(HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp==i,'Replacement']=='Partial')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$HMIN),'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){
      dat$yday=lubridate::yday(dat$date) 
      val1<-dat[match(dates_to_replace,dat$yday),'rhmin']
      val2<-dat[match(dates_to_replace,dat$yday),'rhmax']
    } else{
      val1<-dat[match(dates_to_replace,dat$yday),'nearest_station_values.rhmin']
      val2<-dat[match(dates_to_replace,dat$yday),'nearest_station_values.rhmax']
    }
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$HMIN),'HMIN']<-val1
    daily_weather[daily_weather$Year_Exp==i&is.na(daily_weather$HMAX),'HMAX']<-val2
    
  }
  if(HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp==i,'Replacement']=='Total')
  {
    dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
    dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(HMIN_HMAX_summary[HMIN_HMAX_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){dat$yday=lubridate::yday(dat$date) 
    val1<-dat[match(dates_to_replace,dat$yday),'rhmin']
    val2<-dat[match(dates_to_replace,dat$yday),'rhmax']
    } else{
      val1<-dat[match(dates_to_replace,dat$yday),'nearest_station_values.rhmin']
      val2<-dat[match(dates_to_replace,dat$yday),'nearest_station_values.rhmax']}
    daily_weather[daily_weather$Year_Exp==i,'HMIN']<-val1
    daily_weather[daily_weather$Year_Exp==i,'HMAX']<-val2
    
  }
}




# ------------------------------------------------------------------------------
# Add the derivation of ea based on HMEAN, if no value could be obtained with the formula using RHmax and RHmin
# es-ea: Saturation vapor pressure deficit, kPa
# ------------------------------------------------------------------------------
source('vapor_pressure.R')
index=which(is.na(daily_weather$ea))
daily_weather$ea[index]=get.ea.with.rhmean(tmin=daily_weather$TMIN[index],tmax=daily_weather$TMAX[index],rhmean=daily_weather$HMEAN[index])
daily_weather$vpd=daily_weather$es-daily_weather$ea

# ------------------------------------------------------------------------------
# Replacement Wind
# ------------------------------------------------------------------------------

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation')
wind_summary=read.table('WDSP_summary.txt',header=T,sep="\t")

for (i in set_locations){
  print(i)
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/AWND_GSOD_WeatherCAN')
  dates_to_replace<- daily_weather[daily_weather$Year_Exp==i,'Day.of.Year']
  dat<-read.table(paste0('selected_model_',i,'.txt'),header = T,sep='\t')
    if(wind_summary[wind_summary$Year_Exp==i,'near_or_interpolated_selected']=='interpolated'){
      val<-dat[match(dates_to_replace,dat$yday),'transformed_u2']}
    
    else{val<-dat[match(dates_to_replace,dat$date),'nearest_station_values']}
    
  
  daily_weather[daily_weather$Year_Exp==i,'MEANWINDSPEED']<-val
    
  
}


# ------------------------------------------------------------------------------
# Manually add irrigation values
# ------------------------------------------------------------------------------

irrigation=read.table('IRRIGATION.txt',sep = '\t',header = T)
irrigation$Day.of.Year=yday(as.Date(irrigation$Date,format='%m/%d/%y'))
merged_df=merge(daily_weather,irrigation[,c('Year_Exp','Day.of.Year','Total_mm')],by=c('Year_Exp','Day.of.Year'),all.x=T)
merged_df$total_prcp=rowSums(merged_df[,c("PRCP", "Total_mm")], na.rm=TRUE)

# ------------------------------------------------------------------------------
# Add the reference evapotranspiration rate from the reference surface (ETo) [Penman potential evapotranspiration (mm/day)].
# ------------------------------------------------------------------------------
source('etp.R')
merged_df$net.sum.solar.radiation<-net.sum.solar.radiation.MJ.m2(day_year=merged_df$Day.of.Year,lat=merged_df$lat,elevation=merged_df$elevation,Rs_measured=merged_df$incoming_radiation_MJm2,Tmax=merged_df$TMAX,Tmin=merged_df$TMIN,vpd=merged_df$vpd,with_measurement=T,only_calculated=F)
merged_df$et0=evapotranspiration(elevation = merged_df$elevation,mean.T = merged_df$TMEAN,Rnet = merged_df$net.sum.solar.radiation,daily.mean.wind.speed.ms.s = merged_df$MEANWINDSPEED,vpd=merged_df$vpd)

# ------------------------------------------------------------------------------
# The crop’s water use can be determined by multiplying the reference ETo by
# a crop coefficient (Kc): values reported by FAO are 0.3, 1.2, 0.3–0.6 for the initial, mid-season and late stage, respectively.
# ------------------------------------------------------------------------------
merged_df$etp_initial=0.3*merged_df$et0
merged_df$etp_midseason=1.2*merged_df$et0
merged_df$etp_late=0.5*merged_df$et0
