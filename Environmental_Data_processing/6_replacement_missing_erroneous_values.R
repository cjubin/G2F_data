rm(list = ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)
library(mapdata)
library(lubridate)
library(stringr)
library(maps)
library(maptools)
library(xts)
library(spacetime)
library(rgdal)
library(elevatr)
source('fahrenheit_to_celsius.R')

`%notin%` <- Negate(`%in%`)


library(dplyr)
#library(plyr)
library(lubridate)


# ------------------------------------------------------------------------------
# Load dataset with missing values #
# ------------------------------------------------------------------------------
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)
daily_weather = read.table(
  '2_merged_dataset_with_elevation_withONH2017.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
daily_weather=plyr::arrange(daily_weather,Year_Exp,Day.of.Year)

row.names(daily_weather)=NULL
s=nrow(daily_weather)




# ------------------------------------------------------------------------------
# Solar radiation data
# ------------------------------------------------------------------------------

# 1) Reading the NASA files
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/NASA"
)
NASA_FILES <- list()
for (i in list.files(pattern = ".txt")) {
  NASA_FILES[[i]] <- read.table(i, header = T)
}

NASA_FILES <- do.call('rbind', NASA_FILES)
NASA_FILES$day <- lubridate::yday(NASA_FILES$date)

# 2) Correlation between field solar radiation values and NASA files
cor <- list()
na.count <- list()
for (j in unique(NASA_FILES$Year_Exp)) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == j, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == j, 'Date.Harvested'])
  ind <- day_start:day_end
  if (!all(is.na(daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2']))) {
    na.count[[j]] <-
      length(which(is.na(daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2']))) /
      length(daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2'])
    cor[[j]] <-
      cor(NASA_FILES[NASA_FILES$Year_Exp == j &
                       NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN'],
          daily_weather[daily_weather$Year_Exp == j, 'incoming_radiation_MJm2'],
          method = 'pearson',
          use = 'complete.obs')
  }
  else{
    na.count[[j]] <- 1
    cor[[j]] <- 'no field value'
  }
  
}

# 3) Replacement of all field values for which the correlation with NASA values is < 0.6

for (variable in c(names(which(cor < 0.6)), names(which(cor == 'no field value')))) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
  ind <- day_start:day_end
  daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2'] <-
    NASA_FILES[NASA_FILES$Year_Exp == variable &
                 NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN']
}

# 4) Replacement of missing values for other field trials - but only missing days.

for (variable in unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
                                                c(names(which(cor < 0.6)), names(which(cor == 'no field value')))]) {
  day_start <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
  day_end <-
    unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
  ind <- day_start:day_end
  ind_missing_days <-
    which(is.na(daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2']))
  daily_weather[daily_weather$Year_Exp == variable, 'incoming_radiation_MJm2'][ind_missing_days] <-
    NASA_FILES[NASA_FILES$Year_Exp == variable &
                 NASA_FILES$day %in% ind, 'ALLSKY_SFC_SW_DWN'][ind_missing_days]
}

# 5) Generate data.frame with corr

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING')
cor = as.data.frame(cbind(names(cor),unlist(na.count),unlist(cor)))
colnames(cor) <-
  c('Year_Exp',
    'flagged.NA.percent.growing.season.field.data',
    'pearson.corr'
  )
write.table(
  as.data.frame(cor),
  'cor.daily.solar.radiation.fieldstation.vs.NASA.txt',
  col.names = T,
  row.names = F,
  quote = F,
  sep = '\t'
)



# ------------------------------------------------------------------------------
# Metorological variables: 'TMAX', 'TMIN'
# ------------------------------------------------------------------------------

for (meteo_variable in c('TMAX', 'TMIN')) {
  # 1) Reading the imputed data using kriging and constructing summarized dataset from kriging results before imputation
  
  meteo_in_table = meteo_variable
  
  
  setwd(
    paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  list_files1=list.files(pattern = ".RDS")
  
  
  
  summary = as.data.frame(matrix(ncol = 10, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'model.spatiotemporal.variogram',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'pearson.cor.daily.field.vs.official.stations',
    'weather.network.used',
    'radius'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  predictions <- list()
  
  n = 1
  for (i in list_files1) {
    print(i)
    
    v = str_sub(i, 1, nchar(i) - 4)
    summary[n, 1] <- v
    summary[n, 3] <-
      length(which(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable]))) /
      length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
    
    if (length(readRDS(i))==4){
      summary[n, 8] <- readRDS(i)[[2]]
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, 6] <-
        max(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      summary[n, 7] <-
        min(unlist(lapply(readRDS(i)[[4]], mean)), na.rm = t)
      
      summary[n, 5] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[4]], mean
        ))))
    }
    
    if (length(readRDS(i))==3){
      summary[n, 8] <- 'no field data'
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, 6] <-
        max(unlist(lapply(readRDS(i)[[2]], mean)), na.rm = t)
      summary[n, 7] <-
        min(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      
      summary[n, 5] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[3]], mean
        ))))
    }
    
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    if (length(list.files(pattern = pat, recursive = TRUE))==0&length(list.files(pattern = pat2, recursive = TRUE))>0){ #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 9] <- 'WeatherCAN'
      s=list.files(pattern = pat2,recursive = T)
      summary[n, 4] <- length(unique(read.table(s, header = T,sep = '\t')[, 'station_name']))
      
    }
    
    else{
      s=list.files(pattern = pat, recursive = TRUE)
      summary[n, 9] <- 'GHCND'
      summary[n, 4] <- length(unique(read.table(s, header = T,sep = '\t')[, 'id']))}
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
        meteo_variable,
        '/',
        sep = ""
      )
    )
    n <- n + 1
  }
  
  # 3) Replacement of all field values for which the correlation with values imputed with kriging is < 0.6. Same index will be used for TMIN.
  
  for (variable in summary[summary$pearson.cor.daily.field.vs.GHCNDstations <
                           0.6, 'Year_Exp']) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Harvested'])
    ind <- day_start:day_end
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable] <-
      predictions[[variable]][, 2]
  }
  
  # 4) Replacement of missing values for other field trials - but only missing days.
  
  for (variable in summary[summary$pearson.cor.daily.field.vs.GHCNDstations >=
                           0.6, 'Year_Exp']) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
    ind <- day_start:day_end
    ind_missing_days <-
      which(is.na(daily_weather[daily_weather$Year_Exp == variable, meteo_variable]))
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable][ind_missing_days] <-
      predictions[[variable]][, 2][ind_missing_days]
  }
  write.table(summary,file = paste("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",meteo_variable,'_summary.txt',sep = ''),col.names = T,row.names = F,quote = F,sep = '\t')
}

# ------------------------------------------------------------------------------
# Metorological variables: 'AWND'
# ------------------------------------------------------------------------------

for (meteo_variable in c('AWND')) {
  # 1) Reading the imputed data using kriging and constructing summarized dataset from kriging results before imputation
  meteo_in_table = 'MEANWINDSPEED'
  
  
  setwd(
    paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  list_files1=list.files(pattern = ".RDS")
  llist_files2=v = str_sub(list_files1, 1, nchar(i) - 4)
  
  
  summary = as.data.frame(matrix(ncol = 10, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'model.spatiotemporal.variogram',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'pearson.cor.daily.field.vs.weather.stations',
    'weather.network.used',
    'radius'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  predictions <- list()
  
  n = 1
  for (i in c(list_files1)) {
    
    v = str_sub(i, 1, nchar(i) - 4)
    summary[n, 1] <- v
    summary[n, 3] <-
      length(which(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable]))) /
      length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
    
    summary[n, 8] <- readRDS(i)[[2]]
    predictions[[v]] <- readRDS(i)[[1]]
    
    summary[n, 6] <-
      max(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
    summary[n, 7] <-
      min(unlist(lapply(readRDS(i)[[4]], mean)), na.rm = t)
    
    summary[n, 5] <-
      names(which.min(unlist(lapply(
        readRDS(i)[[4]], mean
      ))))
    
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    if (is.null(list.files(pattern = pat, recursive = TRUE)&!is.null(list.files(pattern = pat2, recursive = TRUE))){ #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 9] <- 'WeatherCAN'
      s=list.files(pattern = pat2, recursive = TRUE)
      summary[n, 4] <- length(unique(read.table(s, header = T)[, 'station_name']))
      
    }) 
      
      else{
        s=list.files(pattern = pat, recursive = TRUE)
        summary[n, 9] <- 'GHCND'
        summary[n, 4] <- length(unique(read.table(s, header = T)[, 'id']))}
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
        meteo_variable,
        '/',
        sep = ""
      )
    )
    n <- n + 1
  }
  
  # 3) Replacement of all field values for which the correlation with values imputed with kriging is < 0.6. Same index will be used for TMIN.
  
  for (variable in summary[summary$pearson.cor.daily.field.vs.GHCNDstations <
                           0.6, 'Year_Exp']) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Harvested'])
    ind <- day_start:day_end
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable] <-
      predictions[[variable]][, 2]
  }
  
  # 4) Replacement of missing values for other field trials - but only missing days.
  
  for (variable in summary[summary$pearson.cor.daily.field.vs.GHCNDstations >=
                           0.6, 'Year_Exp']) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
    ind <- day_start:day_end
    ind_missing_days <-
      which(is.na(daily_weather[daily_weather$Year_Exp == variable, meteo_variable]))
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable][ind_missing_days] <-
      predictions[[variable]][, 2][ind_missing_days]
  }
  write.table(summary,file = paste("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",meteo_variable,'_summary.txt',sep = ''),col.names = T,row.names = F,quote = F,sep = '\t')
}

# ------------------------------------------------------------------------------
# Metorological variables: 'PRCP'
# ------------------------------------------------------------------------------

for (meteo_variable in c('PRCP')) {
  # 1) Reading the imputed data using kriging and constructing summarized dataset from kriging results before imputation
  
  setwd(
    paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  
  summary = as.data.frame(matrix(NA,ncol = 11, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'model.spatiotemporal.variogram',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'difference.in.mm.field.vs.weather.stations',
    'pearson.cor.daily.field.vs.weather.stations',
    'weather.network.used',
    'radius'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  predictions <- list()
  
  n = 1
  for (i in list.files(pattern = ".RDS")) {
    print(i)
    v = str_sub(i, 1, nchar(i) - 4)
    
    if (length(readRDS(i))==3){
      summary[n, 1] <- v
      summary[n, 3] <- 1
      summary[n, 8] <- abs(sum(as.numeric(as.vector(readRDS(i)[[1]][,2])),na.rm = T)-sum(daily_weather[daily_weather$Year_Exp == v, meteo_variable],na.rm = T))
      summary[n, 9] <- 'no field data'
      
      predictions[[v]] <- readRDS(i)[[1]]
      summary[n, 6] <-
        max(unlist(lapply(readRDS(i)[[2]], mean)), na.rm = t)
      summary[n, 7] <-
        min(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      
      summary[n, 5] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[3]], mean
        ))))
      
    }
    
    else{
      
      
      summary[n, 1] <- v
      summary[n, 3] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])) /
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      
      summary[n, 8] <- abs(sum(as.numeric(as.vector(readRDS(i)[[1]][,2])),na.rm = T)-sum(daily_weather[daily_weather$Year_Exp == v, meteo_variable],na.rm = T))
      
      summary[n, 9] <- readRDS(i)[[2]]
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, 6] <-
        max(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      summary[n, 7] <-
        min(unlist(lapply(readRDS(i)[[4]], mean)), na.rm = t)
      
      summary[n, 5] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[4]], mean
        ))))
    }
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    if (length(list.files(pattern = pat, recursive = TRUE))==0&length(list.files(pattern = pat2, recursive = TRUE))>0){ #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 10] <- 'WeatherCAN'
      s=list.files(pattern = pat2,recursive = T)
      summary[n, 4] <- length(unique(read.table(s, header = T,sep = '\t')[, 'station_name']))
      
    }
    
    else{
      s=list.files(pattern = pat, recursive = TRUE)
      summary[n, 10] <- 'GHCND'
      summary[n, 4] <- length(unique(read.table(s, header = T,sep = '\t')[, 'id']))}
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",
        meteo_variable,
        '/',
        sep = ""
      )
    )
    n <- n + 1
  }
  
  # 3) Replacement of all field values for which the abs. difference between total prcp field and imputed values if >100
  complete_removal=unique(summary[summary$difference.in.mm.field.vs.GHCNDstations >
                                    100, 'Year_Exp'],summary[summary$flagged.NA.percent.growing.season.field.data>0.3,'Year_Exp'])
  for (variable in complete_removal) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == meteo_variable, 'Date.Harvested'])
    ind <- day_start:day_end
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable] <-
      predictions[[variable]][, 2]
    
  }
  
  
  
  # 4) Replacement of missing values for other field trials - but only missing days.
  
  for (variable in summary$Year_Exp[summary$Year_Exp%notin%complete_removal]) {
    day_start <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Planted'])
    day_end <-
      unique(daily_weather[daily_weather$Year_Exp == variable, 'Date.Harvested'])
    ind <- day_start:day_end
    ind_missing_days <-
      which(is.na(daily_weather[daily_weather$Year_Exp == variable, meteo_variable]))
    daily_weather[daily_weather$Year_Exp == variable, meteo_variable][ind_missing_days] <-
      predictions[[variable]][, 2][ind_missing_days]
  }
  write.table(summary,file = paste("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/IMPUTED_KRIGING/imputation/",meteo_variable,'_summary.txt',sep = ''),col.names = T,row.names = F,quote = F,sep = '\t')
}

# ------------------------------------------------------------------------------
# Add TMEAN 
# ------------------------------------------------------------------------------
daily_weather$TMEAN=(daily_weather$TMAX+daily_weather$TMIN)/2

# ------------------------------------------------------------------------------
# Add the derivation of ea based on HMEAN, if no value could be obtained with the formula using RHmax and RHmin
# es-ea: Saturation vapor pressure deficit, kPa
# ------------------------------------------------------------------------------
source('vapor_pressure.R')
index=which(is.na(daily_weather$ea))
daily_weather$ea[index]=get.ea.with.rhmean(tmin=daily_weather$TMIN[index],tmax=daily_weather$TMAX[index],rhmean=daily_weather$HMEAN[index])
daily_weather$vpd=daily_weather$es-daily_weather$ea

# ------------------------------------------------------------------------------
# Add the GDD and photothermal time (Photothermal time (PTT) is a product between growing degree-days (GDD) and day length (hours) for each day.)
# ------------------------------------------------------------------------------
daily_weather$GDD=(daily_weather$TMAX+daily_weather$TMIN)/2-10
daily_weather$photothermal_time=daily_weather$daylength*daily_weather$GDD


# ------------------------------------------------------------------------------
# Manually add elevation for stations in Canada
# ------------------------------------------------------------------------------
locations_canada=unique(daily_weather[daily_weather$Year_Exp%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2'),c('Year_Exp','lat','long')])
row.names(locations_canada)=NULL
locations_canada$elevation=c(328,204,328,203,336,204,334,203,203)
daily_weather[daily_weather$Year_Exp%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2'),'elev']=locations_canada[match(daily_weather[daily_weather$Year_Exp%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2'),'Year_Exp'],locations_canada$Year_Exp),'elevation']

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
