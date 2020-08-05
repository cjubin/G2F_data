rm(list = ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)
library(mapdata)
library(lubridate)
library(Metrics)
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
  '2_merged_dataset_with_elevation.txt',
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



# ------------------------------------------------------------------------------
# Include 2017_ONH1 and 2017_ONH2 (no weather stations available)
# ------------------------------------------------------------------------------

canada_2017 = read.table('canada_2017.txt', header = T)
canada_2017$Year_Exp = as.character(canada_2017$Year_Exp)
m = length(canada_2017[1, 5]:canada_2017[1, 6]) + length(canada_2017[2, 5]:canada_2017[2, 6])
daily_weather[nrow(daily_weather) + m, ] <- NA
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Day.of.Year'] <-
  canada_2017[1, 5]:canada_2017[1, 6]
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Year_Exp'] <-
  canada_2017[1, 1]
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Year'] <-
  2017
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Field.Location'] <-
  'ONH1'
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'lat'] <-
  canada_2017[1, 4]
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'long'] <-
  canada_2017[1, 3]
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Date.Planted'] <-
  canada_2017[1, 5]
daily_weather[(s + 1):(s + length(canada_2017[1, 5]:canada_2017[1, 6])), 'Date.Harvested'] <-
  canada_2017[1, 6]

daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Day.of.Year'] <-
  canada_2017[2, 5]:canada_2017[2, 6]
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Year_Exp'] <-
  canada_2017[2, 1]
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Year'] <-
  2017
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Field.Location'] <-
  'ONH2'
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'lat'] <-
  canada_2017[2, 4]
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'long'] <-
  canada_2017[2, 3]
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Date.Planted'] <-
  canada_2017[2, 5]
daily_weather[(s + length(canada_2017[1, 5]:canada_2017[1, 6]) + 1):nrow(daily_weather), 'Date.Harvested'] <-
  canada_2017[2, 6]

daily_weather = plyr::arrange(daily_weather, Year_Exp, Day.of.Year)


write.table(daily_weather,'2_merged_dataset_with_elevation_withONH2017.txt',col.names = T,row.names = F,sep = '\t',quote = F)


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
NASA_FILES[NASA_FILES$ALLSKY_SFC_SW_DWN<0,'ALLSKY_SFC_SW_DWN']<-NA
NASA_FILES = plyr::arrange(NASA_FILES, Year_Exp, day)


# 2) Correlation between field solar radiation values and NASA files
cor <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
names(cor)<-unique(NASA_FILES$Year_Exp)
na.count <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
names(na.count)<-unique(NASA_FILES$Year_Exp)
total_field_station <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
total_NASA_data <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))


for (j in 1:length(unique(NASA_FILES$Year_Exp))) {
  var = unique(NASA_FILES$Year_Exp)[j]
  sub_data=NASA_FILES[NASA_FILES$Year_Exp == var, ]
  matched_values = sub_data[match(daily_weather[daily_weather$Year_Exp ==
                                                  var, 'Day.of.Year'], sub_data$day), 1]
  tryCatch({
    cor[j] <- cor(matched_values,
                  daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'],
                  method = 'pearson',
                  use = 'complete.obs')
  }, error = function(e) 
    return(as.numeric(NA)))
  
  na.count[j] <-
    length(which(is.na(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2']))) /
    length(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'])
  
  total_field_station[j] <- sum(as.numeric(as.vector(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'])),na.rm=T)
  total_NASA_data[j] <- sum(matched_values,na.rm = T)
  
}

# 3) Generate data.frame with corr

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/summary_NASA'
)
cor_table = as.data.frame(cbind(names(cor), na.count, cor,total_field_station,total_NASA_data))
colnames(cor_table) <-
  c('Year_Exp',
    'flagged.NA.percent.growing.season.field.data',
    'pearson.corr','total_field_station_(MJ/m2/day)','total_NASA_data_(MJ/m2/day)')

cor_table$total_NASA_data<-as.numeric(as.vector(cor_table$total_NASA_data))
cor_table$total_field_station<-as.numeric(as.vector(cor_table$total_field_station))

cor_table$difference <-
  cor_table$total_field_station - cor_table$total_NASA_data
cor_table$difference<-as.numeric(as.vector(cor_table$difference))
cor_table$Replacement<-NA
cor_table[abs(cor_table$difference)>300,'Replacement']<-'Total'
cor_table[abs(cor_table$difference)< 0.6,'Replacement']<-'Total'
cor_table[cor_table$Replacement%notin%'Total','Replacement']<-'Partial'




write.table(
  as.data.frame(cor_table),
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
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  list_files1 = list.files(pattern = ".RDS")
  
  
  
  summary = as.data.frame(matrix(ncol = 12, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'model.spatiotemporal.variogram',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'Pearson.correlation.daily.interpolated.with.field.obs',
    'RMSE.interpolated.with.field.obs',
    'weather.network.used',
    'interpolation_method'
   )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  summary$interpolation_method = 'kriging'
  predictions <- list()
  
  n = 1
  for (i in list_files1) {
    print(i)
    
    v = str_sub(i, 1, nchar(i) - 4)
    summary[n, 'Year_Exp'] <- v
    summary[n, 'flagged.NA.percent.growing.season.field.data'] <-
      length(which(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable]))) /
      length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
    
   
    if (length(which(!is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])))>length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])/2){
      summary[n, 'Pearson.correlation.daily.interpolated.with.field.obs'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 2])),
                           daily_weather[daily_weather$Year_Exp == v, meteo_variable],
                           method = 'pearson',
                           use = 'complete.obs')
      summary[n,'RMSE.interpolated.with.field.obs'] <- sqrt( mean( (as.numeric(as.vector(readRDS(i)[[1]][, 2]))- daily_weather[daily_weather$Year_Exp == v, meteo_variable])^2 , na.rm = TRUE ) )
        
                                                                                  }
    else{ summary[n, 'Pearson.correlation.daily.interpolated.with.field.obs'] <- NA
    summary[n,'RMSE.interpolated.with.field.obs'] <-NA}
      
      
      summary[n, 'pearson.cor.5fold.CV.model'] <-
        max(unlist(lapply(readRDS(i)[['5f.cv.kriging.cor']], mean)), na.rm = t)
      summary[n, 'rmse.5fold.CV.model'] <-
        min(unlist(lapply(readRDS(i)[['5f.cv.kriging.rmse']], mean)), na.rm = t)
      
      summary[n, 'model.spatiotemporal.variogram'] <-
        names(which.min(unlist(lapply(
          readRDS(i)[['5f.cv.kriging.rmse']], mean
        ))))
    
    
    values_interpolated=as.numeric(as.vector(readRDS(i)[[1]][, 2]))
    dates_v=as.vector(readRDS(i)[[1]][, 3])
    write.table(as.data.frame(cbind(v,values_interpolated,dates_v)),file=paste(meteo_variable,'_selected_model_',v,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    
      
    
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    
    if (length(list.files(pattern = pat, recursive = TRUE)) == 0 &
        length(list.files(pattern = pat2, recursive = TRUE)) > 0) {
      #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 'weather.network.used'] <- 'WeatherCAN'
      s = list.files(pattern = pat2, recursive = T)
      summary[n, 'nb.weather.stations.used'] <-
        length(unique(read.table(
          s, header = T, sep = '\t'
        )[, 'station_name']))
      
    }
    
    else{
      s = list.files(pattern = pat, recursive = TRUE)
      summary[n, 'weather.network.used'] <- 'GHCND'
      summary[n, 'nb.weather.stations.used'] <-
        length(unique(read.table(
          s, header = T, sep = '\t'
        )[, 'id']))
    }
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
        meteo_variable,
        '/',
        sep = ""
      )
    )
    n <- n + 1
  }
  
  summary$Replacement<-NA
  summary[which(summary$Pearson.correlation.daily.interpolated.with.field.obs<0.65),'Replacement']<-'Total'
  summary[which(summary$flagged.NA.percent.growing.season.field.data>0.5),'Replacement']<-'Total'
  summary[summary$Replacement%notin%'Total','Replacement']<-'Partial'
  
  
  write.table(
    summary,
    file = paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '_summary.txt',
      sep = ''
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}


# ------------------------------------------------------------------------------
# Metorological variables: 'PRCP'
# ------------------------------------------------------------------------------

for (meteo_variable in c('PRCP')) {
  
  # 1) Reading the imputed data using kriging and constructing summarized dataset from kriging results before imputation
  
  setwd(
    paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  
  summary = as.data.frame(matrix(NA, ncol = 19, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'growing_season_length_(days)',
    'model.spatiotemporal.variogram',
    'best_model',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'total_field_station',
    'total_interpolated',
    'diff_interpolated_fieldstation',
    'pearson.cor.daily.field.vs.interpolated',
    'total_1nearest_station',
    'weather.network.used',
    'interpolation_method',
    'distance_1nearest_station_km',
    'name_1nearest_station'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  summary$interpolation_method = 'kriging'
  
  
  predictions <- list()
  
  n = 1
  for (i in list.files(pattern = ".RDS")) {
    print(i)
    v = str_sub(i, 1, nchar(i) - 4)
    
    summary[n, "best_model"]=winning_model[winning_model$Year_Exp==v,'Winning_model']
    
    if (length(readRDS(i)) == 5) {
      summary[n, 'weather.network.used'] <- 'WeatherCAN'
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])) /
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      summary[n, "total_field_station"] <-
        sum(daily_weather[daily_weather$Year_Exp == v, meteo_variable], na.rm = T)
      summary[n, "total_interpolated"] <-
        sum(as.numeric(as.vector(readRDS(i)[[1]][, 2])), na.rm = T)
      
      
      summary[n, "diff_interpolated_fieldstation"] <-
        summary[n, "total_interpolated"]- summary[n, "total_field_station"] 
      
      
      summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 2])),
            daily_weather[daily_weather$Year_Exp == v, meteo_variable],
            method = 'pearson',
            use = 'complete.obs')
      
      
      
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, "pearson.cor.5fold.CV.model"] <-
        max(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      summary[n, "rmse.5fold.CV.model"] <-
        min(unlist(lapply(readRDS(i)[[4]], mean)), na.rm = t)
      
      summary[n,  "model.spatiotemporal.variogram"] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[4]], mean
        ))))
      
      
    } else if (length(readRDS(i)) == 3) {
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])) /
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      summary[n, "total_field_station"] <- NA
      summary[n, "total_interpolated"] <-
        sum(as.numeric(as.vector(readRDS(i)[[1]][, 2])), na.rm = T)
      
      
      summary[n, "diff_interpolated_fieldstation"] <- NA
      
      summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        NA
      
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, "pearson.cor.5fold.CV.model"] <-
        max(unlist(lapply(readRDS(i)[[2]], mean)), na.rm = t)
      summary[n, "rmse.5fold.CV.model"] <-
        min(unlist(lapply(readRDS(i)[[3]], mean)), na.rm = t)
      
      summary[n,  "model.spatiotemporal.variogram"] <-
        names(which.min(unlist(lapply(
          readRDS(i)[[3]], mean
        ))))
      
      
      
    } else{
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])) /
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
      summary[n, "total_field_station"] <-
        sum(daily_weather[daily_weather$Year_Exp == v, meteo_variable], na.rm = T)
      summary[n, "total_interpolated"] <-
        sum(as.numeric(as.vector(readRDS(i)[[1]][, 2])), na.rm = T)
      
      
      summary[n, "diff_interpolated_fieldstation"] <-
         summary[n, "total_interpolated"] - summary[n, "total_field_station"]
      
      
      summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 2])),
            daily_weather[daily_weather$Year_Exp == v, meteo_variable],
            method = 'pearson',
            use = 'complete.obs')
      
      
      predictions[[v]] <- readRDS(i)[[1]]
      
      summary[n, "pearson.cor.5fold.CV.model"] <-
        max(unlist(lapply(readRDS(i)[['5f.cv.kriging.cor']], mean)), na.rm = t)
      summary[n, "rmse.5fold.CV.model"] <-
        min(unlist(lapply(readRDS(i)[['5f.cv.kriging.rmse']], mean)), na.rm = t)
      
      summary[n,  "model.spatiotemporal.variogram"] <-
        names(which.min(unlist(lapply(
          readRDS(i)[['5f.cv.kriging.rmse']], mean
        ))))
    }
    
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    if (length(list.files(pattern = pat, recursive = TRUE)) == 0 &
        length(list.files(pattern = pat2, recursive = TRUE)) > 0) {
      print (paste(v, 'should be a canadian station.'))
      #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 'weather.network.used'] <- 'WeatherCAN'
      s = list.files(pattern = pat2, recursive = T)
      summary[n, "nb.weather.stations.used"] <-
        length(unique(read.table(s, header = T, sep = '\t')[, 'station_name']))
      
    }
    
    else{
      s = list.files(pattern = pat, recursive = TRUE)
      summary[n, 'weather.network.used'] <- 'GHCND'
      summary[n, "nb.weather.stations.used"] <-
        length(unique(read.table(s, header = T, sep = '\t')[, 'id']))
    }
    
    
    setwd(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_1station"
    )
    
    summary[n, "total_1nearest_station"] <-
      sum(
        as.numeric(as.vector(readRDS(
          paste0(v, '_1station.RDS')
        )[[1]][, 2])),na.rm=T)
    
    summary[n, "distance_1nearest_station_km"] <-
      readRDS(
          paste0(v, '_1station.RDS')
        )[[5]]
    
    summary[n, "name_1nearest_station"] <-
      readRDS(
        paste0(v, '_1station.RDS')
      )[[4]]
    
    
        
        
        
        
        
        setwd(
          paste(
            "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
            meteo_variable,
            '/',
            sep = ""
          )
        )
        n <- n + 1
  }
  
  summary$diff_interpolated_neareststation=summary$total_interpolated-summary$total_1nearest_station
  
  summary=summary[,c(1:13,15,16,14,17:19)]
  
  write.table(
    summary,
    file = paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '_summary_with_elevation_allstations.txt',
      sep = ''
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}



# ------------------------------------------------------------------------------
# Metorological variables: 'WDSP'
# ------------------------------------------------------------------------------

for (meteo_variable in c('WDSP')) {
  
  # 1) Reading the imputed data using IDW and constructing summarized dataset from IDW results before imputation
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/AWND_GSOD_WeatherCAN")
  
  summary = as.data.frame(matrix(NA, ncol = 12, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'growing_season_length_(days)',
    'pearson.cor.daily.field.vs.interpolated',
    'weather.network.used',
    'interpolation_method',
    'distance_1nearest_station_km',
    'name_1nearest_station',
    'pearson.cor.interpolated.vs.nearby.stations'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  
  
  predictions <- list()
  
  n = 1
  for (i in list.files(pattern = ".RDS")) {
    print(i)
    v = str_sub(i, 1,9)
    if (length(grep('IDW',i))==1){
      
      summary[n, 'interpolation_method'] <- 'inverse_distance_weighting'
      
      if (v%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')) {
      summary[n, 'weather.network.used'] <- 'WeatherCAN'
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])) /
        length(daily_weather[daily_weather$Year_Exp == v,'MEANWINDSPEED'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])
      
      if (length(which(is.na(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])))>summary[n, "growing_season_length_(days)"]/2){
        summary[n, "pearson.cor.daily.field.vs.interpolated"] <- NA
      }
      
        else{summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 'transformed_u2'])),
            daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'],
            method = 'pearson',
            use = 'complete.obs')}
      
      summary[n,'nb.weather.stations.used'] <- unique(readRDS(i)[[1]][,'nb_stations_used'])
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(i)[['min_dist']]
      
      summary[n,'name_1nearest_station'] <- readRDS(i)[['min_dist_name']]
      
      summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'transformed_u2'])),
                                                                      as.numeric(as.vector(readRDS(i)[[1]][, 'nearest_station_values'])),
                                                                      method = 'pearson',
                                                                      use = 'complete.obs')
      
    }  else{
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])) /
        length(daily_weather[daily_weather$Year_Exp == v,'MEANWINDSPEED'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])
      
      
      if (length(which(is.na(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])))>summary[n, "growing_season_length_(days)"]/2){
        summary[n, "pearson.cor.daily.field.vs.interpolated"] <- NA
      }
      
      else{summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 'transformed_u2'])),
            daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'],
            method = 'pearson',
            use = 'complete.obs')}
      
      summary[n,'nb.weather.stations.used'] <- unique(readRDS(i)[[1]][,'nb_stations_used'])
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(i)[[5]]
      
      summary[n,'name_1nearest_station'] <- readRDS(i)[[4]]
      
      summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'transformed_u2'])),
                                                                      as.numeric(as.vector(readRDS(i)[[1]][, 'nearest_station_values'])),
                                                                      method = 'pearson',
                                                                      use = 'complete.obs')
      
      
    
    }
    } else {
      summary[n, 'interpolation_method'] <- 'only_1_station_available'
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])) /
        length(daily_weather[daily_weather$Year_Exp == v,'MEANWINDSPEED'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'])
      
      
      summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 'transformed_u2'])),
            daily_weather[daily_weather$Year_Exp == v, 'MEANWINDSPEED'],
            method = 'pearson',
            use = 'complete.obs')
      
      summary[n,'nb.weather.stations.used'] <- 1
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(i)[[4]]
      
      if (readRDS(i)[[3]]=='72562194063'){summary[n,'name_1nearest_station'] <- 'OGALLALA SEARLE FIELD AIRPORT, NE US'}
      if (readRDS(i)[[3]]=='72562024023'){summary[n,'name_1nearest_station'] <- 'NORTH PLATTE REGIONAL AIRPORT, NE US'}
      
      
      
      summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- NA
      
      
    }
    
    
    
    n <- n + 1
   
      
}
    
    
    
  
    
    
    
    
   
  
  
  
  
  write.table(
    summary,
    file = paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '_summary.txt',
      sep = ''
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}




# ------------------------------------------------------------------------------
# Metorological variables: 'HMEAN','HMAX','HMIN'
# ------------------------------------------------------------------------------

for (meteo_variable in c('HMEAN')) {
  
  # 1) Reading the imputed data using IDW and constructing summarized dataset from IDW results before imputation
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN")
  
  summary = as.data.frame(matrix(NA, ncol = 16, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'weather.network.used',
    'growing_season_length_(days)',
    'RHMEAN.pearson.cor.daily.field.vs.interpolated',
    'RHMAX.pearson.cor.daily.field.vs.interpolated',
    'RHMIN.pearson.cor.daily.field.vs.interpolated',
    'interpolation_method',
    'distance_1nearest_station_km',
    'name_1nearest_station',
    'RHMEAN.pearson.cor.interpolated.vs.nearby.stations',
    'RHMAX.pearson.cor.interpolated.vs.nearby.stations',
    'RHMIN.pearson.cor.interpolated.vs.nearby.stations'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  summary$meteo_variable='HMEAN;HMIN;HMAX'
  
  
 
  
  n = 1
  for (i in list.files(pattern = ".RDS")) {
    print(i)
    v = str_sub(i, 1,9)
    if (length(grep('IDW',i))==1){
      
      summary[n, 'interpolation_method'] <- 'inverse_distance_weighting'
      
      if (v%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')) {
        summary[n, 'weather.network.used'] <- 'WeatherCAN'
        summary[n, "Year_Exp"] <- v
        
        summary[n, "flagged.NA.percent.growing.season.field.data"] <-
          sum(is.na(daily_weather[daily_weather$Year_Exp == v, c('HMAX',"HMIN")])) /
          (nrow(daily_weather[daily_weather$Year_Exp == v,c('HMAX',"HMIN")])*2)
        
        summary[n, "growing_season_length_(days)"] <-
          nrow(daily_weather[daily_weather$Year_Exp == v,c('HMAX',"HMIN")])
        
      
        
        if (v%in%c('2017_ONH1','2017_ONH2')) {summary[n,  'RHMEAN.pearson.cor.daily.field.vs.interpolated'] <- NA
        summary[n,  'RHMAX.pearson.cor.daily.field.vs.interpolated'] <- NA
        summary[n,  'RHMIN.pearson.cor.daily.field.vs.interpolated'] <- NA
        }

        else{summary[n,  'RHMEAN.pearson.cor.daily.field.vs.interpolated'] <- NA
        summary[n,  'RHMAX.pearson.cor.daily.field.vs.interpolated'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'rhmax'])),
                                                                            daily_weather[daily_weather$Year_Exp == v, 'HMAX'],
                                                                                method = 'pearson',
                                                                                use = 'complete.obs')
        summary[n,  'RHMIN.pearson.cor.daily.field.vs.interpolated'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'rhmin'])),
                                                                             daily_weather[daily_weather$Year_Exp == v, 'HMIN'],
                                                                                 method = 'pearson',
                                                                                 use = 'complete.obs')}
        
        summary[n,  'RHMEAN.pearson.cor.interpolated.vs.nearby.stations'] <- NA
        summary[n,  'RHMAX.pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'rhmax'])),
                                                                                as.numeric(as.vector(readRDS(i)[[1]][, 8][,1])),
                                                                                method = 'pearson',
                                                                                use = 'complete.obs')
        summary[n,  'RHMIN.pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'rhmin'])),
                                                                                as.numeric(as.vector(readRDS(i)[[1]][, 8][,2])),
                                                                                method = 'pearson',
                                                                                use = 'complete.obs')
        summary[n,'nb.weather.stations.used'] <- unique(readRDS(i)[[1]][,'nb_stations_used'])
        
        summary[n,'distance_1nearest_station_km'] <- readRDS(i)[['min_dist']]
        
        summary[n,'name_1nearest_station'] <- readRDS(i)[['min_dist_name']]
        
        
        }else{
        
        summary[n, 'weather.network.used'] <- 'GSOD'
        summary[n, "Year_Exp"] <- v
        summary[n, "flagged.NA.percent.growing.season.field.data"] <-
          sum(is.na(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])) /
          length(daily_weather[daily_weather$Year_Exp == v,'HMEAN'])
        
        summary[n, "growing_season_length_(days)"] <-
          length(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])
          
        if (length(which(is.na(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])))>0.85*length(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])){
          summary[n, 'RHMEAN.pearson.cor.daily.field.vs.interpolated'] <- NA
        } else{summary[n,  'RHMEAN.pearson.cor.daily.field.vs.interpolated'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'var1.pred'])),
                                                                             daily_weather[daily_weather$Year_Exp == v, 'HMEAN'],
                                                                             method = 'pearson',
                                                                             use = 'complete.obs')}
        summary[n,  'RHMAX.pearson.cor.daily.field.vs.interpolated'] <- NA
        summary[n,  'RHMIN.pearson.cor.daily.field.vs.interpolated'] <- NA
          
        
        summary[n,'nb.weather.stations.used'] <- unique(readRDS(i)[[1]][,'nb_stations_used'])
        
        summary[n,'distance_1nearest_station_km'] <- readRDS(i)[['min_dist']]
        
        summary[n,'name_1nearest_station'] <- readRDS(i)[['min_dist_name']]
        
        
        
        summary[n,  'RHMEAN.pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 'var1.pred'])),
                                                                                 as.numeric(as.vector(readRDS(i)[[1]][, 'nearest_station_values'])),
                                                                                 method = 'pearson',
                                                                                 use = 'complete.obs')
        
        summary[n,  'RHMAX.pearson.cor.interpolated.vs.nearby.stations'] <- NA
        summary[n,  'RHMIN.pearson.cor.interpolated.vs.nearby.stations'] <- NA
        
        
      } } else {
      summary[n, 'interpolation_method'] <- 'only_1_station_available'
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- v
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])) /
        length(daily_weather[daily_weather$Year_Exp == v,'HMEAN'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == v, 'HMEAN'])
      
      
      summary[n, "RHMEAN.pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(i)[[1]][, 'HMEAN'])),
            daily_weather[daily_weather$Year_Exp == v, 'HMEAN'],
            method = 'pearson',
            use = 'complete.obs')
      
      summary[n,'nb.weather.stations.used'] <- 1
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(i)[[4]]
      
      if (readRDS(i)[[3]]=='72562194063'){summary[n,'name_1nearest_station'] <- 'OGALLALA SEARLE FIELD AIRPORT, NE US'}
      if (readRDS(i)[[3]]=='72562024023'){summary[n,'name_1nearest_station'] <- 'NORTH PLATTE REGIONAL AIRPORT, NE US'}
      
      
      
      summary[n,  'RHMEAN.pearson.cor.interpolated.vs.nearby.stations'] <- NA
      summary[n,  'RHMAX.pearson.cor.interpolated.vs.nearby.stations'] <- NA
      summary[n,  'RHMIN.pearson.cor.interpolated.vs.nearby.stations'] <- NA
      
      
     
    }
    
    
    n<-n+1
    
    
    
  
 
  }
  write.table(
    summary,
    file = paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '_summary.txt',
      sep = ''
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}




# ------------------------------------------------------------------------------
# Metorological variables: 'PRCP'
# ------------------------------------------------------------------------------



setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/aggregate_prcp')

for (file in list.files()){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/')

daily_weather = read.table(
  '2_merged_dataset_with_elevation_withONH2017.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
`%notin%` <- Negate(`%in%`)
all_experiments = unique(daily_weather$Year_Exp)
setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GHCND/imputation/all_stations_PRCP/more_models/')
american_exp  = unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
                                                 c('2014_ONH1',
                                                   '2014_ONH2',
                                                   '2015_ONH1',
                                                   '2015_ONH2',
                                                   '2016_ONH1',
                                                   '2016_ONH2',
                                                   '2017_ONH1',
                                                   '2017_ONH2',
                                                   '2018_ONH2')]

canadian_exp  = c(
  '2014_ONH1',
  '2014_ONH2',
  '2015_ONH1',
  '2015_ONH2',
  '2016_ONH1',
  '2016_ONH2',
  '2017_ONH1',
  '2017_ONH2',
  '2018_ONH2'
)


winning_model = matrix(NA, nrow = length(all_experiments), ncol = 16)
colnames(winning_model) = c(
  'Year_Exp',
  'Number.days.growing.season',
  'Aggregate_raw_data(mm)',
  'Total_field_after_processing(mm)',
  '%NA.field.season',
  'Winning_model',
  'Total.winning.model(mm)',
  'Total.with.elevation(mm)',
  'RMSE.model.with.elevation',
  'Pearson.Cor.model.with.elevation',
  'Total.no.elevation(mm)',
  'RMSE.model.no.elevation',
  'Pearson.Cor.model.no.elevation',
  'Total.near.station(mm)',
  'Distance.near.station',
  'Name.near.station'
)


winning_model=as.data.frame(winning_model)
for (j in c(2,3,4,6,7,8,9,10,11,12,13,14) ){
  winning_model[,j]<-as.numeric(as.vector(winning_model[,j]))
}
n = 1
for (s in american_exp) {
  if (s=='2017_COH1'){
    n<-n+1 
    next}
  print(s)
  
  winning_model[n, 'Year_Exp'] <- s
  winning_model[n, 'Number.days.growing.season'] <-length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  winning_model[n, 'Aggregate_raw_data(mm)']<-unique(dataset[dataset$Category==s,2])
  winning_model[n, 'Total_field_after_processing(mm)']<-sum(daily_weather[daily_weather$Year_Exp==s,'PRCP'],na.rm = T)
  winning_model[n, '%NA.field.season']<-length(which(is.na(daily_weather[daily_weather$Year_Exp==s,'PRCP'])))/length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  
  
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_1station/all_stations')
  
  
  vector_nearest_station<-as.numeric(as.vector(readRDS(paste(s,'_1station.RDS',sep = ''))[[1]][,2]))
  winning_model[n, 'Total.near.station(mm)']<-sum(as.numeric(as.vector(readRDS(paste(s,'_1station.RDS',sep = ''))[[1]][,2])),na.rm = T)
  winning_model[n,  'Distance.near.station']<-readRDS(paste(s,'_1station.RDS',sep = ''))[['min_dist']]
  winning_model[n, 'Name.near.station']<-readRDS(paste(s,'_1station.RDS',sep = ''))[['min_dist_station_name']]
  
  
  
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GHCND/imputation/all_stations_PRCP/more_models")
  
  elev_model<-readRDS(paste('PRCP_elev',s,'.RDS',sep = ''))[[1]]
  
  elev_model<-cbind(elev_model,vector_nearest_station)
  
  if (!file.exists(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))|!file.exists(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))){
    n <- n + 1
    next
  } else{
    LOOCV<-readRDS(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))
    
    
    LOOCV_min_RMSE_elev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    LOOCV_max_cor_elev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    
    winning_model[n, 'Total.with.elevation(mm)'] <-sum(elev_model[,LOOCV_min_RMSE_elev])
    winning_model[n, 'RMSE.model.with.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    winning_model[n, 'Pearson.Cor.model.with.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    
    no_elev_model<-readRDS(paste('PRCPno_elev',s,'.RDS',sep = ''))[[1]]
    no_elev_model<-cbind(no_elev_model,vector_nearest_station)
    
    
    LOOCV<-readRDS(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))
    names(LOOCV[['kriging_rmse']])[names(LOOCV[['kriging_rmse']])=='summetric_Vgm']<-'summetric_Vgm1'
    names(LOOCV[['kriging_cor']])[names(LOOCV[['kriging_cor']])=='summetric_Vgm']<-'summetric_Vgm1'
    LOOCV_min_RMSE_noelev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    LOOCV_max_cor_noelev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    winning_model[n, 'Total.no.elevation(mm)'] <-sum(no_elev_model[,LOOCV_min_RMSE_noelev])
    winning_model[n, 'RMSE.model.no.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    winning_model[n, 'Pearson.Cor.model.no.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL")
    
    if(winning_model[n,'Distance.near.station']<1.5){
      winning_model[n, 'Winning_model'] <- 'nearest_station'
      winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.near.station(mm)']
      write.table(no_elev_model[,c('predictions.table','vector_nearest_station','as.character(dates)')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    }
    else if (winning_model[n, 'RMSE.model.no.elevation']<winning_model[n, 'RMSE.model.with.elevation']){
      if(max(no_elev_model[,LOOCV_min_RMSE_noelev])<100&var(no_elev_model[,LOOCV_min_RMSE_noelev])>4){
        winning_model[n, 'Winning_model'] <- paste('no_elevation_',LOOCV_min_RMSE_noelev,sep='')
        winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.no.elevation(mm)']
        write.table(no_elev_model[,c('predictions.table',LOOCV_min_RMSE_noelev,'as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }else {winning_model[n, 'Winning_model'] <-'issue_with_model'}
    }else {
      if(max(elev_model[,LOOCV_min_RMSE_elev])<100&var(elev_model[,LOOCV_min_RMSE_elev])>4){
        winning_model[n, 'Winning_model'] <- paste('elevation_',LOOCV_min_RMSE_elev,sep='')
        winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.with.elevation(mm)']
        write.table(elev_model[,c('predictions.table',LOOCV_min_RMSE_elev,'as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      } else {winning_model[n, 'Winning_model'] <-'issue_with_model'}}
    
    
    if (winning_model[n, 'Winning_model']=='issue_with_model'){
      if('summetric_Vgm3'%in%colnames(no_elev_model)){
        winning_model[n, 'Winning_model'] <- paste('elevation_','summetric_Vgm3',sep='')
        winning_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm3'],na.rm=T)
        
        write.table(no_elev_model[,c('predictions.table','summetric_Vgm3','as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }
      if('summetric_Vgm5'%in%colnames(no_elev_model)){
        winning_model[n, 'Winning_model'] <- paste('elevation_','summetric_Vgm5',sep='')
        winning_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm5'],na.rm=T)
        
        write.table(no_elev_model[,c('predictions.table','summetric_Vgm5','as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }
      
    }
    
    
    
    
    
    n <- n + 1
  }
  
  
  
}

########################################
####### CANADIAN LOCATIONS #############
n=105

for (s in canadian_exp) {
  print(s)
  
  winning_model[n, 'Year_Exp'] <- s
  winning_model[n, 'Number.days.growing.season'] <-length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  winning_model[n, 'Aggregate_raw_data(mm)']<-unique(dataset[dataset$Category==s,2])
  winning_model[n, 'Total_field_after_processing(mm)']<-sum(daily_weather[daily_weather$Year_Exp==s,'PRCP'],na.rm = T)
  winning_model[n, '%NA.field.season']<-length(which(is.na(daily_weather[daily_weather$Year_Exp==s,'PRCP'])))/length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  
  
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_1station/all_stations')
  
  
  vector_nearest_station<-as.numeric(as.vector(readRDS(paste(s,'_1station.RDS',sep = ''))[[1]][,2]))
  winning_model[n, 'Total.near.station(mm)']<-sum(as.numeric(as.vector(readRDS(paste(s,'_1station.RDS',sep = ''))[[1]][,2])),na.rm = T)
  winning_model[n,  'Distance.near.station']<-readRDS(paste(s,'_1station.RDS',sep = ''))[['min_dist']]
  winning_model[n, 'Name.near.station']<-readRDS(paste(s,'_1station.RDS',sep = ''))[['min_dist_station_name']]
  
  
  
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GHCND/imputation/all_stations_PRCP/more_models")
  
  elev_model<-readRDS(paste('PRCP_elev',s,'.RDS',sep = ''))[[1]]
  
  elev_model<-cbind(elev_model,vector_nearest_station)
  
  if (!file.exists(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))|!file.exists(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))){
    n <- n + 1
    next
  } else{
    LOOCV<-readRDS(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))
    
    
    LOOCV_min_RMSE_elev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    LOOCV_max_cor_elev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    
    winning_model[n, 'Total.with.elevation(mm)'] <-sum(elev_model[,LOOCV_min_RMSE_elev])
    winning_model[n, 'RMSE.model.with.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    winning_model[n, 'Pearson.Cor.model.with.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    
    no_elev_model<-readRDS(paste('PRCPno_elev',s,'.RDS',sep = ''))[[1]]
    no_elev_model<-cbind(no_elev_model,vector_nearest_station)
    
    
    LOOCV<-readRDS(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))
    names(LOOCV[['kriging_rmse']])[names(LOOCV[['kriging_rmse']])=='summetric_Vgm']<-'summetric_Vgm1'
    names(LOOCV[['kriging_cor']])[names(LOOCV[['kriging_cor']])=='summetric_Vgm']<-'summetric_Vgm1'
    LOOCV_min_RMSE_noelev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    LOOCV_max_cor_noelev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    winning_model[n, 'Total.no.elevation(mm)'] <-sum(no_elev_model[,LOOCV_min_RMSE_noelev])
    winning_model[n, 'RMSE.model.no.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    winning_model[n, 'Pearson.Cor.model.no.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    
    if(winning_model[n,'Distance.near.station']<1.5){
      winning_model[n, 'Winning_model'] <- 'nearest_station'
      winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.near.station(mm)']
      write.table(no_elev_model[,c('predictions.table','vector_nearest_station','as.character(dates)')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    }
    else if (winning_model[n, 'RMSE.model.no.elevation']<winning_model[n, 'RMSE.model.with.elevation']){
      if(max(no_elev_model[,LOOCV_min_RMSE_noelev])<100&var(no_elev_model[,LOOCV_min_RMSE_noelev])>4){
        winning_model[n, 'Winning_model'] <- paste('no_elevation_',LOOCV_min_RMSE_noelev,sep='')
        winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.no.elevation(mm)']
        write.table(no_elev_model[,c('predictions.table',LOOCV_min_RMSE_noelev,'as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }else {winning_model[n, 'Winning_model'] <-'issue_with_model'}
    }else {
      if(max(elev_model[,LOOCV_min_RMSE_noelev])<100&var(elev_model[,LOOCV_min_RMSE_noelev])>4){
        winning_model[n, 'Winning_model'] <- paste('elevation_',LOOCV_min_RMSE_elev,sep='')
        winning_model[n, 'Total.winning.model(mm)']<-winning_model[n, 'Total.with.elevation(mm)']
        write.table(elev_model[,c('predictions.table',LOOCV_min_RMSE_elev,'as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      } else {winning_model[n, 'Winning_model'] <-'issue_with_model'}}
    
    
    
    if (winning_model[n, 'Winning_model']=='issue_with_model'){
      winning_model[n, 'Winning_model'] <- paste('elevation_','summetric_Vgm3',sep='')
      winning_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm3'],na.rm=T)
      
      write.table(no_elev_model[,c('predictions.table','summetric_Vgm3','as.character(dates)','vector_nearest_station')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      
      
    }
    
    
    
    
    
    n <- n + 1
  }
  
  
  
}













#######################################
#######################################


winning_model$`Total.with.elevation(mm)`=as.numeric(as.vector(winning_model$`Total.with.elevation(mm)`))
winning_model$`Total.near.station(mm)`=as.numeric(as.vector(winning_model$`Total.near.station(mm)`))
winning_model$`Total.no.elevation(mm)`=as.numeric(as.vector(winning_model$`Total.no.elevation(mm)`))


winning_model$Difference.selected.model.field.station=winning_model$`Total.winning.model(mm)`-winning_model$`Total_field_after_processing(mm)`
winning_model$Difference.selected.model.nearest.station=winning_model$`Total.winning.model(mm)`-winning_model$`Total.near.station(mm)`

##Total or partial replacement by interpolated/nearest stations based on 1)% NA data from field 2)Difference between field and interpolated or nearest stations

winning_model[which(winning_model$`%NA.field.season`>0.4),'Total_partial_replacement']<-'Total'
winning_model[which(winning_model$`%NA.field.season`<0.4),'Total_partial_replacement']<-'Partial'

winning_model[which(winning_model$Difference.selected.model.field.station>90|winning_model$Difference.selected.model.field.station<(-90)),'Total_partial_replacement']<-'Total'

#######################################
#######################################

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL")
write.table(winning_model,'winning_model_prcp.txt',col.names = T,row.names = F,sep = '\t',quote = F)


