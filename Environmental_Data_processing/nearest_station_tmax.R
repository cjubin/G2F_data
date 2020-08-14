#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{nearest_station_tmax} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable_GHCND Character. GHCND element names (= variable measured at the GHCND station): TMAX, TMIN, PRCP
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested' and @variable_to_impute
#' @param variable_to_impute. Character.
#' @param name_in_table. Character

rm(list = ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(geosptdb)
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
library(weathercan)
source('safeguarding.R')

`%notin%` <- Negate(`%in%`)


library(dplyr)
library(lubridate)

setwd("C:/Users/cathyjubin/Desktop")
daily_weather = read.table(
  'subset_1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")

stations <- ghcnd_stations()
stations <- filter(stations, last_year == 2020)
stations <- filter(stations, first_year <= 2013)



nearest_station_tmax <-
  function(Year_Exp,
           radius = 70,
           meteo_variable_GHCND = 'tmax',
           daily_weather,
           name_in_table = 'TMAX') {
    
    print(Year_Exp)
    #Retrieve information about the experiment
    year = as.numeric(unique(daily_weather[daily_weather$Year_Exp == Year_Exp, 'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp ==
                                                                                                                      Year_Exp, 'Year']))])
    date_start = as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp ==
                                                           Year_Exp, 'Date.Planted'])), origin = paste(year - 1, '12', '31', sep = '-'))
    date_end = as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp ==
                                                         Year_Exp, 'Date.Harvested'])), origin = paste(year - 1, '12', '31', sep = '-'))
    
    dates = seq(as.Date(date_start, tz = "CET"),
                as.Date(date_end, tz = "CET"),
                by = 'days')
    
    longitude = as.numeric(unique(daily_weather[daily_weather$Year_Exp == Year_Exp, 'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp ==
                                                                                                                           Year_Exp, 'long']))])
    latitude = as.numeric(unique(daily_weather[daily_weather$Year_Exp == Year_Exp, 'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp ==
                                                                                                                         Year_Exp, 'lat']))])
    
    #Values recorded for this variable at the field station
    if (!is.null(name_in_table) &
        name_in_table %in% colnames(daily_weather)) {
      field_values = daily_weather[daily_weather$Year_Exp == Year_Exp, name_in_table]
      
    } else{field_values=NA}
    
    
    
    ## Select the data from GHCND stations if they are (case for he majority). Otherwise, take from GSOD data.
    
    if (file.exists(
      paste(
        "C:/Users/cathyjubin/Desktop/wind/ghcnd_files/",
        year,
        '/',
        name_in_table,
        '_',
        Year_Exp,
        '.txt',
        sep = ''
      )
    )) {
      d = read.table(
        file =  paste(
          "C:/Users/cathyjubin/Desktop/wind/ghcnd_files/",
          year,
          '/',
          name_in_table,
          '_',
          Year_Exp,
          '.txt',
          sep = ''
        ),
        header = T,sep = '\t'
      )
      nb_stations = length(unique(d$id))
      stations_close = as.data.frame(rnoaa::meteo_distance(stations, latitude, longitude, radius = radius))
      stations_close = unique(stations_close[stations_close$id %in% unique(d$id),])
      
      
      min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                       min(stations_close$distance)), 'id'])
      min_dist = unique(min(stations_close$distance))
      
      min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                    min(stations_close$distance)), 'name'])
      
      
      nearest_station_values <- d[d$id == min_dist_station,]
      nearest_station_values <-
        nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), 'tmax']
      
      while (length(nearest_station_values) < length(dates) | length(which(is.na(nearest_station_values))) > 8) {
        stations_close = stations_close[-which(stations_close$id == min_dist_station),]
        min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                         min(stations_close$distance)), 'id'])
        min_dist = unique(min(stations_close$distance))
        
        min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                      min(stations_close$distance)), 'name'])
        
        nearest_station_values <- d[d$id == min_dist_station,]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), 'tmax']
        print(sum(nearest_station_values,na.rm = T))
      }
      
      d = d[, c(2, 3, 4, 5)]
      
      colnames(d)[which(colnames(d) == meteo_variable_GHCND)] = 'TMAX'
      colnames(d)[which(colnames(d) == 'date')] = 'dates'
      
      
      d$longitude = as.numeric(as.vector(d$longitude))
      d$latitude = as.numeric(as.vector(d$latitude))
      
      
      d$dates = as.Date(d$dates)
      d$dates = lubridate::yday(d$dates)
      
      
      
      
      
      network = 'GHCND'
      print('Data from GHCND stations prepared')
      
      ########################
      ########################
      
      
      
      
      #d$day<-lubridate::yday(d$date)
      dates = yday(seq(
        as.Date(date_start, tz = "CET"),
        as.Date(date_end, tz = "CET"),
        by = 'days'
      ))
      
      predicted.values <-
        nearest_station_values
      print(predicted.values)
      print(length(predicted.values))
      
      if (length(field_values) > 1 & !all(is.na(field_values))) {
        cors = cor(predicted.values, field_values, use = 'complete.obs')
      } else{
        cors = 'no_field_values'
      }
      
      predictions.table = cbind(Year_Exp, predicted.values, as.character(dates), network,field_values)
      predictions.table = as.data.frame(predictions.table)
      colnames(predictions.table) = c('Year_Exp', 'TMAX', 'date', 'network','field_values')
      
      to_save = list(predictions.table, cors, min_dist_station, min_dist_name,min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist_station_name',
          'min_dist'
        )
      
      saveRDS(to_save, paste('TMAX_nearest/', Year_Exp, '_1station.RDS', sep = ''))
      
      
      
      
    }
    else if (file.exists(
      paste(
        'C:/Users/cathyjubin/Desktop/wind/ghcnd_files/',
        year,
        '/',
        Year_Exp,
        '_weather',
        '.txt',
        sep = ''
      )
    )) {
      d = read.table(
        file =  paste(
          'C:/Users/cathyjubin/Desktop/wind/ghcnd_files/',
          year,
          '/',
          Year_Exp,
          '_weather',
          '.txt',
          sep = ''
        ),header = T,sep = '\t'
        
      )
      print(paste(Year_Exp,'canadian station considered.'))
      #d=d[complete.cases(d[,c(1,2,13)]),]
      nb_stations = length(unique(d$station_name))
      
      stations_close=as.data.frame(weathercan::stations_search(coords = c(latitude,longitude), dist = 70, starts_latest=2014, ends_earliest=2019))
      class(stations_close$start)='numeric'
      class(stations_close$end)='numeric'
      stations_close<-dplyr::filter(stations_close,start<2014&end>2019)
      stations_close=stations_close[stations_close$station_id%in%unique(d$station_id),]
      
      
      min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                       min(stations_close$distance)), 'station_id'])
      min_dist = unique(min(stations_close$distance))
      min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                    min(stations_close$distance)), 'station_name'])
      
      if (nrow(d[d$station_id == min_dist_station,]) > length(dates)) {
        nearest_station_values <- d[d$station_id == min_dist_station,]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), 'tmax']
      }
      while (length(nearest_station_values) < length(dates)|length(which(is.na(nearest_station_values)))>8) {
        stations_close = stations_close[-which(stations_close$station_id == min_dist_station),]
        min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                         min(stations_close$distance)), 'station_id'])
        min_dist = unique(min(stations_close$distance))
        min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                      min(stations_close$distance)), 'station_name'])
        
        nearest_station_values <- d[d$station_id == min_dist_station,]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), 'tmax']
        
      }
      
      
      
      
      
      network = 'WeatherCAN'
      print('Data from WeatherCAN stations prepared')
      
      ########################
      ########################
      
      
      
      
      
      dates = yday(seq(
        as.Date(date_start, tz = "CET"),
        as.Date(date_end, tz = "CET"),
        by = 'days'
      ))
      
      predicted.values <-
        nearest_station_values
      print(predicted.values)
      print(length(predicted.values))
      
      if (length(field_values) > 1 & !all(is.na(field_values))) {
        cors = cor(predicted.values, field_values, use = 'complete.obs')
      } else{
        cors = 'no_field_values'
      }
      
      predictions.table = cbind(Year_Exp, predicted.values, as.character(dates), network,field_values)
      predictions.table = as.data.frame(predictions.table)
      colnames(predictions.table) = c('Year_Exp', 'TMAX', 'date', 'network','field_values')
      
      to_save = list(predictions.table, cors, min_dist_station, min_dist_name,min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist_station_name',
          'min_dist'
        )
      
      saveRDS(to_save, paste('TMAX_nearest/', Year_Exp, '_1station.RDS', sep = ''))
      
      
      
      
    }
  }

setwd("C:/Users/cathyjubin/Desktop")
lapply(unique(daily_weather$Year_Exp),
       function(x)
         safeguarding(
           nearest_station_tmax(
             Year_Exp = x,
             daily_weather = daily_weather ,
             radius = 70
             
           )
         ))

