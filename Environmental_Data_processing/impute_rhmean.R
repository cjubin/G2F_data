#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_rhmean} interpolates values for a specific meteorological variable using Inverse Distance Weighting or nearby station if only 1 station.
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to capture GSOD stations.
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested' and @variable_to_impute
#' @param meteo_variable_GSOD. Character. Name of the variable in the GSOD datasets
#' @param name_in_table. Character. Name of the variable in the original G2F datasets

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


`%notin%` <- Negate(`%in%`)


library(dplyr)
library(lubridate)

setwd(
  "~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Env_data_processing"
)
daily_weather = read.table(
  '2_merged_dataset.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
daily_weather = dplyr::arrange(daily_weather, Year_Exp, Day.of.Year)

daily_weather[daily_weather$Year_Exp == '2017_MNH1', 'lat'] <-
  44.06981
daily_weather[daily_weather$Year_Exp == '2017_MNH1', 'long'] <-
  -93.5338


#stations <- ghcnd_stations()
#stations <- filter(stations, last_year == 2020)
#stations <- filter(stations, first_year <= 2013)



impute_rhmean <-
  function(Year_Exp,
           radius = 70,
           meteo_variable_GSOD = c('TEMP', 'DEWP', 'MAX', 'MIN'),
           daily_weather = daily_weather,
           name_in_table = 'HMEAN') {
    options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
    
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
    
    
    
    ## Dewpoint data taken from GSOD network: more data.
    
    if (Year_Exp %notin% c(
      '2014_ONH1',
      '2015_ONH1',
      '2016_ONH1',
      '2017_ONH1',
      '2014_ONH2',
      '2015_ONH2',
      '2016_ONH2',
      '2017_ONH2',
      '2018_ONH2'
    )) {
      # Take dewpoint and temp data from GSOD stations
      variable_to_impute = 'HMEAN'
      
      stations_close = as.data.frame(rnoaa::isd_stations_search(
        lat = latitude,
        lon = longitude,
        radius = radius
      ))
      stations_close$idstation = paste(stations_close$usaf, stations_close$wban, sep =
                                         '')
      stations_close <- arrange(stations_close, distance)
      
      
      
      
      stations_close$yearstart <-
        substr(stations_close$begin, 0, nchar(stations_close$begin) - 4)
      stations_close$yearend <-
        substr(stations_close$end, 0, nchar(stations_close$end) - 4)
      stations_close <- filter(stations_close, yearstart < 2013 &
                                 yearend > 2019)
      
      #files_stations: path to find weather files for the stations used in kriging (radius 80 km) (if data for theses ISD stations is available)
      
      files_stations <-
        paste0(
          'C:/Users/cathyjubin/Desktop/hmean/gsod_files',
          "/",
          year,
          '/',
          stations_close$idstation,
          ".csv"
        )
      #List of the csv files in the temporary folder
      GSOD_list <-
        list.files(
          paste0('C:/Users/cathyjubin/Desktop/hmean/gsod_files',
                 '/',
                 year),
          pattern = "*\\.csv$",
          full.names = TRUE,
          recursive = TRUE
        )
      
      #Overlap between all initial selected stations and those for which data could be eventually downloaded on the computer.
      GSOD_list2 <-
        subset(GSOD_list, GSOD_list %in% files_stations)
      
      data_frames = lapply(GSOD_list2, function(x) {
        read.csv(x)
      })
      all_data <- plyr::compact(data_frames)
      all_data = as.data.frame(do.call('rbind', all_data))
      
      
      nb_stations = length(unique(all_data$STATION))
      stations_close = stations_close[stations_close$idstation %in% all_data$STATION, ]
      
      
      min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                       min(stations_close$distance)), 'idstation'])
      min_dist = unique(min(stations_close$distance))
      min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                    min(stations_close$distance)), 'station_name'])
      
      #### HMEAN must me computed based on TEMP and DEWP.
      
      if (any(meteo_variable_GSOD %in% c('DEWP'))) {
        index = which(all_data$DEWP > 999)
        if (length(index) > 0) {
          all_data <- all_data[-index, ]
        }
        
      }
      
      for (var in meteo_variable_GSOD) {
        all_data[, eval(var)] <-
          fahrenheit_to_celsius(as.numeric(as.vector(all_data[, eval(var)])))
      }
      
      all_data$HMEAN <- 100 * (exp((17.625 * all_data$DEWP) / (243.04 + all_data$DEWP)) / exp((17.625 * all_data$TEMP) /
                                                                                                (243.04 + all_data$TEMP)))
      
      
      #####
      
      
      if (nrow(all_data[all_data$STATION == min_dist_station, ]) > length(field_values)) {
        nearest_station_values <-
          all_data[all_data$STATION == min_dist_station, ]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$DATE) %in% yday(dates), 'HMEAN']
        
      }
      while (length(nearest_station_values) < length(field_values)) {
        stations_close = stations_close[-which(stations_close$idstation == min_dist_station), ]
        min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                         min(stations_close$distance)), 'idstation'])
        min_dist = unique(min(stations_close$distance))
        min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                      min(stations_close$distance)), 'station_name'])
        
        
        nearest_station_values <-
          all_data[all_data$STATION == min_dist_station, ]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$DATE) %in% yday(dates), 'HMEAN']
        
      }
      
      d = cbind(
        'longitude' = all_data$LONGITUDE,
        'latitude' = all_data$LATITUDE,
        all_data[, colnames(all_data) %in% 'HMEAN'],
        'dates' = as.character(as.vector(all_data$DATE))
      )
      d = as.data.frame(d)
      
      
      colnames(d)[3] = 'HMEAN'
      
      
      d$longitude = as.numeric(as.vector(d$longitude))
      d$latitude = as.numeric(as.vector(d$latitude))
      
      
      d$dates = as.Date(d$dates)
      d$dates = lubridate::yday(d$dates)
      
      for (i in colnames(d)[colnames(d) %notin% c('dates')]) {
        d[, i] <- as.numeric(as.vector(d[, i]))
      }
      
      
      
      network = 'GSOD'
      
      print('Data from GSOD stations prepared')
      
      
    }  else{
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
          nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), c('rhmax','rhmin')]
      }
      while (nrow(nearest_station_values) < length(dates)|length(which(is.na(nearest_station_values)))>20) {
        stations_close = stations_close[-which(stations_close$station_id == min_dist_station),]
        min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                         min(stations_close$distance)), 'station_id'])
        min_dist = unique(min(stations_close$distance))
        min_dist_name = unique(stations_close[which(stations_close$distance ==
                                                      min(stations_close$distance)), 'station_name'])
        
        nearest_station_values <- d[d$station_id == min_dist_station,]
        nearest_station_values <-
          nearest_station_values[yday(nearest_station_values$date) %in% yday(dates), c('rhmax','rhmin')]
        
      }
      
      
      
      d = cbind(
        'longitude' = d$lon,
        'latitude' = d$lat,
        d[, colnames(d) %in% c('rhmax','rhmin')],
        'dates' = as.character(as.vector(d$date))
      )
      d = as.data.frame(d)
      
      colnames(d)[3:4]=c('rhmax','rhmin')
      
      d$longitude = as.numeric(as.vector(d$longitude))
      d$latitude = as.numeric(as.vector(d$latitude))
      
      
      d$dates = as.Date(d$dates)
      d$dates = lubridate::yday(d$dates)
      
      for (i in colnames(d)[colnames(d) %notin% c('dates')]) {
        d[, i] <- as.numeric(as.vector(d[, i]))
      }
      
      
      
      
      network = 'WeatherCAN'
      print('Data from WeatherCAN stations prepared')}
    ########################
    ########################
    
    
    
    if (nb_stations == 1) {
      print(paste(Year_Exp, 'presents only 1 station with humidity data.'))
      print(length(unique(d$id)))
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
      
      predictions.table = cbind(Year_Exp, predicted.values, as.character(dates), network)
      predictions.table = as.data.frame(predictions.table)
      colnames(predictions.table) = c('Year_Exp', 'HMEAN', 'date', 'network')
      
      to_save = list(predictions.table, cors, min_dist_station, min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist'
        )
      
      saveRDS(to_save, paste('HMEAN_GSOD_WeatherCAN/', Year_Exp, '_1station.RDS', sep = ''))
    }
    
    
    
    
    ### IDW ###
    
    if (nb_stations > 1) {
      
      
      if(Year_Exp%in%c(
        '2014_ONH1',
        '2015_ONH1',
        '2016_ONH1',
        '2017_ONH1',
        '2014_ONH2',
        '2015_ONH2',
        '2016_ONH2',
        '2017_ONH2',
        '2018_ONH2'
      )){
        
        
        d = d[complete.cases(d), ]
        pred_grid_field <- expand.grid(
          longitude = longitude,
          latitude = latitude,
          dates = lubridate::yday(substr(dates, 1, 10))
        )
        
        pred1 <- gstat::idw(
          formula = rhmin ~ 1,
          data = d,
          locations = ~ latitude + longitude + dates,
          newdata = pred_grid_field,
          idp = 2
        )
        pred2 <- gstat::idw(
          formula = rhmax ~ 1,
          data = d,
          locations = ~ latitude + longitude + dates,
          newdata = pred_grid_field,
          idp = 2
        )
        pred = cbind(pred1[, 1:4], pred2[, 4])
        colnames(pred)[4:5] = c('rhmin', 'rhmax')
        pred$nb_stations_used = nb_stations
        pred$network = network
        pred$nearest_station_values <- nearest_station_values
        pred$field_Values <- field_values
        
        if (length(field_values) > 1 & !all(is.na(field_values))) {
          cors = cor(pred$var1.pred, field_values, use = 'complete.obs')
        } else{
          cors = 'no_field_values'
        }
        
        
        to_save <-
          list(pred, cors, min_dist_station, min_dist_name, min_dist)
        names(to_save) <-
          c(
            'predictions_YearExp',
            'cors_YearExp_with1station',
            'min_dist_station',
            'min_dist_name',
            'min_dist'
          )
        
        
        saveRDS(to_save,
                paste('HMEAN_GSOD_WeatherCAN/', Year_Exp, '_IDW.RDS', sep = ''))
        
        
        
        
      } else {print(nb_stations)
      
      
      library(gstat)
      
      
      d = d[complete.cases(d),]
      d = d[which(d$HMEAN>0&d$HMEAN<100),]
      
      
      pred_grid_field <- expand.grid(
        longitude = longitude,
        latitude = latitude,
        dates = lubridate::yday(substr(dates, 1, 10))
      )
      
      pred <- gstat::idw(
        formula = HMEAN ~ 1,
        data = d,
        locations = ~ latitude + longitude + dates,
        newdata = pred_grid_field,
        idp = 2
      )
      pred$nb_stations_used = nb_stations
      pred$network = network
      pred$nearest_station_values <- nearest_station_values
      pred$field_Values <- field_values
      
      if (length(field_values) > 1 & !all(is.na(field_values))) {
        cors = cor(pred$var1.pred, field_values, use = 'complete.obs')
      } else{
        cors = 'no_field_values'
      }
      
      
      to_save <- list(pred, cors, min_dist_station, min_dist_name,min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist_name',
          'min_dist'
        )
      
      
      saveRDS(to_save, paste('HMEAN_GSOD_WeatherCAN/', Year_Exp, '_IDW.RDS', sep = ''))
      
      }
    }
    
  }


lapply(unique(daily_weather$Year_Exp),
       function(x)
         safeguarding(
           impute_rhmean(
             Year_Exp = x,
             daily_weather = daily_weather ,
             radius = 70
             
           )
         ))

