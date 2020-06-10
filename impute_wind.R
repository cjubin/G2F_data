#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGHCND} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
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


stations <- ghcnd_stations()
stations <- filter(stations, last_year == 2020)
stations <- filter(stations, first_year <= 2013)



impute_wind2 <-
  function(Year_Exp,
           radius = 70,
           meteo_variable_GHCND = 'awnd',
           meteo_variable_GSOD = 'WDSP',
           daily_weather = daily_weather,
           name_in_table = 'MEANWINDSPEED') {
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
    }
    
    
    
    ## Wind data taken from GSOD network: more data.
    
    
    # Take wind data from GSOD stations
    variable_to_impute = 'WDSP'
    
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
    stations_close=stations_close[stations_close$idstation%in%all_data$STATION,]
    
    
    min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                     min(stations_close$distance)), 'idstation'])
    min_dist = unique(min(stations_close$distance))
    
    
    if (nrow(all_data[all_data$STATION == min_dist_station,]) > length(field_values)) {
      nearest_station_values <-
        all_data[all_data$STATION == min_dist_station,]
      nearest_station_values <-
        nearest_station_values[yday(nearest_station_values$DATE) %in% yday(dates), 'WDSP']
      
    }
    while (length(nearest_station_values) < length(field_values)) {
      stations_close = stations_close[-which(stations_close$idstation == min_dist_station),]
      min_dist_station = unique(stations_close[which(stations_close$distance ==
                                                       min(stations_close$distance)), 'idstation'])
      min_dist = unique(min(stations_close$distance))
      
      nearest_station_values <-
        all_data[all_data$STATION == min_dist_station,]
      nearest_station_values <-
        nearest_station_values[yday(nearest_station_values$DATE) %in% yday(dates), 'WDSP']
      
    }
    
    d = cbind(
      'longitude' = all_data$LONGITUDE,
      'latitude' = all_data$LATITUDE,
      all_data[, colnames(all_data) %in% meteo_variable_GSOD],
      'dates' = as.character(as.vector(all_data$DATE))
    )
    d = as.data.frame(d)
    
    if (length(meteo_variable_GSOD) == 1) {
      colnames(d)[3] = meteo_variable_GSOD
    }
    
    d$longitude = as.numeric(as.vector(d$longitude))
    d$latitude = as.numeric(as.vector(d$latitude))
    
    
    d$dates = as.Date(d$dates)
    d$dates = lubridate::yday(d$dates)
    
    for (i in colnames(d)[colnames(d) %notin% c('dates')]) {
      d[, i] <- as.numeric(as.vector(d[, i]))
    }
    
    
    
    network = 'GSOD'
    
    print('Data from GSOD stations prepared')
    
    ########################
    ########################
    
    
    
    if (nb_stations == 1) {
      print(paste(Year_Exp, 'presents only 1 station with wind data.'))
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
      colnames(predictions.table) = c('Year_Exp', 'WDSP', 'date', 'network')
      predictions.table$transformed_u2 <-
        0.748 * as.numeric(as.vector(predictions.table$WDSP))
      
      to_save = list(predictions.table, cors, min_dist_station, min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist'
        )
      
      saveRDS(to_save, paste('AWND_GSOD/', Year_Exp, '_1station.RDS', sep = ''))
    }
    
    
    
    
    ### IDW ###
    
    if (nb_stations > 1) {
      print(nb_stations)
      
      
      library(gstat)
      
      
      d = d[complete.cases(d),]
      
      
      pred_grid_field <- expand.grid(
        longitude = longitude,
        latitude = latitude,
        dates = lubridate::yday(substr(dates, 1, 10))
      )
      
      pred <- gstat::idw(
        formula = WDSP ~ 1,
        data = d,
        locations = ~ latitude + longitude + dates,
        newdata = pred_grid_field,
        idp = 2
      )
      pred$nb_stations_used = nb_stations
      pred$network = network
      pred$nearest_station_values <- nearest_station_values
      pred$field_Values <- field_values
      pred$transformed_u2 <- 0.748 * pred$var1.pred
      
      if (length(field_values) > 1 & !all(is.na(field_values))) {
        cors = cor(pred$var1.pred, field_values, use = 'complete.obs')
      } else{
        cors = 'no_field_values'
      }
      
      
      to_save <- list(pred, cors, min_dist_station, min_dist)
      names(to_save) <-
        c(
          'predictions_YearExp',
          'cors_YearExp_with1station',
          'min_dist_station',
          'min_dist'
        )
      
      
      saveRDS(to_save, paste('AWND_GSOD/', Year_Exp, '_IDW.RDS', sep = ''))
      
      
    }
    
  }


lapply(unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin% '2018_ONH2'],
       function(x)
         safeguarding(
           impute_wind2(
             Year_Exp = x,
             daily_weather = daily_weather ,
             radius = 70
             
           )
         ))
