#' Download from NCEI the GSOD data files corresponding to the Year_Exp for the complete year (for further modeling of temporal correlation) and for stations present nearby the Year_Exp.
#' Recommend to download locally files (not on server)
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to download data.
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, a column 'long', and column 'lat'.



download_GSOD <- function(Year_Exp,radius=60,daily_weather=daily_weather) {
  library(rnoaa)
  library(data.table)
  source('safeguarding.R')
  options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  # ------------------------------------------------------------------------------
  #Retrieve information about the experiment
  # ------------------------------------------------------------------------------
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  # ------------------------------------------------------------------------------
  #Finding the closest stations in a certain radius and select those for which coverage period 2013-2019 is sure
  # ------------------------------------------------------------------------------
  stations_close=as.data.frame(rnoaa::isd_stations_search(lat =latitude, lon = longitude, radius = radius))
  stations_close$idstation=paste(stations_close$usaf,stations_close$wban,sep='')
  stations_close<-arrange(stations_close,distance)
  stations_close$yearstart<-substr(stations_close$begin,0,nchar(stations_close$begin)-4)
  stations_close$yearend<-substr(stations_close$end,0,nchar(stations_close$end)-4)
  stations_close<-filter(stations_close,yearstart<2013&yearend>2019)
  
  # ------------------------------------------------------------------------------
  #Download individual files for the stations of interest (selected before) from the ncei noaa server 
  # ------------------------------------------------------------------------------
  url_list <-
    CJ(year, stations_close$idstation, sorted = FALSE)[, paste0(
      "https://www.ncei.noaa.gov/data/global-summary-of-the-day/access/",
      year,
      "/",
      stations_close$idstation,
      ".csv"
    )]
  
  
  
  download_station <- function(x) {
    curl::curl_download(
      url = x,
      destfile = paste('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/GSOD/',year,'/', basename(x),sep=''),
      mode = "wb"
    )
  }
  data=lapply(url_list,function(x)safeguarding(download_station(x)))
  
  
}




