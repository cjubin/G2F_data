#' Download from NOAA the GHCNDdata files corresponding to the Year_Exp for the complete year (for further modeling of temporal correlation) and for stations present nearby the Year_Exp.
#' Recommend to download locally files (not on server)
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from GHCND stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, a column 'long', and column 'lat'.
#' @param meteo_variable_GHCND. Character. Name of the variable used in the GHCND meteorological element to retrieve. List of meteorological elements listed in 'GHCND_documentation.pdf'


download_GHCND <- function(Year_Exp,radius=60,daily_weather=daily_weather,meteo_variable_GHCND) {
  library(rnoaa)
  library(data.table)
  source('safeguarding.R')
  options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  #Retrieve information about the experiment
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  
  #Finding the closest stations in a certain radius and select those for which coverage period 2013-2019 is sure
  #Finding the closest stations in a certain radius and select lines (stations) with TEMP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = radius))
  stations_close<-filter(stations_close,element%in%meteo_variable_GHCND)
  stations_close<-arrange(stations_close,distance)
  
  #Download data from stations exhibiting the meteo variable of interest and set data.frame for kriging 
  
  #Download individual files for the stations of interest (selected before) from the NOAA data center
  
  download_data=function(station,datatypeid){
    dat=ghcnd_search(stationid = station,var=datatypeid,date_min = as.Date(paste0(year,'-01-01')),date_max  = as.Date(paste0(year,'-12-31')))[[1]]
    dat=dat[,-which(colnames(dat)%in%c('mflag','cflag','qflag'))]
    if(length(dat)>0){
      
      dat[,2]=dat[,2]/10
      dat$longitude=unique(stations_close[stations_close$id==station,c('longitude')])
      dat$latitude=unique(stations_close[stations_close$id==station,c('latitude')])
      return(dat)}
    
    else{return(NULL)}
  }
  
  
  
  all_data=lapply(stations_close$id,function(x)safeguarding(download_data(x,datatypeid=meteo_variable_GHCND)))
  
  #Organize summary data.frame with all observations of surrounding stations
  all_data<-plyr::compact(all_data)
  all_data=as.data.frame(do.call('rbind',all_data))
  
  all_data$longitude=as.numeric(as.vector(all_data$longitude))  
  all_data$latitude=as.numeric(as.vector(all_data$latitude))  
  all_data$values=as.numeric(as.vector(all_data$values))
  all_data=arrange(all_data,date)
  
  all_data$Year_Exp=Year_Exp
  all_data$variable=meteo_variable_GHCND
  
  
  write.table(all_data,file=paste('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/GHCND/',year,'/',meteo_variable_GHCND,'_',Year_Exp,'.txt',sep = ''))
  
  return(all_data)
  
}




