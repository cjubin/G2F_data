#' Download from NOAA the GHCNDdata files corresponding to the Year_Exp for the complete year (for further modeling of temporal correlation) and for stations present nearby the Year_Exp.
#' Recommend to download locally files (not on server)
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from GHCND stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to download data.
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, a column 'long', and column 'lat'.
#' @param meteo_variable_GHCND. Character. Name of the variable used in the GHCND meteorological element to retrieve. List of meteorological elements listed in 'GHCND_documentation.pdf'
#' @param stations download the initial file containing all the GHCND stations

download_GHCND_allstations_elev <- function(Year_Exp,radius=70,daily_weather=daily_weather,meteo_variable_GHCND,stations) {
  library(rnoaa)
  library(data.table)
  library(elevatr)
  source('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Env_data_processing/safeguarding.R')
  options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  # ------------------------------------------------------------------------------
  #Retrieve information about the experiment
  # ------------------------------------------------------------------------------
  
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  # ------------------------------------------------------------------------------
  #Finding the closest stations in a certain radius and select those with the meteorologcail variables of interest
  # ------------------------------------------------------------------------------
  
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = radius))
  stations_close<-dplyr::filter(stations_close,element%in%meteo_variable_GHCND)
  stations_close<-arrange(stations_close,distance)
  
  
  # ------------------------------------------------------------------------------
  #Download individual files for the stations selected previously from the NOAA data center
  # ------------------------------------------------------------------------------
  
  download_data=function(station,datatypeid){
    dat=ghcnd_search(stationid = station,var=datatypeid,date_min = as.Date(paste0(year,'-01-01')),date_max  = as.Date(paste0(year,'-12-31')))[[1]]
    if(nrow(dat)>0){
      
      dat[,'prcp']=dat[,'prcp']/10
      dat$longitude=unique(stations_close[stations_close$id==station,c('longitude')])
      dat$latitude=unique(stations_close[stations_close$id==station,c('latitude')])
      return(dat)
      return(dat)}
    else{return(NULL)}
  }
  
  
  
  all_data=lapply(stations_close$id,function(x)safeguarding(download_data(x,datatypeid=meteo_variable_GHCND)))
  
  # ------------------------------------------------------------------------------
  #Organize summary data.frame with all observations of surrounding stations
  # ------------------------------------------------------------------------------
  
  all_data<-plyr::compact(all_data)
  all_data=as.data.frame(do.call('rbind',all_data))
  
  all_data$longitude=as.numeric(as.vector(all_data$longitude))  
  all_data$latitude=as.numeric(as.vector(all_data$latitude)) 
  loc1=unique(all_data[,c('id','longitude','latitude')])
  loc=unique(all_data[,c('longitude','latitude')])
  
  
  
  print('done')
  elev<-get_elev_point(locations=loc,units='meters',
                           prj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  elev_df=cbind(loc1,elev@data)
  
  all_data$elev=elev_df[match(all_data$id,loc1$id),'elevation']
  
  
  all_data=arrange(all_data,date)
  
  all_data$Year_Exp=Year_Exp
  all_data$variable=meteo_variable_GHCND
  
  if (length(unique(all_data$id))<6){cat(paste('Less than 6 stations available'))}
  
  write.table(all_data,file=paste('C:/Users/cathyjubin/Desktop/PRCP_allstations_elev/',meteo_variable_GHCND,'_radius70km_',Year_Exp,'.txt',sep = ''),col.names = T,row.names = F,sep = '\t',quote = F)
  if (nrow(all_data)>1){cat(paste('Files for',Year_Exp,'assembled and written'))}
  
  return(all_data)
  
}
daily_weather=read.table(
  'subset_1.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)





lapply(unique(daily_weather$Year_Exp),
       function(x)
         download_GHCND_allstations_elev (
           x,
           radius = 70,
           daily_weather = daily_weather,
           stations=stations,
           meteo_variable_GHCND='PRCP'
         ))
