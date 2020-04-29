#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable_GSOD Character. GSOD element names (= variable measured at the ISD STATION): TEMP, MAX, MIN, DEWP, STP, WDSP, PRCP
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested', and optionally meteo_variable_in_table with its value as column name. Ex: 'TMIN'.
#' @param meteo_variable_in_table. Character with the column name in the table daily_weather of the meteorological variable of interest which has to be imputed from surrounding stations. It does not have to be the same as @meteo_variable_GSOD



impute_kriging_withGSOD <- function(Year_Exp,radius=50,meteo_variable_GSOD,daily_weather=daily_weather,meteo_variable_in_table=NULL,only_download=FALSE) {
  library(rnoaa)
  library(raster)
  library(mapdata)
  library(maps)
  library(maptools)
  library(sp)
  library(gstat)
  library(xts)
  library(spacetime)
  library(raster)
  library(rgdal)
  source('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/fahrenheit_to_celsius.R')
  
  options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  print(Year_Exp)
  #Retrieve information about the experiment
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  #Values recorded for this variable at the field station
  if (!is.null(meteo_variable_in_table)) {
    field_values = daily_weather[daily_weather$Year_Exp == Year_Exp, meteo_variable_in_table]
  }
  
  #Finding the closest stations in a certain radius and select those for which coverage period 2013-2019 is sure
  stations_close=as.data.frame(rnoaa::isd_stations_search(lat =latitude, lon = longitude, radius = radius))
  stations_close$idstation=paste(stations_close$usaf,stations_close$wban,sep='')
  stations_close<-arrange(stations_close,distance)
  stations_close$yearstart<-substr(stations_close$begin,0,nchar(stations_close$begin)-4)
  stations_close$yearend<-substr(stations_close$end,0,nchar(stations_close$end)-4)
  stations_close<-filter(stations_close,yearstart<2013&yearend>2019)
  
  #Download data from stations exhibiting the meteo variable of interest and set data.frame for kriging 
  
  #Download individual files for the stations of interest (selected before) from the ncei noaa server 
  url_list <-
    CJ(year, stations_close$idstation, sorted = FALSE)[, paste0(
      "https://www.ncei.noaa.gov/data/global-summary-of-the-day/access/",
      year,
      "/",
      stations_close$idstation,
      ".csv"
    )]
  
  safe_download <- function(x) {
    tryCatch(
      download_station(x),
      warning = function(w)
        NULL,
      error = function(e)
        NULL
    )
  }
  download_station <- function(x) {
    curl::curl_download(
      url = x,
      destfile = file.path(tempdir(),year, basename(x)),
      mode = "wb"
    )
  }
  data=lapply(url_list,function(x)safe_download(x))
  
  
  #Option: only download data and stop there the function.
  if (only_download==TRUE){
    break}
  else{  
    
  #files_stations contains the path (if data for theses ISD stations is available) where to find on the computer the data for the year ans station(s) of interest.
  files_stations <-
    paste0(tempdir(),
           "/",
           year,
           '/',
           stations_close$idstation,
           ".csv")
  #List of the csv files in the temporary folder
  GSOD_list <-
    list.files(
      tempdir(),
      pattern = "*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
  
  #Overlap between all initial selected stations and those for which data could be eventually downloaded and are on the computer.
  GSOD_list2 <-
    subset(GSOD_list, GSOD_list %in% files_stations)
  
  data_frames=lapply(GSOD_list2,function(x)read.csv(x))
  all_data<-plyr::compact(data_frames)
  all_data=as.data.frame(do.call('rbind',all_data))
    
  d = cbind(
    'station' = all_data$STATION,
    'longitude' = all_data$LONGITUDE,
    'latitude' = all_data$LATITUDE,
    'values' = all_data[, colnames(all_data) %in% meteo_variable_GSOD],
    'dates' = as.character(as.vector(all_data$DATE))
  )
  d=as.data.frame(d) 
  
  d$longitude=as.numeric(as.vector(d$longitude))  
  d$latitude=as.numeric(as.vector(d$latitude)) 
  d$values=as.numeric(as.vector(d$values))
  
  if (meteo_variable_GSOD%in%c('TEMP', 'DEWP', 'MAX', 'MIN')) {
    d$values <- fahrenheit_to_celsius(d$values)
  }
  
  d$dates=as.Date(d$dates)
  d=arrange(d,dates)
  print('Data from GSOD stations prepared')
  
  
  ########################
  ####Ordinary kriging####
  
  sub=d
  sp::coordinates(sub)=c('longitude','latitude')
  proj4string(sub) = "+proj=longlat +datum=WGS84"
  #projection(sub)=CRS("+init=epsg:4326")
  
  #Transform into Mercator Projection
  tempmin.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 
  
  
  tempminSP <- SpatialPoints(tempmin.UTM@coords,CRS("+init=epsg:3395"))
  tempminDF <- data.frame(values=tempmin.UTM$values) 
  tempminTM <- as.POSIXct(date(tempmin.UTM$dates))
  #combine the 3 objects
  timeDF <- STIDF(tempminSP,tempminTM,data=tempminDF) 
  
  
  #variogram
  print('Computation variogram starts:')
  var <- variogramST(values~1,data=timeDF,tunit="days",assumeRegular=F,na.omit=T) 
  
  plot(var,map=F)
  
  #Sum metric model
  sumMetric <-vgmST("sumMetric", space = vgm(psill=5,"Sph", range=500, nugget=0),time = vgm(psill=500,"Sph", range=500, nugget=0), joint = vgm(50,"Mat", range=500, nugget=10), stAni=1) 
  #Automatic fit
  fitted.stvgm=fit.StVariogram(var,sumMetric)
  attr(fitted.stvgm, "MSE")
  par(mfrow=c(2,1))
  plot(var,fitted.stvgm,map=F) 
  
  #Prediction grid: growing season for the field experiment described by Year_Exp
  field=vector(mode = 'numeric',length = 2)
  field<-c(longitude,latitude)
  
  field=as.data.frame(matrix(field,ncol = 2))
  colnames(field)<-c('longitude','latitude')
  sp::coordinates(field)=c('longitude','latitude')
  proj4string(field) = "+proj=longlat +datum=WGS84"
  field_grid<- spTransform(field,CRS("+init=epsg:3395")) 
  tm.grid<-seq(as.POSIXct(date_start,tz="CET"),as.POSIXct(date_end,tz="CET"),by='days')
  grid.ST <- STF(field_grid,tm.grid) 
  
  
  pred<-krigeST(values~1,data=timeDF,modelList =fitted.stvgm,newdata = grid.ST )
  predicted.values=pred@data$var1.pred
  if (!is.null(meteo_variable_in_table)) {
  cors=cor(predicted.values,field_values,use = 'complete.obs')}
  
  
  dates=seq(as.Date(date_start,tz="CET"),as.Date(date_end,tz="CET"),by='days')
  predictions.table=cbind(Year_Exp,predicted.values,as.character(dates))
  colnames(predictions.table)=c('Year_Exp',paste(meteo_variable_GSOD),'dates')
  return(list(predictions.table,cors))
  
  
  }
  
  
  
}


safe_impute_function<- function(y,radius,meteo_variable_GSOD,daily_weather=daily_weather,meteo_variable_in_table,only_download) {
  tryCatch(
    impute_kriging_withGSOD(Year_Exp = y,radius=radius,meteo_variable_GSOD = meteo_variable_GSOD,meteo_variable_in_table=meteo_variable_in_table,daily_weather = daily_weather,only_download = only_download ),
    #warning = function(w) NULL,
    error = function(e) NULL
  )
}



#all_experiments=unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp) %notin%
#                                                 c('2014_ONH1',
#                                                   '2014_ONH2',
#                                                   '2015_ONH1',
#                                                   '2015_ONH2',
#                                                   '2016_ONH1',
#                                                   '2016_ONH2')]

#require(doParallel)
#workers <- makeCluster(3) 
#registerDoParallel(workers)
#results=mclapply(all_experiments[1:3],
#                 function(x)
#                   safe_impute_function(x,radius=50,meteo_variable_GSOD = 'MIN',meteo_variable_in_table  ='TMIN',daily_weather = daily_weather))
