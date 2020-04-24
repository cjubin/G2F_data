#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGHCND} interpolates values for a specific meteorological variable given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable Character. GHCND element names: 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 


impute_kriging_withGHCND <- function(Year_Exp,radius=100,meteo_variable,daily_weather=daily_weather) {
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
  options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  print(Year_Exp)
  #Retrieve information about the experiment
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  #Values recorded for this variable at the field station
  field_values=daily_weather[daily_weather$Year_Exp==Year_Exp,meteo_variable]
  
  #Finding the closest stations in a certain radius and select lines (stations) with TEMP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = radius))
  stations_close<-filter(stations_close,element%in%meteo_variable)
  stations_close<-arrange(stations_close,distance)
  
  #Download data from stations exhibiting the meteo variable of interest
  
  download_data=function(station,datatypeid){
    dat=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid=datatypeid,startdate = date_start,enddate = date_end,limit = 500)$data
    if(length(dat)>0){
      values=dat$value/10
      dates=dat$date
      dist=stations_close[stations_close$id==station,c('distance')]
      longitude=unique(stations_close[stations_close$id==station,c('longitude')])
      latitude=unique(stations_close[stations_close$id==station,c('latitude')])
      
      
      d=cbind('station'=station,'longitude'=longitude,'latitude'=latitude,values,dates)
      
      
      return(d)}
  }
  
  
  
  all_data=lapply(stations_close$id, function(x)download_data(x,datatypeid = meteo_variable))
  
  all_data<-plyr::compact(all_data)
  all_data=as.data.frame(do.call('rbind',all_data))
  
  all_data$longitude=as.numeric(as.vector(all_data$longitude))  
  all_data$latitude=as.numeric(as.vector(all_data$latitude))  
  all_data$values=as.numeric(as.vector(all_data$values))
  all_data=arrange(all_data,dates)
  print('Data from GHCND stations prepared')
  
  ########################
  ####Ordinary kriging####
  
  sub=all_data
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
  
  #Simple Sum metric model
  simplesumMetric <- vgmST("simpleSumMetric",space = vgm(5,"Sph", 500, 0),time = vgm(500,"Sph", 500, 0), joint = vgm(1,"Sph", 500, 0), nugget=1, stAni=500) 
  
  #Automatic fit
  fitted.stvgm=fit.StVariogram(var,simplesumMetric,tunit="days",method="L-BFGS-B")
  attr(fitted.stvgm, "MSE")
  
  #Prediction grid
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
  cors=cor(predicted.values,field_values,use = 'complete.obs')
  
  
  dates=seq(as.Date(date_start,tz="CET"),as.Date(date_end,tz="CET"),by='days')
  predictions.table=cbind(Year_Exp,predicted.values,as.character(dates))
  colnames(predictions.table)=c('Year_Exp',paste(meteo_variable),'dates')
  return(list(predictions.table,cors))
  
  
  
  
  
  
}
