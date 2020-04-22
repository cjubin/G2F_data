rm(list=ls())
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)

library(countyweather)


library(dplyr)
#library(plyr)
library(lubridate)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
stations <- ghcnd_stations() 
stations<-filter(stations,last_year==2020)
stations<-filter(stations,first_year<=2013)

weather = read.table(
  '~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

daily_weather=read.table(
  '~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/daily_weather.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

##TEMPERATURE: flagged values to be imputed

##Control range semi-hourly or 15-, 20-minutes data values
weather$flagged_temp_instant=NA
weather[which(weather$Temperature..C.>60),'flagged_temp_instant']<-'flagged'
weather[which(weather$Temperature..C.>60),'Temperature..C.']<-NA
weather[which(weather$Temperature..C.<(-40)),'flagged_temp_instant']<-'flagged'
weather[which(weather$Temperature..C.<(-40)),'Temperature..C.']<-NA

#First flagged value assignment based on ratio temp values present/nb obs per day + if weird number of daily obs (ex: 50, 100 etc) then the value is flagged too.
weather$Year_Exp=as.factor(weather$Year_Exp)
weather$Day.of.Year=as.factor(weather$Day.of.Year)
total_obs_daily<-as.data.frame(weather %>% group_by(Year_Exp,Day.of.Year) %>% tally())

nomissing_temp<-weather %>% 
  aggregate(Temperature..C.~Day.of.Year+Year_Exp,data=.,FUN = (function(x) {sum(!is.na(x))})) 
colnames(nomissing_temp)[3]='number_obs_present'

alldaymissing_temp<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all(is.na(Temperature..C.)))
colnames(alldaymissing_temp)[3]='all_missing_temp'



temp<-merge(total_obs_daily,nomissing_temp,by=c('Day.of.Year','Year_Exp'),all.x=T)
temp<-merge(temp,alldaymissing_temp,by=c('Day.of.Year','Year_Exp'),all.x=T)

rm(total_obs_daily,nomissing_temp,alldaymissing_temp)

temp<-temp%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=number_obs_present/n)%>%mutate(flagged_temp=case_when(ratio<0.9 ~'flagged',
                                                                                                                 ratio>=0.9 ~'OK'))

temp[temp$n<24,'flagged_temp']='flagged'
temp[temp$n>96,'flagged_temp']='flagged'
temp[temp$n>24&temp$n<48,'flagged_temp']='flagged'
temp[temp$n>48&temp$n<72,'flagged_temp']='flagged'
temp[temp$n>72&temp$n<96,'flagged_temp']='flagged'
temp[temp$all_missing_temp==TRUE,'flagged_temp']='flagged'
temp[temp$all_missing_temp==TRUE,'ratio']=0

weather=merge(weather,temp[,c(1:2,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)

daily_weather=merge(daily_weather,temp[,c(1:3,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)



#Daily mean,max and min computed for unflagged values:


maxT<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(max_temp=max(Temperature..C.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,max_temp)
maxT<-unique(maxT)
minT<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(min_temp=min(Temperature..C.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,min_temp)
minT<-unique(minT)
temperatures=merge(maxT,minT,by=c('Day.of.Year','Year_Exp'),all.x=T)

temperatures=arrange(temperatures,Year_Exp,Day.of.Year)



#Internal consistency test: Tmax(d) > Tmin(d-1) + Tmin(d) < Tmax(d-1)
library(data.table)
nm1 <- c('max_temp','min_temp')
nm2 <- paste("lag", nm1, sep=".")

library(dplyr)
test  <- 
  temperatures %>%
  group_by(Year_Exp) %>%
  mutate(lag = dplyr::lag(min_temp, n = 1, default = NA))%>%
  mutate(diff = max_temp-lag)%>%
  mutate(lag2 = dplyr::lag(max_temp, n = 1, default = NA))%>%
  mutate(diff2 = min_temp-lag2)%>%
  filter(diff>0)%>%
  filter(diff2<0)%>%
  select(Year_Exp,Day.of.Year,min_temp,max_temp)

temperatures<-test


#Mean temperature: WMO (2010) recommends use of this estimator:'Even though this method is not the best statistical approximation, its consistent use satisfies the comparative purpose of normal'

temperatures$mean_temp=(temperatures$max_temp+temperatures$min_temp)/2



#Add this daily temp sum to the daily_weather table
daily_weather<-merge(daily_weather,temperatures,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)





######IMPUTE temp MISSING VALUES #####



daily_weather$stationID_NOAA=NA
daily_weather$dist=NA
cc=1
par(mfrow=c(2,2))
`%notin%` <- Negate(`%in%`)



impute_temp_idw=function(x=Year_Exp,radius=50){
  print(x)
  #Retrieve information about the experiment
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==x,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==x,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==x,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==x,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==x,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==x,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==x,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==x,'lat']))])
  
  #Finding the closest stations in a certain radius and select lines (stations) with TEMP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = radius))
  stations_close<-filter(stations_close,element%in%c('TMAX','TMIN'))
  stations_close<-arrange(stations_close,distance)
  
  #Some stations do not exhibit data for all days requested (days of the growing season)
  find_station=function(x){
    id_stations=stations_close[x,'id']
    values_tmax= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values_tmin= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values=data.frame(values_tmax,values_tmin)
    if(nrow(values)==length(yday(date_start):yday(date_end))){
      return(id_stations)
    }
    else{return(NULL)}
  }
  
  find_station2=function(x){
    id_stations=stations_close[x,'id']
    values_tmax= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values_tmin= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values=data.frame(values_tmax,values_tmin)
    if(nrow(values)>length(yday(date_start):yday(date_end))-5){
      return(id_stations)
    }
    else{return(NULL)}
  }
  
  safe_find_function<- function(x) {
    tryCatch(
      find_station(x),
      warning = function(w) NULL,
      error = function(e) NULL
    )
  }
  
  safe_find_function2<- function(x) {
    tryCatch(
      find_station2(x),
      warning = function(w) NULL,
      error = function(e) NULL
    )
  }
  
  list_stations_no_missingdates=unique(unlist(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))))
  
  
  if(length(list_stations_no_missingdates)!=0) {
    download_data=function(station,datatypeid){
      values=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid=datatypeid,startdate = date_start,enddate = date_end,limit = 500)$data$value/10
      dates=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid=datatypeid,startdate = date_start,enddate = date_end,limit = 500)$data$date
      dist=stations_close[stations_close$id==station,c('distance')]
      longitude=unique(stations_close[stations_close$id==station,c('longitude')])
      latitude=unique(stations_close[stations_close$id==station,c('latitude')])
      
      
      d=matrix(c('station'=station,'longitude'=longitude,'latitude'=latitude,values),nrow = 1)
      colnames(d)<-c('station','longitude','latitude',dates)
      
      return(d)
    }
      
    
    tmaxdata=do.call('rbind',lapply(list_stations_no_missingdates, function(x)download_data(x,datatypeid = 'TMAX')))
    tmindata=do.call('rbind',lapply(list_stations_no_missingdates, function(x)download_data(x,datatypeid = 'TMIN')))
    
  }
  
  
  if(length(list_stations_no_missingdates)==0) {
    list_stations_possible_missingdates=unique(unlist(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function2))))
    
    download_data=function(station,datatypeid){
      values=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid=datatypeid,startdate = date_start,enddate = date_end,limit = 500)$data$value/10
      dates=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid=datatypeid,startdate = date_start,enddate = date_end,limit = 500)$data$date
      dist=stations_close[stations_close$id==station,c('distance')]
      longitude=unique(stations_close[stations_close$id==station,c('longitude')])
      latitude=unique(stations_close[stations_close$id==station,c('latitude')])
      
      
      d=matrix(c('station'=station,'longitude'=longitude,'latitude'=latitude,values),nrow = 1)
      colnames(d)<-c('station','longitude','latitude',dates)
      
      return(d)
    }
    
    
    tmaxdata=do.call('rbind',lapply(list_stations_possible_missingdates, function(x)download_data(x,datatypeid = 'TMAX')))
    tmindata=do.call('rbind',lapply(list_stations_possible_missingdates, function(x)download_data(x,datatypeid = 'TMIN')))
                     
  }
  
  
  ####Ordinary kriging####
  
  tmindata=as.data.frame(tmindata)
  tmaxdata=as.data.frame(tmaxdata)
  tmindata$longitude=as.numeric(tmindata$longitude)
  tmindata$latitude=as.numeric(tmindata$latitude)
  tmaxdata$longitude=as.numeric(tmaxdata$longitude)
  tmaxdata$latitude=as.numeric(tmaxdata$latitude)
  
  #tmindata
  sp::coordinates(tmindata)=c('longitude','latitude')
  proj4string(tmindata) = "+proj=longlat +datum=WGS84"
  
  
  k <- gstat(formula=OZDLYAV~1, locations=aq, model=fve)
  kp <- predict(k, g)
  
  
  if (length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function)))!=0){
    list_stations=unique(unlist(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))))
    
    
    station=list_stations[1]
    station_coords=stations_close[stations_close$id==station,c('longitude','latitude')]
    if (list_stations>1){
      station2= unique(list_stations[list_stations%notin%station][1])
      station2_coords=unique(stations_close[stations_close$id==station2,c('longitude','latitude')])
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
}
lapply(unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp)%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2')], safe_find_function)




















for (j in unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp)%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2')]) {
  
  
  
  station=vector()
  
  
  find_station=function(x){
    id_stations=stations_close[x,'id']
    values_tmax= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values_tmin= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values=data.frame(values_tmax,values_tmin)
    if(nrow(values)==length(yday(date_start):yday(date_end))){
      return(id_stations)
    }
    else{return(NULL)}
  }
  
  find_station2=function(x){
    id_stations=stations_close[x,'id']
    values_tmax= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values_tmin= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
    values=data.frame(values_tmax,values_tmin)
    if(nrow(values)>length(yday(date_start):yday(date_end))-5){
      return(id_stations)
    }
    else{return(NULL)}
  }
  
  safe_find_function<- function(x) {
    tryCatch(
      find_station(x),
      warning = function(w) NULL,
      error = function(e) NULL
    )
  }
  
  safe_find_function2<- function(x) {
    tryCatch(
      find_station2(x),
      warning = function(w) NULL,
      error = function(e) NULL
    )
  }
  
  if (length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function)))!=0){
    list_stations=unique(unlist(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))))
    
    station=list_stations[1]
    station_coords=stations_close[stations_close$id==station,c('longitude','latitude')]
    if (list_stations>1){
      station2= unique(list_stations[list_stations%notin%station][1])
      station2_coords=unique(stations_close[stations_close$id==station2,c('longitude','latitude')])
    }
    
  }
  
  
  
  if (length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function)))==0){
    list_stations=unique(unlist(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function2))))
    
    station=list_stations[1]
    station_coords=unique(stations_close[stations_close$id==station,c('longitude','latitude')])
    if (list_stations>1){
      station2= list_stations[list_stations%notin%station][1]
      station2_coords=unique(stations_close[stations_close$id==station2,c('longitude','latitude')])
      
    }
    
  }
  
  #dist=stations_close[stations_close$id==station,'distance']
  #dist2=stations_close[stations_close$id==station2,'distance']
  
  
  
  
  
  
  #Downloading the data for the growing season for two nearby stations
  ncdc_data_station1=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
  ncdc_data2_station1=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
  ncdc_data_station_dates=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$date
  ncdc_data_station_lonlat=station_coords
  values=data.frame( ncdc_data_station1, ncdc_data2_station1,ncdc_data_station_dates)
  download_data_station1<-data.frame('TMAX'=values$ncdc_data_station1,'TMIN'=values$ncdc_data2_station1,'YDAY'=yday(values$ncdc_data_station_dates),'month'=month(values$ncdc_data_station_dates))
  
  ncdc_data_station2=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station2,sep = ''),datatypeid='TMAX',startdate = date_start,enddate = date_end,limit = 500)$data$value
  ncdc_data2_station2=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station2,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$value
  ncdc_data_station2_dates=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station2,sep = ''),datatypeid='TMIN',startdate = date_start,enddate = date_end,limit = 500)$data$date
  values2=data.frame( ncdc_data_station2, ncdc_data2_station2,ncdc_data_station2_dates)
  download_data_station2<-data.frame('TMAX'=values2$ncdc_data_station2,'TMIN'=values2$ncdc_data2_station2,'YDAY'=yday(values2$ncdc_data_station2_dates),'month'=month(values2$ncdc_data_station2_dates))
  
  d=cbind(download_data_station1,download_data_station2)
  
  
  
  
  d$prec <- rowSums(d[, c(6:17)])
  
  dsp <- SpatialPoints(d[,4:3], proj4string=CRS("+proj=longlat +datum=NAD83"))
  dsp <- SpatialPointsDataFrame(dsp, d)
  
  
  
  gs<-gstat(formula=prec~1, locations=dta)
  
  
  for (s in 1:nrow(daily_weather[daily_weather$Year_Exp==j,])) {
    if (daily_weather[daily_weather$Year_Exp==j,'flagged_rain'][s]=='flagged'){
      day=as.numeric(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'][s])
      if (length(download_data_station1[download_data_station1$YDAY==day,'prcp'])!=0){
        daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'sum_rainfall']=download_data_station1[download_data_station1$YDAY==day,'prcp']/10
        daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'stationID_NOAA']=station
        daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'dist']=round(dist,2)
      }
      else{}
      
    }
  }
  
  data_plot=as.data.frame(cbind('month'=daily_weather[daily_weather$Year_Exp==j,'month'],'station_NCDC_1'=as.vector(download_data_station1[match(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'],download_data_station1$YDAY),'prcp'])/10,'station_NCDC_2'=as.vector(download_data_station2[match(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'],download_data_station2$YDAY),'prcp'])/10,'field_station'=daily_weather[daily_weather$Year_Exp==j,'sum_rainfall']))
  
  
  
  
  data_plot$month=as.factor(data_plot$month)
  data_sum <- plyr::ddply(data_plot, "month",
                          transform, total_station1_NOAA=sum(station_NCDC_1,na.rm = T))
  data_sum <- plyr::ddply(data_sum, "month",
                          transform, total_station2_NOAA=sum(station_NCDC_2,na.rm = T))
  data_sum <- plyr::ddply(data_sum, "month",
                          transform, total_station_field=sum(field_station))
  
  data_sum<-unique(data_sum[,c(1,5:7)])
  barplot(
    t(as.matrix(data_sum[, c(2, 3,4)])),
    beside = T,
    names.arg = data_sum$month,
    legend.text = T,
    ylim = c(0, max(data_sum[, c(2, 3)],na.rm = T) + 50),
    ylab = "Rainfall (mm)",
    xlab = "Month",
    cex.main=0.7,
    args.legend = list(x = 'topright', bty='n'),
    main = paste(j , '\n', 'Distance field to NOAA station 1: ', round(dist, 2),'km','\n','Distance field to NOAA station 2: ', round(dist2, 2),'km','\n','Nb flagged values: ',length(which(daily_weather[daily_weather$Year_Exp==j,'flagged_rain']=='flagged')), sep = '')
  )
  
  
  cc=cc+1
  
}
















