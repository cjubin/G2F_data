rm(list=ls())
library(rnoaa)
library(GSODR)
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

##RAINFALL: flagged values to be imputed

##Control range semi-hourly or 15-, 20-minutes data values
weather$flagged_rain_instant=NA
weather[which(weather$Rainfall..mm.>120),'flagged_rain_instant']<-'flagged'
weather[which(weather$Rainfall..mm.>120),'Rainfall..mm.']<-NA
weather[which(weather$Rainfall..mm.<0),'flagged_rain_instant']<-'flagged'
weather[which(weather$Rainfall..mm.<0),'Rainfall..mm.']<-NA

#First flagged value assignment based on ratio rain values present/nb obs per day + if weird number of daily obs (ex: 50, 100 etc) then the value is flagged too.
weather$Year_Exp=as.factor(weather$Year_Exp)
weather$Day.of.Year=as.factor(weather$Day.of.Year)
total_obs_daily<-as.data.frame(weather %>% group_by(Year_Exp,Day.of.Year) %>% tally())

nomissing_rain<-weather %>% 
  aggregate(Rainfall..mm.~Day.of.Year+Year_Exp,data=.,FUN = (function(x) {sum(!is.na(x))})) 
colnames(nomissing_rain)[3]='number_obs_present'

alldaymissing_rain<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all(is.na(Rainfall..mm.)))
colnames(alldaymissing_rain)[3]='all_missing_rain'

#detach('package:plyr')
already_imputed<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  mutate(already_imputed=case_when(Cleaning.Method%in%'Imputed' ~'already_imputed'))%>%
  select(Day.of.Year,Year_Exp,already_imputed)
already_imputed<-unique(already_imputed)


rain<-merge(total_obs_daily,nomissing_rain,by=c('Day.of.Year','Year_Exp'),all.x=T)
rain<-merge(rain,alldaymissing_rain,by=c('Day.of.Year','Year_Exp'),all.x=T)

rm(total_obs_daily,nomissing_rain,alldaymissing_rain)

rain<-rain%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=number_obs_present/n)%>%mutate(flagged_rain=case_when(ratio<0.9 ~'flagged',
                                                                                                                 ratio>=0.9 ~'OK'))

rain[rain$n<24,'flagged_rain']='flagged'
rain[rain$n>96,'flagged_rain']='flagged'
rain[rain$n>24&rain$n<48,'flagged_rain']='flagged'
rain[rain$n>48&rain$n<72,'flagged_rain']='flagged'
rain[rain$n>72&rain$n<96,'flagged_rain']='flagged'
rain[rain$all_missing_rain==TRUE,'flagged_rain']='flagged'
rain[rain$all_missing_rain==TRUE,'ratio']=0

weather=merge(weather,rain[,c(1:2,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)

daily_weather=merge(daily_weather,rain[,c(1:3,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)
daily_weather=merge(daily_weather,already_imputed,by=c('Day.of.Year','Year_Exp'),all.x=T)



#Daily sum computed for not flagged values:

j<-weather%>%
  filter(flagged_rain%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(sum_rainfall=sum(Rainfall..mm.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,sum_rainfall)


#Add this daily rainfall sum to the daily_weather table
daily_weather<-merge(daily_weather,unique(j),by=c('Day.of.Year','Year_Exp'),all.x = T)


#Second Control range: daily values

daily_weather[which(daily_weather$sum_rainfall>=500),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall>=500),'sum_rainfall']=NA
daily_weather[which(daily_weather$sum_rainfall<0),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall<0),'sum_rainfall']=NA

daily_weather[is.na(daily_weather$flagged_rain),'flagged_rain']='flagged'
colnames(daily_weather)[29]<-'PRCP'

write.table(daily_weather,'daily_weather_prcp_processed1.txt',col.names = T,row.names = F,sep = '\t',quote = F)

######IMPUTE RAIN MISSING VALUES #####

daily_weather=arrange(daily_weather,Year,Field.Location,Day.of.Year)

daily_weather$stationID_NOAA=NA
daily_weather$dist=NA
cc=1
par(mfrow=c(2,1))
`%notin%` <- Negate(`%in%`)

for (j in unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp)%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2')]) {
  print(j)
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'lat']))])
  
  #Finding the closest stations in a radius of 50 km and select those with PRCP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = 50))
  stations_close<-filter(stations_close,element=='PRCP')
  
  #Choosing the closest to the field location with all the dates available
  stations_close<-arrange(stations_close,distance)
  station=vector()
  
  
  find_station=function(x){
    id_stations=stations_close[x,'id']
    values= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data$value
    if(length(values)==length(yday(date_start):yday(date_end))){
      return(id_stations)
    }
    else{return(NULL)}
  }
  
  find_station2=function(x){
    id_stations=stations_close[x,'id']
    values= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data$value
    if(length(values)>length(yday(date_start):yday(date_end))-5){
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
    station=Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))[[1]]
    if (length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))[[2]])!=0){
      station2=Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function))[[2]]
    }
    
  }
  
  

  if(length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function)))==0){
    station=Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function2))[[1]]
    if (length(Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function2))[[2]])!=0){
      station2=Filter(function(x) !is.null(x), x = lapply(1:nrow(stations_close), safe_find_function2))[[2]]
    }
  }
  
  dist=stations_close[stations_close$id==station,'distance']
  dist2=stations_close[stations_close$id==station2,'distance']

  
  
  
  
  
  #Downloading the data for the growing season
  ncdc_data_station1=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data
  download_data_station1<-data.frame('prcp'=ncdc_data_station1$value,'YDAY'=yday(ncdc_data_station1$date),'month'=month(ncdc_data_station1$date))
  
  ncdc_data_station2=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station2,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data
  download_data_station2<-data.frame('prcp'=ncdc_data_station2$value,'YDAY'=yday(ncdc_data_station2$date),'month'=month(ncdc_data_station2$date))
  
  
  
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







