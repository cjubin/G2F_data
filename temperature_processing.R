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
daily_weather=merge(daily_weather,already_imputed,by=c('Day.of.Year','Year_Exp'),all.x=T)

#Minimum/ maximum  allowed  variability  of  an  instantaneous  value: step test and persistence test

<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%



#Daily mean,max and min computed for not flagged values:


maxT<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(max_temp=max(Temperature..C,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,max_temp)
minT<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(min_temp=min(Temperature..C,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,min_temp)




#Add this daily temp sum to the daily_weather table
daily_weather<-merge(daily_weather,unique(j),by=c('Day.of.Year','Year_Exp'),all.x = T)


#Second Control range: daily values



######IMPUTE temp MISSING VALUES #####

daily_weather=arrange(daily_weather,Year,Field.Location,Day.of.Year)

daily_weather$stationID_NOAA=NA
daily_weather$dist=NA
cc=1
par(mfrow=c(2,2))
`%notin%` <- Negate(`%in%`)

for (j in unique(daily_weather$Year_Exp)[unique(daily_weather$Year_Exp)%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2')]) {
  print(j)
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'lat']))])
  
  #Finding the closest stations in a radius of 30 km and select those with PRCP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = 50))
  stations_close<-filter(stations_close,element=='TEMP')
  
  #Choosing the closest to the field location with all the dates available
  stations_close<-arrange(stations_close,distance)
  station=vector()
  for (l in 1:nrow(stations_close)) {
    id_stations=stations_close[l,'id']
    values= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data$value
    if(length(values)==length(yday(date_start):yday(date_end))){
      station= id_stations
      break}
  }
  
  if(length(station)==0){
    for (l in 1:nrow(stations_close)) {
      id_stations=stations_close[l,'id']
      values= ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data$value
      if(length(values)>length(yday(date_start):yday(date_end))-5){
        station= id_stations
        break}
    }
  }
  
  dist=stations_close[stations_close$id==station,'distance']
  
  daily_weather[daily_weather$Year_Exp==j,'stationID_NOAA']=station
  daily_weather[daily_weather$Year_Exp==j,'dist']=round(dist,2)
  
  
  
  
  #Downloading the data for the growing season
  ncdc_data=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',station,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data
  download_data<-data.frame('prcp'=ncdc_data$value,'YDAY'=yday(ncdc_data$date),'month'=month(ncdc_data$date))
  
  
  
  for (s in 1:nrow(daily_weather[daily_weather$Year_Exp==j,])) {
    if (daily_weather[daily_weather$Year_Exp==j,'flagged_temp'][s]=='flagged'){
      day=as.numeric(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'][s])
      #print(download_data[download_data$YDAY==day,'prcp'])
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'sum_temp']=download_data[download_data$YDAY==day,'prcp']/10
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'stationID_NOAA']=id_stations
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'dist']=dist
    }
  }
  
  
  
  data_plot=as.data.frame(cbind('month'=daily_weather[daily_weather$Year_Exp==j,'month'],'station_NCDC'=as.vector(download_data[match(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'],download_data$YDAY),'prcp'])/10,'field_station'=daily_weather[daily_weather$Year_Exp==j,'sum_temp']))
  
  data_plot$month=as.factor(data_plot$month)
  
  data_sum <- plyr::ddply(data_plot, "month",
                          transform, total_station_NOAA=sum(station_NCDC))
  data_sum <- plyr::ddply(data_sum, "month",
                          transform, total_station_field=sum(field_station))
  data_sum<-unique(data_sum[,c(1,4,5)])
  barplot(
    t(as.matrix(data_sum[, c(2, 3)])),
    beside = T,
    names.arg = data_sum$month,
    legend.text = T,
    ylim = c(0, max(data_sum[, c(2, 3)],na.rm = T) + 50),
    ylab = "temp (mm)",
    xlab = "Month",
    main = paste(j , '\n', 'Distance field to NOAA station: ', round(dist, 2),'km','\n','Nb flagged values: ',length(which(daily_weather[daily_weather$Year_Exp==j,'flagged_temp']=='flagged')), sep = '')
  )
  
  
  cc=cc+1
}

















#Step test vary according to the measurements interval: at most stations: semi-hourly, but not always the case
if (daily_weather$daily_interval_measurements==48)
  
  
  
  daily_weather$missing_obs_temperature=aggregate(weather[,18],by=list(weather$Day.of.Year, weather$Year_Exp),FUN=function(x)length(which(is.na(x))))[,3]


#weather$growth_TEMP <- with(weather, ave(Temperature..C., Year_Exp,
#                          FUN=function(x) c(NA, diff(x) / tail(x, -1))))
weather$growth_TEMP <- with(weather, ave(Temperature..C., Year_Exp,
                                         FUN=function(x) c(NA, diff(x) )))


