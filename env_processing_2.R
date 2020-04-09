##############################################################################
################################ WEATHER DATA PROCESSING #######################
##############################################################################
rm(list=ls())

library(rnoaa)
library(GSODR)
library(countyweather)
stations <- ghcnd_stations() 
stations<-filter(stations,last_year==2020)
stations<-filter(stations,first_year<=2014)

library(dplyr)
library(plyr)
library(lubridate)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")


weather = read.table(
  '~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

##CREATE A DAILY WEATHER DATASET


daily_weather=aggregate(weather[,18],by=list(weather$Day.of.Year, weather$Year_Exp),FUN=function(x)length(x))
colnames(daily_weather)=c('Day.of.Year','Year_Exp','daily_interval_measurements')


##Add other variables not related to the weather (lon,lat, coty, county, soil data, previous crop)
daily_weather$Year=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'Year']
daily_weather$Field.Location=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'Field.Location']
daily_weather$lat=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'lat']
daily_weather$long=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'long']
daily_weather$lonlat_added=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'lonlat_added']
daily_weather$City=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'City']
daily_weather$City_revgeo=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'city_revgeo']
daily_weather$state=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'state']
daily_weather$county=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'county']
daily_weather$zip=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'zip']
daily_weather$X..Silt=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'X..Silt']
daily_weather$X..Sand=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'X..Sand']
daily_weather$X..Clay=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'X..Clay']
daily_weather$soil_data_imputed=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'imputed']
daily_weather$Texture=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'Texture']
daily_weather$OM=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'OM']
daily_weather$Previous.crop=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'Previous.crop']
daily_weather$Station.ID=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'Station.ID']
daily_weather$NWS.Network=weather[match(daily_weather$Year_Exp,weather$Year_Exp),'NWS.Network']




##Add beginning and end of the growing season to weather table

## 2015_ILH2, 2015_IAH1, 2015_IAH2, 2015_IAH3, 2015_IAH4 removed because no phenotypic data present in the final hybrid pheno files

growingseason=read.table("~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/planting_harvest_dates/growingseason_extendeddates.txt",header = T,sep = '\t')
daily_weather$Date.Planted=growingseason[match(daily_weather$Year_Exp,growingseason$Year_Exp),2]
daily_weather$Date.Harvested=growingseason[match(daily_weather$Year_Exp,growingseason$Year_Exp),3]
daily_weather$Date.Planted2=as.POSIXct(as.Date(daily_weather$Date.Planted, "%m/%d/%Y"))
daily_weather$Date.Harvested2=as.POSIXct(as.Date(daily_weather$Date.Harvested, "%m/%d/%Y"))
daily_weather$Date.Planted3=NA
daily_weather$Date.Harvested3=NA

for (j in unique(daily_weather$Year)) {
  date=paste('01/01/',substring(j,3),sep = '')
  daily_weather[daily_weather$Year==j,'Date.Planted3']=as.numeric(difftime(daily_weather[daily_weather$Year==j,]$Date.Planted2,as.POSIXct(as.Date(date,"%m/%d/%Y")),units='days')+1)
  daily_weather[daily_weather$Year==j,'Date.Harvested3']=difftime(daily_weather[daily_weather$Year==j,]$Date.Harvested2,as.POSIXct(as.Date(date,"%m/%d/%Y")),units='days')+1
}
daily_weather=daily_weather[,-which(colnames(daily_weather)%in%c("Date.Planted","Date.Harvested","Date.Planted2","Date.Harvested2"))]
colnames(daily_weather)[c(23,24)]<-c("Date.Planted","Date.Harvested")

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

dd=list()
s=1
for (v in unique(daily_weather[!is.na(daily_weather$Date.Planted),'Year_Exp'])) {
  dd[[s]]=as.data.frame(matrix(NA,ncol = 2,nrow = length(unique(daily_weather[daily_weather$Year_Exp==v,"Date.Planted"]):unique(daily_weather[daily_weather$Year_Exp==v,"Date.Harvested"]))))
  dd[[s]][,2]=v
  dd[[s]][,1]=unique(daily_weather[daily_weather$Year_Exp==v,"Date.Planted"]):unique(daily_weather[daily_weather$Year_Exp==v,"Date.Harvested"])
  dd[[s]]=cbind(dd[[s]],rep.row(unique(daily_weather[daily_weather$Year_Exp==v,4:24]),nrow(dd[[s]])))
  s=s+1
}


df<-plyr::ldply(dd,data.frame)
colnames(df)<-colnames(daily_weather)[c(1:2,4:24)]
daily_weather=df
rm(df,dd)

#Add months and day
daily_weather$month <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%m'))
daily_weather$day <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%d'))



#########################################
####CHECK INSTANTANEOUS WEATHER DATA#####


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
  summarise(all(is.na(Rainfall..mm.)))
colnames(alldaymissing_rain)[3]='all_missing_rain'

detach('package:plyr')
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
  mutate(sum_rainfall=sum(Rainfall..mm.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,sum_rainfall)


#Add this daily rainfall sum to the daily_weather table
daily_weather<-merge(daily_weather,unique(j),by=c('Day.of.Year','Year_Exp'),all.x = T)


#Second Control range: daily values

daily_weather[which(daily_weather$sum_rainfall>=500),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall>=500),'sum_rainfall']=NA
daily_weather[which(daily_weather$sum_rainfall<0),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall<0),'sum_rainfall']=NA

daily_weather[is.na(daily_weather$flagged_rain),'flagged_rain']='flagged'


######IMPUTE RAIN MISSING VALUES #####



list_yearexp=list()
daily_weather$stationID_NOAA=NA
daily_weather$dist=NA
cc=1

for (j in unique(daily_weather$Year_Exp)) {
  print(j)
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'Year']))])
   
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==j,'lat']))])
  
  #Finding the closest stations in a radius of 30 km and select those with PRCP data
  stations_close=as.data.frame(meteo_distance(stations,latitude,longitude,radius = 30))
  stations_close<-filter(stations_close,element=='PRCP')
  #Choosing the closest to the field location
  dist=unique(stations_close[stations_close$distance==min(stations_close$distance),'distance'])
  
  id_stations=unique(stations_close[stations_close$distance==min(stations_close$distance),'id'])
  
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==j,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))

  
  #Downloading the data for the growing season
  download_data<-data.frame('prcp'=ncdc(datasetid = 'GHCND',stationid = paste('GHCND:',id_stations,sep = ''),datatypeid='PRCP',startdate = date_start,enddate = date_end,limit = 500)$data$value,'YDAY'=NA)
  download_data$YDAY<-yday(date_start):yday(date_end)
  
  
  
  for (s in 1:nrow(daily_weather[daily_weather$Year_Exp==j,])) {
    if (daily_weather[daily_weather$Year_Exp==j,'flagged_rain'][s]=='flagged'){
      day=as.numeric(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'][s])
      #print(download_data[download_data$YDAY==day,'prcp'])
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'sum_rainfall']=download_data[download_data$YDAY==day,'prcp']/10
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'stationID_NOAA']=id_stations
      daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'dist']=dist
    }
    if (daily_weather[daily_weather$Year_Exp==j,'flagged_rain'][s]=='OK'){
      day=as.numeric(daily_weather[daily_weather$Year_Exp==j,'Day.of.Year'][s])
      #print(download_data[download_data$YDAY==day,'prcp'])
      print(paste(day,daily_weather[daily_weather$Year_Exp==j&daily_weather$Day.of.Year==day,'sum_rainfall'],download_data[download_data$YDAY==day,'prcp']/10))
      
    }
  }
  
  #list_prcp[[cc]]=meteo_pull_monitors(monitors = list_yearexp[[cc]]$id,date_min = year,date_max = year,var='PRCP')
  #ut <- ncdc(datasetid='NORMAL_DLY', stationid='GHCND:USW00014895', datatypeid='dly-tmax-normal', startdate = '2010-05-01', enddate = '2010-05-10')
  
  cc=cc+1
}



####COMPARE MONTHLY DAILY SUMMARIES TO ENSURE RAINFALL DATA OK#######








list_counties=list()
#for (j in unique(daily_weather$zip)) {
  print(j)
  #list_counties[[j]]=daily_fips(fips=j,date_min='2014-01-01',
                                date_max='2018-12-15',var='prcp')
  
}














#The aim of the step test is to verify the rate of change of instantaneous data (detection of unrealistic jumps in values or dead band caused by blocked sensors). 

##TEMPERATURE
#Range test
weather[which(weather$Temperature..C.<(-30)),]<-NA
weather[which(weather$Temperature..C.>50),]<-NA

#Missing values per day should not be above 80%
aggregate()

aggregate(weather[,'Temperature..C.'],by=list(weather$Day.of.Year, weather$Year_Exp),FUN=function(x)length(which(is.na(x))))[,3]


#Step test vary according to the measurements interval: at most stations: semi-hourly, but not always the case
if (daily_weather$daily_interval_measurements==48)
  
  
  
  daily_weather$missing_obs_temperature=aggregate(weather[,18],by=list(weather$Day.of.Year, weather$Year_Exp),FUN=function(x)length(which(is.na(x))))[,3]


#weather$growth_TEMP <- with(weather, ave(Temperature..C., Year_Exp,
#                          FUN=function(x) c(NA, diff(x) / tail(x, -1))))
weather$growth_TEMP <- with(weather, ave(Temperature..C., Year_Exp,
                                         FUN=function(x) c(NA, diff(x) )))








colnames(daily_weather)=c('Day.of.Year','Year_Exp','mean.temperature')
daily_weather$=NA
