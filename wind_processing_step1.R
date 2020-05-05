rm(list = ls())
# ------------------------------------------------------------------------------
#Load packages
# ------------------------------------------------------------------------------
library(rnoaa)
library(GSODR)
library(gstat)
library(Rcpp)
library(raster)
library(sp)
library(mapdata)
library(maps)
library(maptools)
library(xts)
library(spacetime)
library(rgdal)

`%notin%` <- Negate(`%in%`)


library(dplyr)
#library(plyr)
library(lubridate)

# ------------------------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------------------------
weather = read.table(
  'weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

daily_weather = read.table(
  'daily_weather.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
# ------------------------------------------------------------------------------
##HUMIDITY: flagged values to be imputed
# ------------------------------------------------------------------------------
##Control range semi-hourly or 15-, 20-minutes data values
weather$flagged_humidity_instant=NA
weather[which(weather$Relative.Humidity....>100),'flagged_humidity_instant']<-'flagged'
weather[which(weather$Relative.Humidity....>100),'Relative.Humidity....']<-NA
weather[which(weather$Relative.Humidity....<0),'flagged_humidity_instant']<-'flagged'
weather[which(weather$Relative.Humidity....<0),'Relative.Humidity....']<-NA

#First flagged value assignment based on ratio (humidity values present)/(nb obs per day) >0.9 
weather$Year_Exp=as.factor(weather$Year_Exp)
weather$Day.of.Year=as.factor(weather$Day.of.Year)
total_obs_daily<-as.data.frame(weather %>% group_by(Year_Exp,Day.of.Year) %>% tally())

nomissing_humidity<-weather %>% 
  aggregate(Relative.Humidity....~Day.of.Year+Year_Exp,data=.,FUN = (function(x) {sum(!is.na(x))})) 
colnames(nomissing_humidity)[3]='number_obs_present'

alldaymissing_humidity<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all(is.na(Relative.Humidity....)))
colnames(alldaymissing_humidity)[3]='all_missing_humidity'



humidity<-merge(total_obs_daily,nomissing_humidity,by=c('Day.of.Year','Year_Exp'),all.x=T)
humidity<-merge(humidity,alldaymissing_humidity,by=c('Day.of.Year','Year_Exp'),all.x=T)

rm(total_obs_daily,nomissing_humidity,alldaymissing_humidity)

humidity<-humidity%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=number_obs_present/n)%>%mutate(flagged_humidity=case_when(ratio<0.9 ~'flagged',
                                                                                                                             ratio>=0.9 ~'OK'))
# ------------------------------------------------------------------------------
#Weird number of daily observations recorded (ex: 50, 100 etc and not 24, 48, 72) then the value is flagged too.
# ------------------------------------------------------------------------------

humidity[humidity$n<24,'flagged_humidity']='flagged'
humidity[humidity$n>96,'flagged_humidity']='flagged'
humidity[humidity$n>24&humidity$n<48,'flagged_humidity']='flagged'
humidity[humidity$n>48&humidity$n<72,'flagged_humidity']='flagged'
humidity[humidity$n>72&humidity$n<96,'flagged_humidity']='flagged'
humidity[humidity$all_missing_humidity==TRUE,'flagged_humidity']='flagged'
humidity[humidity$all_missing_humidity==TRUE,'ratio']=0

weather=merge(weather,humidity[,c(1:2,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)

daily_weather=merge(daily_weather,humidity[,c(1:3,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)


# ------------------------------------------------------------------------------
#Daily mean,max and min humidity computed for unflagged values:
# ------------------------------------------------------------------------------

maxH<-weather%>%
  filter(flagged_humidity%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(HMAX=max(Relative.Humidity....,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,HMAX)
maxH<-unique(maxH)
minH<-weather%>%
  filter(flagged_humidity%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(HMIN=min(Relative.Humidity....,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,HMIN)
minH<-unique(minH)
allhumidity=merge(maxH,minH,by=c('Day.of.Year','Year_Exp'),all.x=T)
allhumidity$HMEAN=(allhumidity$HMAX+allhumidity$HMIN)/2
allhumidity=arrange(allhumidity,Year_Exp,Day.of.Year)

# ------------------------------------------------------------------------------
#Add not flagged daily computed HMIN, HMAX, HMEAN to the daily_weather table (based on data from the field station)
# ------------------------------------------------------------------------------
daily_weather<-merge(daily_weather,allhumidity,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)






##Note (to add): The preferred method to obtain an average RH over a given interval is to average the vapor pressure and saturation vapor pressures, then calculate RH from those averages




#Write the table which will be used for comparison with interpolated values for non missing values
write.table(daily_weather,'daily_weather_humidity_processed1.txt',col.names=T,row.names=F,sep='\t',quote=F)

