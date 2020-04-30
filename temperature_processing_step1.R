rm(list=ls())
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



weather = read.table(
  'weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

daily_weather=read.table(
  'daily_weather.txt',
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
  dplyr::mutate(TMAX=max(Temperature..C.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,TMAX)
maxT<-unique(maxT)
minT<-weather%>%
  filter(flagged_temp%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(TMIN=min(Temperature..C.,na.rm = T))%>%
  select(Day.of.Year,Year_Exp,TMIN)
minT<-unique(minT)
temperatures=merge(maxT,minT,by=c('Day.of.Year','Year_Exp'),all.x=T)

temperatures=arrange(temperatures,Year_Exp,Day.of.Year)



#Internal consistency test: Tmax(d) > Tmin(d-1) + Tmin(d) < Tmax(d-1)
library(data.table)
nm1 <- c('TMAX','TMIN')
nm2 <- paste("lag", nm1, sep=".")

library(dplyr)
test  <- 
  temperatures %>%
  group_by(Year_Exp) %>%
  mutate(lag = dplyr::lag(TMIN, n = 1, default = NA))%>%
  mutate(diff = TMAX-lag)%>%
  mutate(lag2 = dplyr::lag(TMAX, n = 1, default = NA))%>%
  mutate(diff2 = TMIN-lag2)%>%
  filter(diff>0)%>%
  filter(diff2<0)%>%
  select(Year_Exp,Day.of.Year,TMIN,TMAX)

temperatures<-test


#Mean temperature: WMO (2010) recommends use of this estimator:'Even though this method is not the best statistical approximation, its consistent use satisfies the comparative purpose of normal'

temperatures$TMEAN=(temperatures$TMAX+temperatures$TMIN)/2



#Add not flagged daily computed TMIN,TMAX, TMEAN to the daily_weather table (based on data from the field station)
daily_weather<-merge(daily_weather,temperatures,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)


#Write the table which will be used for comparison with interpolated values for non missing values
write.table(daily_weather,'daily_weather_temp_processed1.txt',col.names=T,row.names=F,sep='\t',quote=F)

