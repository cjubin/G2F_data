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
##Radiation: flagged values to be imputed
# ------------------------------------------------------------------------------
##Control range semi-hourly or 15-, 20-minutes data values
weather$flagged_radiation_instant=NA
weather$Solar.Radiation..W.m2.=as.numeric(weather$Solar.Radiation..W.m2.)
weather[which(weather$Solar.Radiation..W.m2.>1500),'flagged_radiation_instant']<-'flagged'
weather[which(weather$Solar.Radiation..W.m2.>1500),'Solar.Radiation..W.m2.']<-NA
weather[which(weather$Solar.Radiation..W.m2.<0),'flagged_radiation_instant']<-'flagged'
weather[which(weather$Solar.Radiation..W.m2.<0),'Solar.Radiation..W.m2.']<-NA



#First flagged value assignment based on ratio (solar rad values present)/(nb obs per day) >0.9 
weather$Year_Exp=as.factor(weather$Year_Exp)
weather$Day.of.Year=as.factor(weather$Day.of.Year)
total_obs_daily<-as.data.frame(weather %>% group_by(Year_Exp,Day.of.Year) %>% tally())

nomissing_solarrad<-weather %>% 
  aggregate(Solar.Radiation..W.m2.~Day.of.Year+Year_Exp,data=.,FUN = (function(x) {sum(!is.na(x))})) 
colnames(nomissing_solarrad)[3]='number_obs_present'

alldaymissing_solarrad<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all(is.na(Solar.Radiation..W.m2.)))
colnames(alldaymissing_solarrad)[3]='all_missing_solarrad'



solarrad<-merge(total_obs_daily,nomissing_solarrad,by=c('Day.of.Year','Year_Exp'),all.x=T)
solarrad<-merge(solarrad,alldaymissing_solarrad,by=c('Day.of.Year','Year_Exp'),all.x=T)

rm(total_obs_daily,nomissing_solarrad,alldaymissing_solarrad)

solarrad<-solarrad%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=number_obs_present/n)%>%mutate(flagged_solarrad=case_when(ratio<0.9 ~'flagged',
                                                                                                                 ratio>=0.9 ~'OK'))


# ------------------------------------------------------------------------------
#Weird number of daily observations recorded (ex: 50, 100 etc and not 24, 48, 72) then the value is flagged too.
# ------------------------------------------------------------------------------

solarrad[solarrad$n<24,'flagged_solarrad']='flagged'
solarrad[solarrad$n>96,'flagged_solarrad']='flagged'
solarrad[solarrad$n>24&solarrad$n<48,'flagged_solarrad']='flagged'
solarrad[solarrad$n>48&solarrad$n<72,'flagged_solarrad']='flagged'
solarrad[solarrad$n>72&solarrad$n<96,'flagged_solarrad']='flagged'
solarrad[solarrad$all_missing_solarrad==TRUE,'flagged_solarrad']='flagged'
solarrad[solarrad$all_missing_solarrad==TRUE,'ratio']=0

weather=merge(weather,solarrad[,c(1:2,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)

daily_weather=merge(daily_weather,solarrad[,c(1:3,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)



# ------------------------------------------------------------------------------
#Daily total solar rad speed computed for unflagged values:
# ------------------------------------------------------------------------------

total_solarrad<-weather%>%
  filter(flagged_solarrad%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(incoming_radiation=sum(Solar.Radiation..W.m2.,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,incoming_radiation)
total_solarrad<-unique(total_solarrad)

total_solarrad=arrange(total_solarrad,Year_Exp,Day.of.Year)

# ------------------------------------------------------------------------------
#Step test: absolute difference between two consecutive days should be different from 0:
# ------------------------------------------------------------------------------

library(dplyr)
test  <- 
  total_solarrad %>%
  group_by(Year_Exp) %>%
  mutate(lag = dplyr::lag(incoming_radiation, n = 1, default = NA))%>%
  mutate(diff = abs(incoming_radiation-lag))%>%
  filter(diff>0)%>%
  dplyr::select(Year_Exp,Day.of.Year,incoming_radiation)

total_solarrad<-test

# ------------------------------------------------------------------------------
#Add not flagged daily computed incoming_radiation to the daily_weather table (based on data from the field station)
# ------------------------------------------------------------------------------
daily_weather<-merge(daily_weather,total_solarrad,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)

# ------------------------------------------------------------------------------
#Conversion from  watts per square meter perday to megajoules per square meter per day 
# ------------------------------------------------------------------------------
daily_weather$incoming_radiation_MJm2<- daily_weather$incoming_radiation* 0.0864/daily_weather$n


#Write the table which will be used for comparison with interpolated values for non missing values
write.table(daily_weather,'daily_weather_solarrad_processed1.txt',col.names=T,row.names=F,sep='\t',quote=F)

