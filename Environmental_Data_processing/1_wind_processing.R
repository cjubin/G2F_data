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
library(plyr)
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
##WIND: flagged values to be imputed
# ------------------------------------------------------------------------------
##Control range semi-hourly or 15-, 20-minutes data values
##Plausible value check according to Guidelines on Quality Control Procedures for Data from Automatic Weather Stations  (WORLD METEOROLOGICAL ORGANIZATION)
weather$flagged_WIND_instant=NA
weather$Wind.Speed..m.s.=as.numeric(weather$Wind.Speed..m.s.)
weather[which(weather$Wind.Speed..m.s.>75),'flagged_WIND_instant']<-'flagged'
weather[which(weather$Wind.Speed..m.s.>75),'Wind.Speed..m.s.']<-NA
weather[which(weather$Wind.Speed..m.s.<0),'flagged_WIND_instant']<-'flagged'
weather[which(weather$Wind.Speed..m.s.<0),'Wind.Speed..m.s.']<-NA

# ------------------------------------------------------------------------------
#Internal consistency test 1)wind speed = 00 and wind direction = 00 and 2)wind speed ≠ 00 and wind direction ≠ 00:
# ------------------------------------------------------------------------------
weather[which(weather$Wind.Direction..degrees.==0&weather$Wind.Speed..m.s.!=0),'flagged_WIND_instant']<-'flagged'
weather[which(weather$Wind.Direction..degrees.==0&weather$Wind.Speed..m.s.!=0),'Wind.Speed..m.s.']<-NA

weather[which(weather$Wind.Direction..degrees.!=0&weather$Wind.Speed..m.s.==0),'flagged_WIND_instant']<-'flagged'
weather[which(weather$Wind.Direction..degrees.!=0&weather$Wind.Speed..m.s.==0),'Wind.Speed..m.s.']<-NA



#First flagged value assignment based on ratio (WIND values present)/(nb obs per day) >0.9 
weather$Year_Exp=as.factor(weather$Year_Exp)
weather$Day.of.Year=as.factor(weather$Day.of.Year)
total_obs_daily<-as.data.frame(weather %>% group_by(Year_Exp,Day.of.Year) %>% tally())

nomissing_WIND<-weather %>% 
  aggregate(Wind.Speed..m.s.~Day.of.Year+Year_Exp,data=.,FUN = (function(x) {sum(!is.na(x))})) 
colnames(nomissing_WIND)[3]='number_obs_present'

alldaymissing_WIND<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all(is.na(Wind.Speed..m.s.)))
colnames(alldaymissing_WIND)[3]='all_missing_WIND'



WIND<-merge(total_obs_daily,nomissing_WIND,by=c('Day.of.Year','Year_Exp'),all.x=T)
WIND<-merge(WIND,alldaymissing_WIND,by=c('Day.of.Year','Year_Exp'),all.x=T)

rm(total_obs_daily,nomissing_WIND,alldaymissing_WIND)

WIND<-WIND%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=number_obs_present/n)%>%mutate(flagged_WIND=case_when(ratio<0.9 ~'flagged',
                                                                                                                             ratio>=0.9 ~'OK'))


# ------------------------------------------------------------------------------
#Weird number of daily observations recorded (ex: 50, 100 etc and not 24, 48, 72) then the value is flagged too.
# ------------------------------------------------------------------------------

WIND[WIND$n<24,'flagged_WIND']='flagged'
WIND[WIND$n>96,'flagged_WIND']='flagged'
WIND[WIND$n>24&WIND$n<48,'flagged_WIND']='flagged'
WIND[WIND$n>48&WIND$n<72,'flagged_WIND']='flagged'
WIND[WIND$n>72&WIND$n<96,'flagged_WIND']='flagged'
WIND[WIND$all_missing_WIND==TRUE,'flagged_WIND']='flagged'
WIND[WIND$all_missing_WIND==TRUE,'ratio']=0

weather=merge(weather,WIND[,c(1:2,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)

daily_weather=merge(daily_weather,WIND[,c(1:3,7)],by=c('Day.of.Year','Year_Exp'),all.x=T)



# ------------------------------------------------------------------------------
#Daily average wind speed computed for unflagged values:
# ------------------------------------------------------------------------------

mean_wind<-weather%>%
  filter(flagged_WIND%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(MEANWINDSPEED=mean(Wind.Speed..m.s.,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,MEANWINDSPEED)
mean_wind<-unique(mean_wind)

mean_wind=arrange(mean_wind,Year_Exp,Day.of.Year)

daily_weather<-merge(daily_weather,mean_wind,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)

# ------------------------------------------------------------------------------
# Step test: absolute difference between two consecutive days should be <10 m/s.
# Persistence test: daily mean wind speed should not be strictly equal to the value of the day before or to the day even before:
# U(d) != U(d-1) != U(d-2);
# ------------------------------------------------------------------------------


test  <- 
  daily_weather[,c('Year_Exp','Day.of.Year','flagged_WIND','MEANWINDSPEED')] %>%
  group_by(Year_Exp) %>%
  mutate(lag_1 = dplyr::lag(MEANWINDSPEED, n = 1, default = NA))%>%
  mutate(lag_2 = dplyr::lag(MEANWINDSPEED, n = 2, default = NA))%>%
  mutate(diff_1 = abs(MEANWINDSPEED-lag_1))%>%
  mutate(diff_2 = abs(MEANWINDSPEED-lag_2))%>%
  filter(diff_1<10&diff_1!=0&diff_2!=0)%>%
  mutate(flagged_WIND = 'OK')%>%
  dplyr::select(Year_Exp,Day.of.Year,MEANWINDSPEED,flagged_WIND)

mean_wind<-test





# ------------------------------------------------------------------------------
#Add not flagged daily computed MEANWINDSPEED to the daily_weather table (based on data from the field station)
# ------------------------------------------------------------------------------
daily_weather$MEANWINDSPEED<-NA
daily_weather<-merge(daily_weather,mean_wind,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)
daily_weather=daily_weather[,c(1:26,29,30)]



#Write the table which will be used for comparison with interpolated values for non missing values
write.table(daily_weather,'daily_weather_wind_processed1.txt',col.names=T,row.names=F,sep='\t',quote=F)

