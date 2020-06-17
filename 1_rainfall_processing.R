rm(list=ls())


library(dplyr)
library(lubridate)

options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")


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



#Daily sum computed for not flagged values:

j<-weather%>%
  filter(flagged_rain%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(sum_rainfall=sum(Rainfall..mm.,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,sum_rainfall)


#Add this daily rainfall sum to the daily_weather table
daily_weather<-merge(daily_weather,unique(j),by=c('Day.of.Year','Year_Exp'),all.x = T)


#Second Control range: daily values

daily_weather[which(daily_weather$sum_rainfall>=160),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall>=160),'sum_rainfall']=NA
daily_weather[which(daily_weather$sum_rainfall<0),'flagged_rain']='flagged'
daily_weather[which(daily_weather$sum_rainfall<0),'sum_rainfall']=NA

daily_weather[is.na(daily_weather$flagged_rain),'flagged_rain']='flagged'

daily_weather=arrange(daily_weather,Year_Exp,Day.of.Year)


write.table(daily_weather,'daily_weather_prcp_processed1.txt',col.names = T,row.names = F,sep = '\t',quote = F)

