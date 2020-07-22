rm(list=ls())


library(dplyr)
library(data.table)
library(plyr)
library(lubridate)
`%notin%` <- Negate(`%in%`)
options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")


weather = read.table(
  'weather_intermediate.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

prcp_aggregate_raw=aggregate(Rainfall..mm. ~ Year_Exp,weather, sum)
write.table(prcp_aggregate_raw,'prcp_aggregate_raw.txt',col.names = T,row.names = F,sep = '\t',quote = F)

daily_weather=read.table(
  'daily_weather_afterstep0.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA,''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

##RAINFALL: flagged values to be imputed

##Control range semi-hourly or 15-, 20-minutes data values

weather[which(weather$Rainfall..mm.>120),'Rainfall..mm.']<-NA

weather[which(weather$Rainfall..mm.<0),'Rainfall..mm.']<-NA

#First flagged value assignment based on ratio rain values present/nb obs per day + if weird number of daily obs (ex: 50, 100 etc) then the value is flagged too.
total_obs_daily<-weather %>% group_by(Year_Exp,Day.of.Year) %>% tally()
total_obs_daily$Year_Exp=as.factor(total_obs_daily$Year_Exp)
total_obs_daily<-data.table(total_obs_daily)
m1<-total_obs_daily[, .(number_of_distinct_n = uniqueN(n)), by = Year_Exp]
m2<-setDT(total_obs_daily)[, .(Freq = .N), by = .(Year_Exp, n)]

#Identification of cases where the number of daily obs is not round (24,48,72,96...) but still usable because constant during growing season
subset<-m2[m2$Freq>60,c('Year_Exp','n')]
#Only the case for 2015_WIH2

## Stations for which too much irregularities in the total number of semi-hourly observations per day (i.e. when the number of distinct n (n=number of semi-hourly obs. per day) per Year_Exp is superior to 12)
subset2<-m1[m1$number_of_distinct_n>18,c('Year_Exp')]

## Get the number of missing values per day, as well as a logical vector for days where all values are missing (values automatically assigned as flagged)

nomissing_rain<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(sum_notNA=sum(!is.na(Rainfall..mm.)))
 

alldaymissing_rain<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all_missing_rain=all(is.na(Rainfall..mm.)))
nomissing_rain<-merge(nomissing_rain,alldaymissing_rain,by=c('Year_Exp','Day.of.Year'))

total_obs_daily<-merge(total_obs_daily,nomissing_rain,by=c('Year_Exp','Day.of.Year'))



## Flagged Day.of.Year x Year_Exp values based on:

# 1) all records for that day are missing --> flagged

total_obs_daily<-total_obs_daily%>%group_by(Day.of.Year,Year_Exp)%>%mutate(flagged_rain=case_when(all_missing_rain==TRUE ~'flagged'))

# 2) 90% of the records for that day are not NA --> OK

total_obs_daily<-total_obs_daily%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=sum_notNA/n)%>%mutate(flagged_rain=case_when(ratio<0.9 ~'flagged',
                                                                                                                 ratio>=0.9 ~'OK'))

# 3) Less than 23 observations per day, even if no NA present, should be removed.

total_obs_daily[total_obs_daily$n<23,'flagged_rain']='flagged'

# 4) Except 2018_WIH1 and 2015_WIH2 for which many n==47 (manually checked and actually OK),
# we apply the following filters on n: when n is not a round value, something usually went wrong and data are not certain (unobserved values at unknown periods of the day for instance). Example For 2018_OHH1, from Day.of.Year 151 only round values for temp! Uncertainty therefore about the rest of the weather variables 
# + Some Year locations (identified in subset 2 above) very difficult to tell the reliability because of too many irregularities in the number of daily records. For these, we prefer to exclude the weather records.

total_obs_daily[total_obs_daily$n>96,'flagged_rain']='flagged'
total_obs_daily[total_obs_daily$n>24&total_obs_daily$n<48,'flagged_rain']='flagged'
total_obs_daily[total_obs_daily$n>48&total_obs_daily$n<72,'flagged_rain']='flagged'
total_obs_daily[total_obs_daily$n>72&total_obs_daily$n<96,'flagged_rain']='flagged'
total_obs_daily[total_obs_daily$Year_Exp%in%c('2018_WIH1','2015_WIH2')&total_obs_daily$n==47,'flagged_rain']='OK'
total_obs_daily[total_obs_daily$Year_Exp%in%subset2$Year_Exp,'flagged_rain']='flagged'

weather=merge(weather,total_obs_daily[,c('Day.of.Year','Year_Exp','flagged_rain')],by=c('Day.of.Year','Year_Exp'),all.x=T)


write.table(total_obs_daily,'flagged_prcp.txt',col.names = T,row.names = F,sep = '\t',quote = F)



#Daily sum computed for not flagged values:

j<-weather%>%
  filter(flagged_rain%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(sum_rainfall=sum(Rainfall..mm.,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,sum_rainfall)
j<-unique(j)

#Add this daily rainfall sum to the daily_weather table
daily_weather<-merge(daily_weather,j,by=c('Day.of.Year','Year_Exp'),all.x = T)


#Second Control range: daily values

daily_weather[which(daily_weather$sum_rainfall>=140),'sum_rainfall']=NA
daily_weather[which(daily_weather$sum_rainfall<0),'sum_rainfall']=NA


daily_weather=arrange(daily_weather,Year_Exp,Day.of.Year)
prcp_aggregate_processed=aggregate(sum_rainfall ~ Year_Exp,daily_weather, sum)
prcp_aggregate_processed<-merge(prcp_aggregate_processed,prcp_aggregate_raw,by='Year_Exp',all.y = T)
colnames(prcp_aggregate_processed)[2:3]<-c('after_processing_and_reduced_to_growingseason_length','before_any_processing')

write.table(daily_weather,'daily_weather_prcp_processed1.txt',col.names = T,row.names = F,sep = '\t',quote = F)

