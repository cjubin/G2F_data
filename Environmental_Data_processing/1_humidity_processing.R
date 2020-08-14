rm(list = ls())
# ------------------------------------------------------------------------------
#Load packages
# ------------------------------------------------------------------------------


packages = c("xts","data.table","rgdal","readxl","opencage","intrval","dplyr","USAboundaries","ggplot2","rnoaa","stringr","ggmap","revgeo","lubridate")

## Load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

`%notin%` <- Negate(`%in%`)



# ------------------------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------------------------
weather = read.table(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/weather_intermediate.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

daily_weather = read.table(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/daily_weather_afterstep0.txt',
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

#First flagged value assignment based on ratio humidity values present/nb obs per day + if weird number of daily obs (ex: 50, 100 etc) then the value is flagged too.
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

nomissing_humidity<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(sum_notNA=sum(!is.na(Relative.Humidity....)))


alldaymissing_humidity<-weather %>% 
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::summarise(all_missing_humidity=all(is.na(Relative.Humidity....)))
nomissing_humidity<-merge(nomissing_humidity,alldaymissing_humidity,by=c('Year_Exp','Day.of.Year'))

total_obs_daily<-merge(total_obs_daily,nomissing_humidity,by=c('Year_Exp','Day.of.Year'))



## Flagged Day.of.Year x Year_Exp values based on:

# 1) all records for that day are missing --> flagged

total_obs_daily<-total_obs_daily%>%group_by(Day.of.Year,Year_Exp)%>%mutate(flagged_humidity=case_when(all_missing_humidity==TRUE ~'flagged'))

# 2) 90% of the records for that day are not NA --> OK

total_obs_daily<-total_obs_daily%>%group_by(Day.of.Year,Year_Exp)%>%mutate(ratio=sum_notNA/n)%>%mutate(flagged_humidity=case_when(ratio<0.9 ~'flagged',
                                                                                                                              ratio>=0.9 ~'OK'))

# 3) Less than 24 observations per day, even if no NA present, should be removed.

total_obs_daily[total_obs_daily$n<24,'flagged_humidity']='flagged'

# 4) Except 2018_WIH1 and 2015_WIH2 for which many n==47 (manually checked and actually OK),
# we apply the following filters on n: when n is not a round value, something usually went wrong and data are not certain (unobserved values at unknown periods of the day for instance). Example For 2018_OHH1, from Day.of.Year 151 only round values for temp! Uncertainty therefore about the rest of the weather variables 
# + Some Year locations (identified in subset 2 above) very difficult to tell the reliability because of too many irregularities in the number of daily records. For these, we prefer to exclude the weather records.

total_obs_daily[total_obs_daily$n>96,'flagged_humidity']='flagged'
total_obs_daily[total_obs_daily$n>24&total_obs_daily$n<48,'flagged_humidity']='flagged'
total_obs_daily[total_obs_daily$n>48&total_obs_daily$n<72,'flagged_humidity']='flagged'
total_obs_daily[total_obs_daily$n>72&total_obs_daily$n<96,'flagged_humidity']='flagged'
#table(total_obs_daily[total_obs_daily$ratio>0.95&total_obs_daily$n>24&total_obs_daily$flagged_humidity=='flagged',1])
total_obs_daily[total_obs_daily$Year_Exp%in%c('2018_WIH1','2015_WIH2')&total_obs_daily$n==47,'flagged_humidity']='OK'
total_obs_daily[total_obs_daily$Year_Exp%in%c('2017_NEH4')&total_obs_daily$n%in%c(69,70,71)&total_obs_daily$ratio>0.95,'flagged_humidity']='OK'
total_obs_daily[total_obs_daily$Year_Exp%in%c('2018_DEH1')&total_obs_daily$ratio>0.95,'flagged_humidity']='OK'
total_obs_daily[total_obs_daily$Year_Exp%in%subset2$Year_Exp,'flagged_humidity']='flagged'

weather=merge(weather,total_obs_daily[,c('Day.of.Year','Year_Exp','flagged_humidity')],by=c('Day.of.Year','Year_Exp'),all.x=T)


write.table(total_obs_daily,'flagged_humidity.txt',col.names = T,row.names = F,sep = '\t',quote = F)
print(sum(total_obs_daily$flagged_humidity=='flagged')/nrow(total_obs_daily))

# ------------------------------------------------------------------------------
#Daily mean,max and min humidity computed for unflagged values:
# ------------------------------------------------------------------------------

maxH<-weather%>%
  filter(flagged_humidity%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(HMAX=max(Relative.Humidity....,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,HMAX,flagged_humidity)
maxH<-unique(maxH)
minH<-weather%>%
  filter(flagged_humidity%in%'OK')%>%
  group_by(Day.of.Year,Year_Exp)%>%
  dplyr::mutate(HMIN=min(Relative.Humidity....,na.rm = T))%>%
  dplyr::select(Day.of.Year,Year_Exp,HMIN)
minH<-unique(minH)
allhumidity=merge(maxH,minH,by=c('Day.of.Year','Year_Exp'),all.x=T)
allhumidity$HMEAN=(allhumidity$HMAX+allhumidity$HMIN)/2
allhumidity=arrange(allhumidity,Year_Exp,Day.of.Year)

# ------------------------------------------------------------------------------
#Add not flagged daily computed HMIN, HMAX, HMEAN to the daily_weather table (based on data from the field station)
# ------------------------------------------------------------------------------
daily_weather<-merge(daily_weather,allhumidity,by=c('Day.of.Year','Year_Exp'),all.x = T)
daily_weather=arrange(daily_weather,Year,Year_Exp,Day.of.Year)
print(nrow(daily_weather))
print(sum(daily_weather$flagged_humidity%notin%'OK')/nrow(daily_weather))

#Exclude for % Year_Exp without any weather station
locations_no_weather_station_used=c('2017_ONH1','2017_ONH2','2018_TXH1- Dry','2018_TXH1- Early','2018_TXH1- Late','2018_SCH1','2015_SDH1','2018_NEH2','2018_ILH1','2016_NEH1','2016_NEH4','2018_MIH1','2018_KSH1')
subset=daily_weather[daily_weather$Year_Exp%notin%locations_no_weather_station_used,]
print(nrow(subset))
print(sum(subset$flagged_humidity%notin%'OK')/nrow(subset))
daily_weather%>%
  group_by(Year_Exp)%>%
  summarize(allNA=all(is.na(HMIN)))


#Write the table which will be used for comparison with interpolated values for non missing values
write.table(daily_weather,'daily_weather_humidity_processed1.txt',col.names=T,row.names=F,sep='\t',quote=F)

