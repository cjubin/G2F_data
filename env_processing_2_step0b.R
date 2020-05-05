##############################################################################
################################ WEATHER DATA PROCESSING #######################
##############################################################################
rm(list=ls())


library(dplyr)
library(plyr)
library(lubridate)



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


##Add other variables not related to the weather (lon,lat, county, county, soil data, previous crop)
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

library(plyr)
rep.row <- function(r, n){
  colwise(function(x) rep(x, n))(r)
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
colnames(df)[c(1,2)]<-colnames(daily_weather)[c(1:2)]
daily_weather=df
rm(df,dd)

#Add months and day
daily_weather$month <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%m'))
daily_weather$day <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%d'))

setwd('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1')
write.table(daily_weather,file='daily_weather.txt',col.names = T,sep = '\t',quote = F,row.names = F)
