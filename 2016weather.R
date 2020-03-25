rm(list=ls())
library(elevatr)
library(rgdal)
library(tidyr)
library(intrval)
library(dplyr)
library(ggplot2)
library(nasapower)
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/2016")
weather=read.csv('g2f_2016_weather_data.csv',header = T,sep = ',')
weather$Field.Location=as.character(weather$Field.Location)
weather$day_year=paste0(weather$Month,'/',weather$Day,'/','16')
weather$day_year=as.POSIXct(as.Date(weather$day_year, "%m/%d/%Y"))
weather$day_year=difftime(weather$day_year,as.POSIXct(as.Date("01/01/16","%m/%d/%Y")),units='days')+1
weather$day_year=as.numeric(weather$day_year)

weather <- separate_rows(weather,Field.Location)

##Locations part of 'special lcoations' to be removed. See file.
##GAH2 not eliminated but remember that 4 harvest dates in a range of 5 days
weather=weather[weather$Field.Location%ni%c("IAH2","ILH2","NYH1",'ILH1','TXH2','ARH1','ARH2','GAH1','NYH2'),]


#Sum of solar radiation during each day and total rainfall/day
daily_weather=aggregate(weather[,12:13],by=list(weather$day_year, weather$Field.Location),FUN=sum,na.rm=T)
colnames(daily_weather)=c('day_year','location','sum.solar.radiation.W.m2','rainfall.mm')
daily_weather[daily_weather$sum.solar.radiation.W.m2==0,'sum.solar.radiation.W.m2']=NA

#################################################################
############# Add growing season start/end per environment#######
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/pheno_hybrids/2016")
pheno2016 <- read.csv("g2f_2016_hybrid_data_clean.csv", header=T, stringsAsFactors=F)
data=unique(pheno2016[,c("Field.Location","Date.Planted","Date.Harvested")])
data=data[data$Date.Harvested!='',]
data=data[data$Field.Location%in%unique(weather$Field.Location),]
row.names(data)=NULL

data=data[-c(14,21,22,13,24,26,17,18),]


daily_weather$Date.Planted=data[match(daily_weather$location,data$Field.Location),2]
daily_weather$Date.Harvested=data[match(daily_weather$location,data$Field.Location),3]


daily_weather$Date.Planted=as.POSIXct(as.Date(daily_weather$Date.Planted, "%m/%d/%Y"))
daily_weather$Date.Harvested=as.POSIXct(as.Date(daily_weather$Date.Harvested, "%m/%d/%Y"))

daily_weather$Date.Planted1=difftime(daily_weather$Date.Planted,as.POSIXct(as.Date("01/01/16","%m/%d/%Y")),units='days')+1
daily_weather$Date.Harvested1=difftime(daily_weather$Date.Harvested,as.POSIXct(as.Date("01/01/16","%m/%d/%Y")),units='days')+1


###############################
###############################


##Amount of measurements/day and amount of hourly missing data according to each weather feature

daily_weather$interval.measurement=aggregate(weather[,13],by=list(weather$day_year, weather$Field.Location),FUN=length)[,3]
daily_weather$sum.solar.radiation.MJ.m2=daily_weather$sum.solar.radiation.W.m2*24*3600/daily_weather$interval.measurement*1e-6
daily_weather$missing_data=aggregate(weather[,12:20],by=list(weather$day_year, weather$Field.Location),FUN=function(x)sum(is.na(x)))[,c(3:11)]



#Mean, max and min temperatures
daily_weather$max.temperature=aggregate(weather[,9],by=list(weather$day_year, weather$Field.Location),FUN=max,na.rm=T)[,3]
daily_weather$min.temperature=aggregate(weather[,9],by=list(weather$day_year, weather$Field.Location),FUN=min,na.rm=T)[,3]
daily_weather$mean.temperature=aggregate(weather[,9],by=list(weather$day_year, weather$Field.Location),FUN=mean,na.rm=T)[,3]


##Add linear GDD

daily_weather$GDD=(daily_weather$max.temperature+daily_weather$min.temperature)/2-10
daily_weather[which(daily_weather$GDD<0),'GDD']=0


#max and min relative humidity (%)
daily_weather$max.relative.humidity=aggregate(weather[,11],by=list(weather$day_year, weather$Field.Location),FUN=max,na.rm=T)[,3]
daily_weather$min.relative.humidity=aggregate(weather[,11],by=list(weather$day_year, weather$Field.Location),FUN=min,na.rm=T)[,3]

##Weird values WIH2 between 145 and 196 days

daily_weather[daily_weather$location=='WIH2'&daily_weather$day_year>144&daily_weather$day_year<196,]$max.relative.humidity=NA
daily_weather[daily_weather$location=='WIH2'&daily_weather$day_year>144&daily_weather$day_year<196,]$min.relative.humidity=NA

daily_weather[daily_weather%in%Inf]=NA


#vapour pressure deficit: (es - ea)
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/Rfunctions_weather/get_vpd.R")
daily_weather$vpd=get.vpd(rhmin = daily_weather$min.relative.humidity,rhmax = daily_weather$max.relative.humidity,tmin = daily_weather$min.temperature,tmax =daily_weather$max.temperature )


#daily dewpoint 
#daily_weather$dewpoint=temp - ((100 - daily_weather$max.relative.humidity)/5.)


##Wind speed average
daily_weather$wind.speed.ms.s=aggregate(weather[,14],by=list(weather$day_year, weather$Field.Location),FUN=mean,na.rm=T)[,3]


#Photothermal time & Daily mean photoperiod
#Daily mean photoperiod
daily_weather$photoperiod.hrs=unlist(aggregate(weather[,22],by=list(weather$day_year, weather$Field.Location),FUN=mean,na.rm=T)[,3])
daily_weather$photothermal.time=daily_weather$GDD*daily_weather$photoperiod.hrs



# Heat stress (threshold of 32 Celsius degrees)
# Binary feature : 1 if max. daily temperature superior to  32 celsius degrees and 0 if not
daily_weather$Sup32binary=NA
daily_weather[which(daily_weather$max.temperature>=32),'Sup32binary']=1
daily_weather[which(daily_weather$max.temperature<32),'Sup32binary']=0




# Crop Reference Evapotranspiration (ETref) - ETo from meteorological data - Unit: mm/day
# FAO Penman-Monteith method : determines the evapotrasnpiration from the hypothetical grass reference surface.
# This measure is usually considered equivalent to the computation of potential evapotranspiration (amount of evaporation and transpiration that would occur if a sufficient water source were available)

# g is the psychrometric constant
# elevation is needed (use of a package which gets these data from elevatr package - checked with real values for 2016)
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/2016")
info_fields=read.csv('g2f_2016_field_metadata.csv',sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
colnames(info_fields)[10]='x'
colnames(info_fields)[9]='y'
colnames(info_fields)[1]='Experiment'
info_fields=info_fields[info_fields$Experiment%in%unique(daily_weather$location),]
row.names(info_fields)=NULL
info_fields$x[15]=-82.728803
info_fields$y[15]=34.624744
elevation_data=get_elev_point(locations = info_fields[!is.na(info_fields$x),c(10,9)],units="meters", prj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
elevation=as.data.frame(elevation_data[1])
##Waterloo, Ontario
elevation[13,1]=328
elevation=as.data.frame(elevation)

info_fields=merge(info_fields,elevation,by=c('x','y'),all.x=T)
info_fields$Experiment=as.character(info_fields$Experiment)


daily_weather$elevation=info_fields[match(daily_weather$location,info_fields$Experiment),'elevation']
daily_weather$lon=info_fields[match(daily_weather$location,info_fields$Experiment),'x']
daily_weather$lat=info_fields[match(daily_weather$location,info_fields$Experiment),'y']

source('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/Rfunctions_weather/etp.R')
daily_weather$Rnet=NA
daily_weather$Rnet=net.sum.solar.radiation.MJ.m2(day_year = daily_weather$day_year,elevation = daily_weather$elevation,lat = daily_weather$lat,Rs_measured = daily_weather$sum.solar.radiation.MJ.m2,Tmax = daily_weather$max.temperature,Tmin = daily_weather$min.temperature,vpd = daily_weather$vpd,with_measurement = T)
daily_weather$ref.evapotranspiration=evapotranspiration(elevation = daily_weather$elevation,mean.T = daily_weather$mean.temperature,Rnet = daily_weather$Rnet,daily.mean.wind.speed.ms.s = daily_weather$wind.speed.ms.s,vpd = daily_weather$vpd)



daily_weather$diff.PminusETP=daily_weather$rainfall.mm-daily_weather$ref.evapotranspiration
daily_weather$diff.PminushalfETP=daily_weather$rainfall.mm-daily_weather$ref.evapotranspiration/2
daily_weather$dry=NA
daily_weather$dry[daily_weather$diff.PminushalfETP<0]='dry'
daily_weather$dry[daily_weather$diff.PminushalfETP>0&daily_weather$diff.PminusETP<0]='intermediate'
daily_weather$dry[daily_weather$diff.PminusETP>0]='wet'




##Treatment of missing values

#Add missing days start/during/end growing season
j=list()
loc=list()
lon=list()
lat=list()
elevation=list()
date.planted=list()
date.harvested=list()
h=1
for (i in unique(daily_weather$location)){
  dat=daily_weather[daily_weather$location==i,]
  vector_length=dat[1,'Date.Planted1']:dat[1,'Date.Harvested1']
  location=rep(i,length(vector_length))
  elevation[[h]]=rep(unique(daily_weather[daily_weather$location==i,'elevation']),length(vector_length))
  lon[[h]]=rep(unique(daily_weather[daily_weather$location==i,'lon']),length(vector_length))
  lat[[h]]=rep(unique(daily_weather[daily_weather$location==i,'lat']),length(vector_length))
  j[[h]]=vector_length
  loc[[h]]=location
  date.planted[[h]]=rep(dat[1,'Date.Planted1'],length(vector_length))
  date.harvested[[h]]=rep(dat[1,'Date.Harvested1'],length(vector_length))
  
  h=h+1
}
jj=unlist(j)
loc=unlist(loc)
date.planted=unlist(date.planted)
date.harvested=unlist(date.harvested)
lon=unlist(lon)
lat=unlist(lat)
elevation=unlist(elevation)

dat=cbind(loc,jj,date.planted,date.harvested,elevation,lon,lat)
colnames(dat)=c('location','day_year','Date.Planted1','Date.Harvested1','elevation','lon','lat')
dat=as.data.frame(dat)
daily_weather$day_year=as.factor(daily_weather$day_year)
daily_weather$location=as.factor(daily_weather$location)


daily_weather2=full_join(dat,daily_weather,by=c('location','day_year'))
daily_weather2$day_year=as.numeric(daily_weather2$day_year)





##
m1=as.numeric(which(daily_weather2$interval.measurement<24))
m2=as.numeric(which(daily_weather2$interval.measurement>24&daily_weather2$interval.measurement<48))
m3=as.numeric(which(daily_weather2$interval.measurement>48&daily_weather2$interval.measurement<96))
m4=as.numeric(which(daily_weather2$interval.measurement>96))
m=sort(c(m1,m2,m3,m4))
daily_weather2[m,c(8:ncol(daily_weather2))]=NA

missing.rain=as.numeric(which(daily_weather2$missing_data$Rainfall..mm.>3))
missing.rain2=as.numeric(which(daily_weather2$rainfall.mm>90))
daily_weather2$rainfall.mm[c(missing.rain,missing.rain2)]=NA

missing.temp=as.numeric(which(daily_weather2$missing_data$Temperature..C.>5))
daily_weather2$max.temperature[missing.temp]=NA
daily_weather2$min.temperature[missing.temp]=NA
daily_weather2$mean.temperature[missing.temp]=NA
daily_weather2$GDD[missing.temp]=NA
daily_weather2$Sup32binary[missing.temp]=NA


missing.sum.solar.rad=as.numeric(which(daily_weather2$missing_data$Solar.Radiation..W.m2.>3))
missing.sum.solar.rad2=as.numeric(which(daily_weather2$sum.solar.radiation.MJ.m2>38))
missing.sum.solar.rad3=as.numeric(which(daily_weather2$sum.solar.radiation.MJ.m2<8))
daily_weather2$sum.solar.radiation.W.m2[missing.sum.solar.rad]=NA
daily_weather2$sum.solar.radiation.MJ.m2[missing.sum.solar.rad]=NA
daily_weather2$sum.solar.radiation.W.m2[missing.sum.solar.rad2]=NA
daily_weather2$sum.solar.radiation.MJ.m2[missing.sum.solar.rad2]=NA

missing.wind.speed=as.numeric(which(daily_weather2$missing_data$Wind.Speed..m.s.>5))
daily_weather2$wind.speed.ms.s[missing.wind.speed]=NA

missing.relative.humidity=as.numeric(which(daily_weather2$missing_data$Relative.Humidity....>5))
daily_weather2$max.relative.humidity[missing.relative.humidity]=NA
daily_weather2$min.relative.humidity[missing.relative.humidity]=NA
daily_weather2$vpd[missing.relative.humidity]=NA
daily_weather2$ref.evapotranspiration[missing.relative.humidity]=NA



daily_weather2$photoperiod.hrs=NA
daily_weather2$photothermal.time=NA

daily_weather2=daily_weather2 %>% mutate_if(is.numeric, list(~na_if(., Inf)))
daily_weather2=daily_weather2 %>% mutate_if(is.numeric, list(~na_if(.,-Inf)))
daily_weather2=daily_weather2 %>% mutate_if(is.numeric, list(~na_if(., NaN)))

############################################################
############################################################
##OBTAIN DATA FROM NASA METEOROLOGICAL DATASETS

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/2016")


agdtotal=read.table('agd2016.txt',header = T,fill = T,check.names = F)
colnames(agdtotal)[c(1,7)]=c('location','day_year')

########################################################
##Adding values from NASA to replace NA weather values

d1=daily_weather2[which(is.na(daily_weather2$rainfall.mm)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,13)],by=c('location','day_year'))
daily_weather2[which(is.na(daily_weather2$rainfall.mm)),'rainfall.mm']=d2$PRECTOT

d1=daily_weather2[which(is.na(daily_weather2$min.temperature)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,12)],by=c('location','day_year'))
daily_weather2[which(is.na(daily_weather2$min.temperature)),'min.temperature']=d2$T2M_MIN

d1=daily_weather2[which(is.na(daily_weather2$max.temperature)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,11)],by=c('location','day_year'))
daily_weather2[which(is.na(daily_weather2$max.temperature)),'max.temperature']=d2$T2M_MAX

d1=daily_weather2[which(is.na(daily_weather2$mean.temperature)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,10)],by=c('location','day_year'))
daily_weather2[which(is.na(daily_weather2$mean.temperature)),'mean.temperature']=d2$T2M

d1=daily_weather2[which(is.na(daily_weather2$wind.speed.ms.s)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,14)],by=c('location','day_year'))
daily_weather2[which(is.na(daily_weather2$wind.speed.ms.s)),'wind.speed.ms.s']=d2$WS2M


d1=daily_weather2[which(is.na(daily_weather2$sum.solar.radiation.MJ.m2)),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,15)],by=c('location','day_year'),all.x=T)
daily_weather2[which(is.na(daily_weather2$sum.solar.radiation.MJ.m2)),'sum.solar.radiation.MJ.m2']=d2$ALLSKY_SFC_SW_DWN

d1=daily_weather2[which(daily_weather2$sum.solar.radiation.MJ.m2<0),c('location','day_year')]
d2=merge(d1,agdtotal[,c(1,7,16)],by=c('location','day_year'),all.x=T)
daily_weather2[which(daily_weather2$sum.solar.radiation.MJ.m2<0),'sum.solar.radiation.MJ.m2']=d2$ALLSKY_SFC_LW_DWN


d1=daily_weather2[which(is.na(daily_weather2$vpd)),c('location','day_year','min.temperature','max.temperature')]
d2=merge(d1,agdtotal[,c(1,7,9)],by=c('location','day_year'))
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/Rfunctions_weather/get_vpd_withmeanRH.R")
daily_weather2[which(is.na(daily_weather2$vpd)),'vpd']=get.vpd.with.meanRH(rhmean=d2$RH2M,tmin = d2$min.temperature,tmax =d2$max.temperature)

##################
d1=daily_weather2[which(is.na(daily_weather2$Date.Planted1.x)),c('location','day_year')]
row.names(d1)=NULL
d1$Date.Planted1.x=NA
d1$Date.Harvested1.x=NA
d1$elevation.x=NA
d1$lon.x=NA
d1$lat.x=NA
for (i in 1:nrow(d1)){
  m=unique(daily_weather2[daily_weather2$location==d1$location[i],'Date.Planted1.x'])
  s=unique(daily_weather2[daily_weather2$location==d1$location[i],'Date.Harvested1.x'])
  e=unique(daily_weather2[daily_weather2$location==d1$location[i],'elevation.x'])
  lon=unique(daily_weather2[daily_weather2$location==d1$location[i],'lon.x'])
  lat=unique(daily_weather2[daily_weather2$location==d1$location[i],'lat.x'])
  d1$elevation.x[i]<-as.vector(e[!is.na(e)])
  d1$lon.x[i]<-as.vector(lon[!is.na(lon)])
  d1$lat.x[i]<-as.vector(lat[!is.na(lat)])
  d1$Date.Planted1.x[i]<-as.vector(m[!is.na(m)])
  d1$Date.Harvested1.x[i]<-as.vector(s[!is.na(s)])
}
v=which(is.na(daily_weather2$Date.Planted1.x))
daily_weather2[v,'Date.Planted1.x']=d1$Date.Planted1.x
daily_weather2[v,'Date.Harvested1.x']=d1$Date.Harvested1.x
daily_weather2[v,'lon.x']=d1$lon.x
daily_weather2[v,'lat.x']=d1$lat.x
daily_weather2[v,'elevation.x']=d1$elevation.x

daily_weather2$elevation.x=as.numeric(as.vector(daily_weather2$elevation.x))
daily_weather2$lat.x=as.numeric(as.vector(daily_weather2$lat.x))
daily_weather2$lon.x=as.numeric(as.vector(daily_weather2$lon.x))



################################################
################################################
####Recomputing some specific values (ETP,GDD,photoperiod and photothermal time)


#photoperiod and photothermal time
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/Rfunctions_weather/daylength.R")
v=which(is.na(daily_weather2$photoperiod.hrs))
daily_weather2$day_year=as.numeric(daily_weather2$day_year)
daily_weather2$lat.x=as.numeric(daily_weather2$lat.x)
daily_weather2$photoperiod.hrs[v]=daylength(lat = daily_weather2$lat.x[v],day_of_year = daily_weather2$day_year[v])
daily_weather2$GDD=(daily_weather2$max.temperature+daily_weather2$min.temperature)/2-10
daily_weather2[which(daily_weather2$GDD<0),'GDD']=0
daily_weather2$photothermal.time=daily_weather2$GDD*daily_weather2$photoperiod.hrs


#ref. evapotranspiration
daily_weather2$elevation.x=as.numeric(daily_weather2$elevation.x)
daily_weather2$Rnet=net.sum.solar.radiation.MJ.m2(day_year = daily_weather2$day_year,elevation = daily_weather2$elevation.x,lat = daily_weather2$lat.x,Rs_measured = daily_weather2$sum.solar.radiation.MJ.m2,Tmax = daily_weather2$max.temperature,Tmin = daily_weather2$min.temperature,vpd = daily_weather2$vpd,with_measurement = T)
daily_weather2$ref.evapotranspiration=evapotranspiration(elevation = daily_weather2$elevation.x,mean.T = daily_weather2$mean.temperature,Rnet = daily_weather2$Rnet,daily.mean.wind.speed.ms.s = daily_weather2$wind.speed.ms.s,vpd = daily_weather2$vpd)
daily_weather2$diff.PminusETP=daily_weather2$rainfall.mm-daily_weather2$ref.evapotranspiration
daily_weather2$diff.PminushalfETP=daily_weather2$rainfall.mm-daily_weather2$ref.evapotranspiration/2
daily_weather2$dry=NA
daily_weather2$dry[daily_weather2$diff.PminushalfETP<0]='dry'
daily_weather2$dry[daily_weather2$diff.PminushalfETP>0&daily_weather2$diff.PminusETP<0]='intermediate'
daily_weather2$dry[daily_weather2$diff.PminusETP>0]='wet'







# Heat stress (threshold of 32 Celsius degrees)
# Binary feature : 1 if max. daily temperature superior to  32 celsius degrees and 0 if not
daily_weather2$Sup32binary=NA
daily_weather2[which(daily_weather2$max.temperature>=32),'Sup32binary']=1
daily_weather2[which(daily_weather2$max.temperature<32),'Sup32binary']=0


## Subsetting by variables which will be used further as covariates

daily_weather_final=daily_weather2[,c(1:7,9,15,17:27,31:35)]
daily_weather_final=cbind(daily_weather_final,year='2016')



library(dplyr)
df=daily_weather_final %>%
  group_by(location)  %>%
  mutate(cum_rainfall_2016 = cumsum(rainfall.mm))



df=df %>%
  group_by(location) %>%
  mutate(cum_GDD_2016 = cumsum(GDD))

##Subset to the beginning to the end of the growing season in each location.
##Remove values after the harvest
vv=vector()
j=1

for (i in 1:nrow(df)){
  if (df$day_year[i]<(as.numeric(as.vector(df$Date.Planted1.x[i])))|df$day_year[i]>(as.numeric(as.vector(df$Date.Harvested1.x[i])))){
    vv[j]=i
    j=j+1}
}
df=df[-vv,]

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/2016")
write.table(df,'daily_weather2016.txt',col.names=T,row.names=F,sep='\t',quote=F)



v=ggplot(df, aes(x=day_year, y=cum_rainfall_2016,col=location)) + geom_line() + geom_point()
s=ggplot(df, aes(x=day_year, y=cum_GDD_2016,col=location)) + geom_line() + geom_point()
vv=list(v,s)
pdf('cumulative_rain.pdf',onefile=T,height = 12,width=15)
print(v)
dev.off()
pdf('cumulative_GDD.pdf',onefile=T,height = 12,width=15)
print(s)
dev.off()


ggplot(df, aes(x=day_year, y=ref.evapotranspiration,col=location)) + geom_line() + geom_point()
ggsave('et0_reference.pdf',width = 20,height = 25)

df[,c(1,2,8:16)] %>%
  gather(-day_year, -location, key = "var", value = "value") %>% 
  ggplot(aes(x = day_year, y = value, color = location)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
ggsave('allweathervariablesplot1.pdf',width = 20,height = 15)



df[,c(1,2,17:24)] %>%
  gather(-day_year, -location, key = "var", value = "value") %>% 
  ggplot(aes(x = day_year, y = value, color = location)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
ggsave('allweathervariablesplot2.pdf',width = 20,height = 15)

df[,c(1,2,10:13,16,17,18,21,22)] %>%
  gather(-day_year, -location, key = "var", value = "value") %>% 
  ggplot(aes(x = day_year, y = value, color = location)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
ggsave('allweathervariablesplot3.pdf',width = 20,height = 15)




sapply(df, function(x)length(which(is.na(x))))
