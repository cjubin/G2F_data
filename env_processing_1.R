### Load useful packages

rm(list=ls())
library(elevatr)
library(rgdal)
library(tidyr)
library(readxl)
library(intrval)
library(dplyr)
library(ggplot2)
library(rnoaa)
library(stringr)
library(ggmap)
library(revgeo)
library(pracma)
library(FedData)
library(tidyverse)
source('latlong2county.R')
`%notin%` <- Negate(`%in%`)
options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")

##Import 5 years of weather data
setwd("~/Final_datasets_G2F/ALL_WEATHER")
weather2014=read.csv('g2f_2014_weather.csv',header = T,sep = ',')
weather2015=read.csv('g2f_2015_weather.csv',header = T,sep = ',')
weather2016=read.csv('g2f_2016_weather_data.csv',header = T,sep = ',')
weather2017=read.csv('g2f_2017_weather_data.csv',header = T,sep = ',')
weather2018=read.csv('g2f_2018_weather_clean.csv',header = T,sep = ',')[1:28]

####Setting the exact same columns names across the 5 datasets: changing names, removal of some columns + Add Year column if not present
colnames(weather2014)[which(colnames(weather2014)%notin%colnames(weather2015))]
colnames(weather2014)[c(6,7,8,9)]=c('Day','Month','Year','Day.of.Year')
colnames(weather2015)[c(3)]='Field.Location'
colnames(weather2015)[21]=c('Soil.Moisture...VWC.')
colnames(weather2016)[c(13)]=c('Rainfall..mm.')
colnames(weather2016)[22]=c("Photoperiod..hours.")
colnames(weather2018)[24]=c("Photoperiod..hours.")
weather2015=weather2015[,-11]
colnames(weather2015)[10]='Time.Local.'
weather2014=weather2014[,-11]
colnames(weather2014)[10]='Time.Local.'
colnames(weather2017)[8]='Time.Local.'
colnames(weather2018)[9]='Time.Local.'
colnames(weather2018)[21]='UV.Light..uM.m2s.'

##Add day of year column for the annual datasets which do not have it
weather2016$Day.of.Year=paste0(weather2016$Month,'/',weather2016$Day,'/','16')
weather2016$Day.of.Year=as.POSIXct(as.Date(weather2016$Day.of.Year, "%m/%d/%Y"))
weather2016$Day.of.Year=difftime(weather2016$Day.of.Year,as.POSIXct(as.Date("01/01/16","%m/%d/%Y")),units='days')+1
weather2016$Day.of.Year=as.numeric(weather2016$Day.of.Year)
weather2017$Day.of.Year=paste0(weather2017$Month,'/',weather2017$Day,'/','17')
weather2017$Day.of.Year=as.POSIXct(as.Date(weather2017$Day.of.Year, "%m/%d/%Y"))
weather2017$Day.of.Year=difftime(weather2017$Day.of.Year,as.POSIXct(as.Date("01/01/17","%m/%d/%Y")),units='days')+1
weather2017$Day.of.Year=as.numeric(weather2017$Day.of.Year)
weather2018$Day.of.Year=paste0(weather2018$Month,'/',weather2018$Day,'/','18')
weather2018$Day.of.Year=as.POSIXct(as.Date(weather2018$Day.of.Year, "%m/%d/%Y"))
weather2018$Day.of.Year=difftime(weather2018$Day.of.Year,as.POSIXct(as.Date("01/01/18","%m/%d/%Y")),units='days')+1
weather2018$Day.of.Year=as.numeric(weather2018$Day.of.Year)

##Merging the annual weather files based on a subset of columns
names_col=c("Field.Location" ,"Station.ID","NWS.Network","Month","Day","Year",'Day.of.Year','Time.Local.','Temperature..C.','Dew.Point..C.','Relative.Humidity....',"Solar.Radiation..W.m2.","Rainfall..mm.",'Wind.Speed..m.s.','Wind.Gust..m.s.',  "Photoperiod..hours.","Column.Altered" ,"Altered.Column.Names", "Cleaning.Method" ,'Comment')              
weather2014bis=weather2014[,names_col]
weather2015bis=weather2015[,names_col]
weather2016bis=weather2016[,names_col]
weather2017bis=weather2017[,names_col]
weather2018bis=weather2018[,names_col]
weather_all=rbind(weather2014bis,weather2015bis,weather2016bis,weather2017bis,weather2018bis)
weather_all <- separate_rows(weather_all,'Field.Location')
#write.table(weather_all,'weather_all.txt',sep = '\t',quote = F,col.names = T,row.names = F)
rm(weather2014,weather2014bis,weather2015,weather2015bis,weather2016,weather2016bis,weather2017,weather2017bis,weather2018,weather2018bis)

##################################################################################
####Loading the field metadata files: geographical coordinates and name cities etc
setwd("~/Final_datasets_G2F/ALL_WEATHER")
field2014=read.csv('g2f_2014_field_characteristics.csv',header = T,sep = ',')
field2014=field2014[-which(field2014$Type=='inbred'),]
field2015=read.csv('g2f_2015_field_metadata.csv',header = T,sep = ',')
field2016=read.csv('g2f_2016_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2017=read.csv('g2f_2017_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2018=read.csv('g2f_2018_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2014$Year=2014
field2015$Year=2015
field2016$Year=2016
field2017$Year=2017
field2018$Year=2018

#Retrieving Lon/Lat data for each field experiment --> use of the field corner instead of the WS coordinates
#For retrieving soil data: more accurate to use field than weather location, which can be very different
#If the Weather station had dubious data, then we will use a nearby location to the field instead of one near to the weather station.
colnames(field2015)[c(19,20)]=c('lat','long')
colnames(field2016)[c(23,24)]=c('lat','long')
colnames(field2016)[1]='Experiment'
colnames(field2017)[c(26,27)]=c('lat','long')
field2017[29,c('lat','long')]=c(30.54535,	-96.43258)
colnames(field2017)[1]='Experiment'
colnames(field2018)[c(26,27)]=c('lat','long')
field2018[6,c('lat','long')]=c(41.19870,	-91.48618)
field2018[7,c('lat','long')]=c(42.06593,	-94.72742)
field2018[8,c('lat','long')]=c(41.98745,	-92.26014)
field2018[9,c('lat','long')]=c(41.99775,	-93.69625)

colnames(field2018)[1]='Experiment'


geofield2014=field2014[,c('lat','long','Year','Experiment','City')]
geofield2015=field2015[,c('lat','long','Year','Experiment','City')]
geofield2016=field2016[,c('lat','long','Year','Experiment','City')]
geofield2017=field2017[,c('lat','long','Year','Experiment','City')]
geofield2018=field2018[,c('lat','long','Year','Experiment','City')]
geofield2014$Experiment=as.character(as.vector(geofield2014$Experiment))
geofield2014[22,'Experiment']=as.character('WIH1')
geofield2018[c(15,16),'Experiment']='MOH1'
geo_data_fields=rbind(geofield2014,geofield2015,geofield2016,geofield2017,geofield2018)
colnames(geo_data_fields)[4]='Field.Location'


#Missing coordinates 
#For those with at least the city,nand a unique field location in the other years --> imputation possible
geo_data_fields$lonlat_added='NO'
geo_data_fields[82,c('lat','long','lonlat_added')]=c(41.19869,	-91.48620,'YES')
geo_data_fields[c(30),c('lat','long','lonlat_added')]=c(42.06591,	-94.72745,'YES')
geo_data_fields[c(83),c('lat','long','lonlat_added')]=c(42.06591,	-94.72745,'YES')
geo_data_fields[c(84),c('lat','long','lonlat_added')]=c(41.98738,	-92.26016,'YES')
geo_data_fields[c(85),c('lat','long','lonlat_added')]=c(41.99750,	-93.69990,'YES')
geo_data_fields[c(93),c('lat','long','lonlat_added')]=c(41.16636,	-96.41726,'YES')#NEH1 2016: assume to be the same field as the year before
##TXH2: except 2014 no info on geographical coordinates AT ALL
geo_data_fields[c(101),c('lat','long','lonlat_added')]=c(34.62261,	-82.73796,'YES')#SCH1 2016: suppose same location as 2017
geo_data_fields[c(116),c('lat','long','lonlat_added')]=c(40.06119,	-88.23327,'YES')#ILH1 2017: assume same location as previous years too
geo_data_fields[c(118),c('lat','long','lonlat_added')]=c(40.47835,	-86.99013,'YES')#INH1 2017: assume same location as previous years too
geo_data_fields[c(167),c('lat','long','lonlat_added')]=c(43.30480,	-89.38520,'YES')#change WIH2 2018
geo_data_fields[c(94),c('lat','long','lonlat_added')]=c(41.167022,	-96.417192,'YES')#change NEH4 based on 2015, same (900 m distance) in 2017 too
geo_data_fields=geo_data_fields[!is.na(geo_data_fields$lat),]
geo_data_fields$lat=as.numeric(geo_data_fields$lat)
geo_data_fields$long=as.numeric(geo_data_fields$long)



#Use revgeo to add county, ZIP code, and city
results_couties<-latlong2county(data.frame(geo_data_fields$long[!is.na(geo_data_fields$long)],geo_data_fields$lat[!is.na(geo_data_fields$lat)]))
geo_data_fields$state=sapply(strsplit(results_couties,","), `[`, 1)
geo_data_fields$county=sapply(strsplit(results_couties,","), `[`, 2)
results_cities=revgeo(longitude=geo_data_fields$long[!is.na(geo_data_fields$long)], latitude=geo_data_fields$lat[!is.na(geo_data_fields$lat)], provider = 'photon',item='city', output='frame')

geo_data_fields$city_revgeo=results_cities$city
geo_data_fields$city_revgeo[is.na(geo_data_fields$city_revgeo)]=geo_data_fields$City[is.na(geo_data_fields$city_revgeo)]
geo_data_fields$zip=results_cities$zip



#Merge weather data with geographical coordinates information
weather=merge(geo_data_fields,weather_all,by=c('Year','Field.Location'), all.y = T)##Merging the weather data with Lon/lat + name cities 
rm(weather_all)



###################ELIMINATE SOME EXPERIMENTS###################################
#Eliminate those corresponding to inbred lines experiments, disease trials
##NYH1 disease trial with no info, TXH2 missing --> need to be removed, absolutely no info on this trial
##The rest are inbred field trials
missing_lonlat=unique(weather[which(is.na(weather$long)),c('Year','Field.Location')])
weather=weather[-which(is.na(weather$long)),]

to_remove_inbred_trials=c("AZI1", "AZI2", "DEI1", "GAI1", "IAI1", "IAI2", "IAI3", "IAI4","ILI1", "INI1","KSI1" ,"MNI1" ,"MOI1", "MOI2", "NCI1" ,"NYI2", "PAI1" ,"PAI2" ,"TXI1" ,"TXI2" , "WII1", "WII2" )
weather=weather[-which(weather$Field.Location%in%to_remove_inbred_trials),]

weather$Year_Exp=paste(weather$Year,weather$Field.Location,sep = '_')





##############################################################
################## SOIL INFORMATION ##########################
#######Info very limited for years 2014 + 2015 --> need to use external package : feddata
###Load the soil data files to add the soil texture (at least)
soil2014=as.data.frame(cbind('Year'=2014,'Experiment'=as.character(field2014$Experiment),'Soil.test.type'=as.character(as.vector(field2014$Soil.test.type)),'Soil.texture'=as.character(as.vector(field2014$Soil.texture)),'Soil.pH'=as.character(as.vector(field2014$Soil.pH))))
soil2014$Year_Exp=paste(soil2014$Year,soil2014$Experiment,sep = '_')
soil2015=read.csv('g2f_2015_soil_data.csv',header = T,sep = ',')
soil2015$X..Sand=NA
soil2015$X..Silt=NA
soil2015$X..Clay=NA
soil2015$Texture=NA
soil2015$Year=2015
soil2015$Year_Exp=paste(soil2015$Year,soil2015$Experiment,sep = '_')
soil2016=read.csv('g2f_2016_soil_data_clean.csv',header = T,sep = ',',stringsAsFactors = F)
soil2016$Year=2016
colnames(soil2016)[10]='OM'
soil2016$Year_Exp=paste(soil2016$Year,soil2016$Location,sep = '_')
soil2017=read.csv('g2f_2017_soil_data_clean.csv',header = T,sep = ',',stringsAsFactors = F)
soil2017$Year=2017
colnames(soil2017)[9]='OM'
soil2017$Year_Exp=paste(soil2017$Year,soil2017$Location,sep = '_')
soil2018=as.data.frame(readxl::read_xlsx('g2f_2018_soil_data.xlsx'))
soil2018$Year=2018
soil2018$Year_Exp=paste(soil2018$Year,soil2018$`Field ID`,sep = '_')
colnames(soil2018)[c(25,26,27)]=c('X..Sand','X..Silt','X..Clay')
colnames(soil2018)[10]='OM'
#Soil texture (no data for 2015, no composition percentage for 2014)
colnames(soil2014)[4]='Texture'
soil2014$X..Sand=NA
soil2014$X..Silt=NA
soil2014$X..Clay=NA
soil2014$OM=NA


soil_1=rbind(soil2014[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2015[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2016[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2017[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2018[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')])
soildata=as.data.frame(unique(cbind('Year_Exp'=weather$Year_Exp,'Year'=weather$Year,'City'=as.character(weather$city_revgeo),'lat'=weather$lat,'long'=weather$long,'lonlat_imputed'=weather$lonlat_added,'Field.Location'=as.character(weather$Field.Location))))
soildata=merge(soildata,soil_1,by='Year_Exp',all.x=T)

#Removing duplicated lines for same Year_Loc
soildata=soildata[-which(duplicated(soildata[,1])),]
#Create column to indicate that the data were imputed form other years of experiments
soildata$imputed=NA
soildata[!is.na(soildata$X..Sand),'imputed']='NO'

soildata[soildata$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay')]=c(11,30,59)

soildata[soildata$Year_Exp=='2015_ONH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_ONH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp=='2015_ONH1','imputed']='close_location_other_years'
soildata[soildata$Year_Exp=='2018_DEH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_DEH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp=='2018_DEH1','imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_GAH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2016_IAH1','2018_IAH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2017_IAH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2016_IAH1','2018_IAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_IAH3'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2017_IAH3',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_IAH3'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_MOH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2018_MOH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2014_MOH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_MOH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_MOH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_MOH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_ILH1','2015_ILH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_ILH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2014_ILH1','2015_ILH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_MNH1','2015_MNH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_MNH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2014_MNH1','2015_MNH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2017_NCH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_NCH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2016_NCH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2014_NCH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_OHH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2017_OHH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_OHH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_NYH2','2014_NYH2'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=soildata[soildata$Year_Exp=='2017_NYH2',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
soildata[soildata$Year_Exp%in%c('2015_NYH2','2014_NYH2'),'imputed']='close_location_other_years'



#Looking for missing soil information
setwd("~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1")
write.table(soildata,'soildata.txt',sep = '\t',quote = F,row.names = F,col.names = T)

soil_ssurgo=readxl::read_xlsx('soil_imputed_wss.xlsx')

#



###Add previous crop on the field (agronomic feature) or rotation info




##Check name cities Metadata file with real cities names from revgeo

##Remove inbred-related field trials weather information

################################
################################
##CREATE A DAILY WEATHER DATASET

