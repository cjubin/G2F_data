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
library(revgeo)
library(FedData)
`%notin%` <- Negate(`%in%`)

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
####Attribute Lon/Lat data using field metadata files + name cities + name counties for each field experiment 
setwd("~/Final_datasets_G2F/ALL_WEATHER")
field2014=read.csv('g2f_2014_field_characteristics.csv',header = T,sep = ',')
field2015=read.csv('g2f_2015_field_metadata.csv',header = T,sep = ',')
field2016=read.csv('g2f_2016_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2017=read.csv('g2f_2017_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2018=read.csv('g2f_2018_field_metadata.csv',header = T,sep = ',',stringsAsFactors=FALSE, fileEncoding="latin1")
field2014$Year=2014
field2015$Year=2015
field2016$Year=2016
field2017$Year=2017
field2018$Year=2018

#Retrieving Lon/Lat data for each field experiment
names_col2=c('lat','long')
colnames(field2015)[c(6,7)]=c('lat','long')
colnames(field2016)[c(9,10)]=c('lat','long')
colnames(field2016)[1]='Experiment'
colnames(field2017)[c(9,10)]=c('lat','long')
colnames(field2017)[1]='Experiment'
colnames(field2018)[c(9,10)]=c('lat','long')
colnames(field2018)[1]='Experiment'
field2018[field2018$Experiment=='WIH2',c('lat','long')]=c('43.324231'	,'-89.33564')
field2016[field2016$Experiment=='INH1',c('lat','long')]=c('40.47925'	,'-86.99013')
geofield2014=field2014[,c('lat','long','Year','Experiment','City')]
geofield2015=field2015[,c('lat','long','Year','Experiment','City')]##NYH1 disease trial with no info, TXH2 missing --> need to be removed, absolutely no info on this trial
geofield2016=field2016[,c('lat','long','Year','Experiment','City')]
geofield2017=field2017[,c('lat','long','Year','Experiment','City')]
geofield2018=field2018[,c('lat','long','Year','Experiment','City')]
geofield2018[c(15,16),'Experiment']='MOH1'
geo_data_fields=rbind(geofield2014,geofield2015,geofield2016,geofield2017,geofield2018)
colnames(geo_data_fields)[4]='Field.Location'

##Merge weather data with geographical coordinates
weather=merge(geo_data_fields,weather_all,by=c('Year','Field.Location'), all.y = T)##Merging the weather data with Lon/lat + name cities 
rm(weather_all)

##If missing global Lon/Lat for some locations, eliminate trial from analyses depending on the case (Cf 'specific locations' document)
## Or if just problem with the columns: add geographical coordinates values from corner lower left
#Table missing 
missing_lonlat=unique(weather[which(is.na(weather$long)),c('Year','Field.Location')])
missing_lonlat$Year_Exp=paste(missing_lonlat$Year,missing_lonlat$Field.Location,sep = '_')
row.names(missing_lonlat)=NULL
#Create column Year_Field.Location and eliminate those corresponding to inbred lines experiments, disease trials, or trials which should be removed
#MNH1_2017: bad phenotypic quality + weird location regarding geographical position in 2017 (MADISON ?)
#2014_NEH3 ,2017_NEH3 and 2017_NEH4 have no silking date
#2015_NEH1 has no info on irrigation, 2015_NEH2 has no plant height, 2015_NEH3 has no plant height and  2015_NEH4 has no silking date
#SCH1: no grain yield data
weather$Year_Exp=paste(weather$Year,weather$Field.Location,sep = '_')
to_remove=missing_lonlat$Year_Exp[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30)]
weather=weather[-which(weather$Year_Exp%in%to_remove),]
weather=weather[-which(weather$Year_Exp%in%c('2017_MNH1','2014_NEH3','2015_NEH1','2015_NEH2','2015_NEH3','2015_NEH4','2017_NEH3',"2015_GAI1","2015_INI1","2015_TXI2","2015_IAI1","2015_NYI2","2015_KSI1","2015_IAI2","2015_IAI3","2015_MNI1","2015_WII1","2015_AZI1","2015_IAI4","2015_MOI1","2015_PAI1","2015_WII2", "2015_PAI2", "2015_MOI2","2015_AZI2","2015_DEI1","2015_ILI1","2015_NCI1","2015_TXI1")),]
weather=weather[-which(weather$Year_Exp%in%c('2014_NYH1','2015_NYH1','2016_NYH1','2017_NYH1','2018_NYH1')),]


##############################################################
################## SOIL INFORMATION ##########################
#######Info very limited for years 2014 + 2015 --> need to use external package : feddata
###Load the soil data files to add the soil texture (at least)
soil2014=as.data.frame(cbind('Year'=2014,'Experiment'=as.character(field2014$Experiment),'Soil.test.type'=as.character(as.vector(field2014$Soil.test.type)),'Soil.texture'=as.character(as.vector(field2014$Soil.texture)),'Soil.pH'=as.character(as.vector(field2014$Soil.pH))))
soil2014$Year_Exp=paste(soil2014$Year,soil2014$Experiment,sep = '_')
soil2015=read.csv('g2f_2015_soil_data.csv',header = T,sep = ',')
soil2016=read.csv('g2f_2016_soil_data_clean.csv',header = T,sep = ',',stringsAsFactors = F)
soil2016$Year=2016
soil2016$Year_Exp=paste(soil2016$Year,soil2016$Location,sep = '_')
soil2017=read.csv('g2f_2017_soil_data_clean.csv',header = T,sep = ',',stringsAsFactors = F)
soil2017$Year=2017
soil2017$Year_Exp=paste(soil2017$Year,soil2017$Location,sep = '_')
soil2018=as.data.frame(readxl::read_xlsx('g2f_2018_soil_data.xlsx'))
soil2018$Year=2018
soil2018$Year_Exp=paste(soil2018$Year,soil2018$`Field ID`,sep = '_')
#Soil texture (no data for 2015, no composition percentage for 2014)
colnames(soil2014)[4]='Texture'
soil2014$X..Sand=NA
soil2014$X..Silt=NA
soil2014$X..Clay=NA

colnames(soil2018)[c(25,26,27)]=c('X..Sand','X..Silt','X..Clay')
soil_1=rbind(soil2014[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture')],soil2016[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture')],soil2017[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture')],as.data.frame(soil2018[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture')]))
soildata=as.data.frame(unique(cbind('Year_Exp'=weather$Year_Exp,'Year'=weather$Year,'City'=as.character(weather$City),'lat'=weather$lat,'long'=weather$long,'Field.Location'=as.character(weather$Field.Location))))
soildata=merge(soildata,soil_1,by='Year_Exp',all.x=T)

#Removing duplicated lines
soildata=soildata[-which(duplicated(soildata[,1])),]
soildata$imputed=NA
soildata[!is.na(soildata$X..Sand),'imputed']='NO'

soildata[soildata$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay')]=c(11,30,59)
soildata[soildata$Year_Exp=='2014_TXH1','imputed']='NO'
soildata[soildata$Year_Exp=='2015_ONH1',c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_ONH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp=='2015_ONH1','imputed']='close_location_other_years'
soildata[soildata$Year_Exp=='2018_DEH1',c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_DEH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp=='2018_DEH1','imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_GAH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_IAH2','2016_IAH2','2018_IAH2'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_IAH2',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_IAH2','2016_IAH2','2018_IAH2'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2016_IAH1','2018_IAH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_IAH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2016_IAH1','2018_IAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_IAH3'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_IAH3',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2016_IAH1','2018_IAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_IAH3'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_IAH3',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_IAH3'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_MOH1','2014_MOH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_MOH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_MOH1','2014_MOH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_ILH1','2015_ILH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_ILH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_ILH1','2015_ILH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_INH1','2015_INH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_INH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_INH1','2015_INH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_MNH1','2015_MNH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_MNH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_MNH1','2015_MNH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_NCH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_NCH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2014_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2016_NCH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2014_NCH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2016_IAH4','2018_IAH4','2014_IAH1','2015_IAH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_IAH4',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2016_IAH4','2018_IAH4','2014_IAH1','2015_IAH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_WIH2'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_WIH2',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_WIH2'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_OHH1','2018_OHH1'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_OHH1',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_OHH1','2018_OHH1'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_NYH3'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2018_NYH3',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_NYH3'),'imputed']='close_location_other_years'
soildata[soildata$Year_Exp%in%c('2015_NYH2','2014_NYH2'),c('X..Sand','X..Silt','X..Clay','Texture')]=soildata[soildata$Year_Exp=='2017_NYH2',c('X..Sand','X..Silt','X..Clay','Texture')]
soildata[soildata$Year_Exp%in%c('2015_NYH2','2014_NYH2'),'imputed']='close_location_other_years'



#Looking for missing soil information
soil_missing=soildata[is.na(soildata$X..Sand),]
s1=get_ssurgo(template=c('TX051'),label = 'TXH1')
daymet_tiles_data=daymet_tiles[]

#



###Add previous crop on the field (agronomic feature) or rotation info




##Check name cities Metadata file with real cities names from revgeo

##Remove inbred-related field trials weather information

################################
################################
##CREATE A DAILY WEATHER DATASET

