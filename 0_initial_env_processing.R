### Load useful packages

rm(list=ls())
library(elevatr)
library(rgdal)
library(tidyr)
library(readxl)
library(opencage)
library(intrval)
library(plyr)
library(dplyr)
library(USAboundaries)
library(ggplot2)
library(rnoaa)
library(stringr)
library(ggmap)
library(revgeo)
library(lubridate)
library(tidyverse)
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing")
`%notin%` <- Negate(`%in%`)


##Import 5 years of weather data
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/original_data_files")
weather2014=read.csv('g2f_2014_weather.csv',header = T,sep = ',',na.strings=c("","NA"))
weather2015=read.csv('g2f_2015_weather.csv',header = T,sep = ',',na.strings=c("","NA"))
weather2016=read.csv('g2f_2016_weather_data.csv',header = T,sep = ',',na.strings=c("","NA"))
weather2017=read.csv('g2f_2017_weather_data.csv',header = T,sep = ',',na.strings=c("","NA"))
weather2018=read.csv('g2f_2018_weather_clean.csv',header = T,sep = ',',na.strings=c("","NA","T"))[,1:28]

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
weather2016$Day.of.Year=as.numeric(as.vector(weather2016$Day.of.Year))
weather2017$Day.of.Year=paste0(weather2017$Month,'/',weather2017$Day,'/','17')
weather2017$Day.of.Year=as.POSIXct(as.Date(weather2017$Day.of.Year, "%m/%d/%Y"))
weather2017$Day.of.Year=difftime(weather2017$Day.of.Year,as.POSIXct(as.Date("01/01/17","%m/%d/%Y")),units='days')+1
weather2017$Day.of.Year=as.numeric(as.vector(weather2017$Day.of.Year))
weather2018$Day.of.Year=paste0(weather2018$Month,'/',weather2018$Day,'/','18')
weather2018$Day.of.Year=as.POSIXct(as.Date(weather2018$Day.of.Year, "%m/%d/%Y"))
weather2018$Day.of.Year=difftime(weather2018$Day.of.Year,as.POSIXct(as.Date("01/01/18","%m/%d/%Y")),units='days')+1
weather2018$Day.of.Year=as.numeric(as.vector(weather2018$Day.of.Year))

## Remove dulicated lines (present for 2017! --> NYH2)
weather2014<-weather2014[!duplicated(weather2014), ]
weather2015<-weather2015[!duplicated(weather2015), ]
weather2016<-weather2016[!duplicated(weather2016), ]
weather2017<-weather2017[!duplicated(weather2017), ]
weather2018<-weather2018[!duplicated(weather2018), ]


## Merging the annual weather files based on a subset of columns
names_col=c("Field.Location" ,"Station.ID","NWS.Network","Month","Day","Year",'Day.of.Year','Time.Local.','Temperature..C.','Dew.Point..C.','Relative.Humidity....',"Solar.Radiation..W.m2.","Rainfall..mm.",'Wind.Speed..m.s.','Wind.Gust..m.s.','Wind.Direction..degrees.' , "Photoperiod..hours.","Column.Altered" ,"Altered.Column.Names", "Cleaning.Method" ,'Comment')              
weather2014bis=weather2014[,names_col]
weather2015bis=weather2015[,names_col]
weather2016bis=weather2016[,names_col]
weather2017bis=weather2017[,names_col]
weather2018bis=weather2018[,names_col]


weather_all=rbind(weather2014bis,weather2015bis,weather2016bis,weather2017bis,weather2018bis)
weather_all <- separate_rows(weather_all,'Field.Location')
#write.table(weather_all,'weather_all.txt',sep = '\t',quote = F,col.names = T,row.names = F)
#rm(weather2014,weather2014bis,weather2015,weather2015bis,weather2016,weather2016bis,weather2017,weather2017bis,weather2018,weather2018bis)

##################################################################################
##################################################################################
####Loading the field metadata files: geographical coordinates and name cities etc

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

## Adapt col names and some experiment names
colnames(field2016)[1]='Experiment'
colnames(field2017)[1]='Experiment'
colnames(field2018)[1]='Experiment'
colnames(field2015)[c(19,20)]=c('lat','long')
colnames(field2016)[c(23,24)]=c('lat','long')
colnames(field2017)[c(26,27)]=c('lat','long')
colnames(field2018)[c(26,27)]=c('lat','long')
field2014$Experiment=as.character(as.vector(field2014$Experiment))
field2014[field2014$Experiment=="G2FWI-HYB",'Experiment']=as.character('WIH1')

#Retrieving Lon/Lat data for each field experiment --> use of the field corner (lower left) 
#For retrieving soil data: more accurate to use field than weather location, which can be very different
#If the Weather station had dubious data, then we will use a nearby location (<1.5 km) or interpolated data from weather stations in the area instead of field weather data.

## Note: correction needed for obvious mistakes in the field Lon/Lat columns to correct (case 1) + or data not reported in the field coordinates columns (only in the weather station columns) (case 2) + data completely absent for some Year_Exp: case 3 (if city indicated --> use of coordinates from previous year on corresponding experiment)


field2017[field2017$Experiment=='TXH1',c('lat','long')]=c(30.54535,	-96.43258) # case 1
field2017[field2017$Experiment=='MNH1',c('lat','long')]=c(44.06981,-93.5338) # case 1
field2018[field2018$Experiment=='IAH1',c('lat','long')]=c(41.19870,	-91.48618) # case 1
field2018[field2018$Experiment=='IAH2',c('lat','long')]=c(42.06593,	-94.72742) # case 1
field2018[field2018$Experiment=='IAH3',c('lat','long')]=c(41.98745,	-92.26014) # case 1
field2018[field2018$Experiment=='IAH4',c('lat','long')]=c(41.99775,	-93.69625) # case 1
field2018[field2018$Experiment=='WIH2',c('lat','long')]=c(43.30480,	-89.38520) # case 1: mistakes with some field coordinates corners (but not all: probable correct ones used) with 43.40490 as latitude impossible --> located in Dekorra (not Arlington as specified) + would be in forest!

field2016[field2016$Experiment=='IAH1',c('lat','long')]=c(41.19869,	-91.48620) # case 2
field2015[field2015$Experiment=='IAH2',c('lat','long')]=c(42.0675, -93.618 ) # case 2
field2016[field2016$Experiment=='IAH2',c('lat','long')]=c(42.06591,	-94.72745) # case 2
field2016[field2016$Experiment=='IAH3',c('lat','long')]=c(41.98738,	-92.26016) # case 2
field2016[field2016$Experiment=='IAH4',c('lat','long')]=c(41.99750,	-93.69990) # case 2

field2016[field2016$Experiment=='NEH1',c('lat','long')]=c(41.16636,	-96.41726) #Case 3: NEH1 2016 (no coord. data): assume to be the same field as 2015
field2016[field2016$Experiment=='NEH4',c('lat','long')]=c(41.16702,	-96.41719) #Case 3: NEH4 2016 (no coord. data): assume to be the same field as 2015
field2016[field2016$Experiment=='SCH1',c('lat','long')]=c(34.62261,	-82.73796) #Case 3: SCH1 2016 (no coord. data): assume to be the same field as 2017
field2017[field2017$Experiment=='ILH1',c('lat','long')]=c(40.06119,	-88.23327) #Case 3: ILH1 2017 (no coord. data): assume to be the same field as 2016
field2017[field2017$Experiment=='INH1',c('lat','long')]=c(40.47835,	-86.99013) #Case 3: INH1 2017 (no coord. data): assume to be the same field as 2016


## 2 reps for MOH1 next to each other --> take only 1 field Lon/Lat as reference for the two reps
field2018[field2018$Experiment=='MOH1- rep 1','Experiment']='MOH1'
field2018=field2018[-which(field2018$Experiment=='MOH1- rep 2'),]
field2018[field2018$Experiment=='MOH1',c('lat','long')]=c(38.89852,-92.20918)


## Binding all field data from 5 years
geofield2014=field2014[,c('lat','long','Year','Experiment','City')]
geofield2015=field2015[,c('lat','long','Year','Experiment','City')]
geofield2016=field2016[,c('lat','long','Year','Experiment','City')]
geofield2017=field2017[,c('lat','long','Year','Experiment','City')]
geofield2018=field2018[,c('lat','long','Year','Experiment','City')]
geo_data_fields=rbind(geofield2014,geofield2015,geofield2016,geofield2017,geofield2018)
colnames(geo_data_fields)[4]='Field.Location'



###################ELIMINATE SOME EXPERIMENTS###################################

## Eliminate those corresponding to inbred lines experiments and disease trials
## NYH1 disease trial with no info + TXH2 no info except for 2014 --> need to be removed, absolutely no info on this trials
## removed because no phenotypic data present in the final hybrid pheno files: 2015_ILH2, 2015_IAH1, 2015_IAH2, 2015_IAH3, 2015_IAH4, 2017_ILH2, 2017_NEH1, 2017_NEH2, 2017_SCH1 (no yield data) 


to_remove_inbred_trials=c("AZI1", "AZI2", "DEI1", "GAI1", "GAI2","IAI1", "IAI2","SDI1", "IAI3", "IAI4","ILI1", "INI1","KSI1" ,"MNI1" ,"MNI2","MOI1", "MOI2","MOI3", "NCI1" ,"NEI1","NYI1","NYI2", "PAI1" ,"PAI2" ,"TXI1" ,"TXI2","TXI3" , "WII1", "WII2" )
geo_data_fields=geo_data_fields[-which(geo_data_fields$Field.Location%in%to_remove_inbred_trials),]

geo_data_fields=geo_data_fields[-which(geo_data_fields$Field.Location=='NYH1'&geo_data_fields$Year%in%c(2015,2016,2017)),]

geo_data_fields=geo_data_fields[-which(geo_data_fields$Field.Location=='TXH2'&geo_data_fields$Year%in%c(2015,2016,2017,2018)),]

geo_data_fields=geo_data_fields[-which(geo_data_fields$Field.Location%in%c("ILH2","IAH1","IAH2","IAH3","IAH4")&geo_data_fields$Year%in%c(2015)),]
geo_data_fields=geo_data_fields[-which(geo_data_fields$Field.Location%in%c("ILH2","NEH1","NEH2","SCH1")&geo_data_fields$Year%in%c(2017)),]


geo_data_fields$lat=as.numeric(as.vector(geo_data_fields$lat))
geo_data_fields$long=as.numeric(as.vector(geo_data_fields$long))

#################################################
#################################################

# Use opencage to add county, state and city
# Check that coordinates correspond to meta information.
# 2014 DEH1 OK ? --> coordinates look weird in comparison to following years, but no modifications done.


state=vector()
county=vector()
city=vector()
village=vector()
hamlet=vector()
for (j in 1:nrow(geo_data_fields)) {
  output <- opencage_reverse(latitude = geo_data_fields$lat[j], 
                              longitude = geo_data_fields$long[j],key='144c9a0f9b394bd78730e12d3127fe8c')
  state[j]<-as.character(output$results$components.state)
  county[j]<-as.character(output$results$components.county)
  tryCatch({city[j]<-as.character(output$results$components.city)},warning=function(e) NULL,error = function(e)NULL)
  tryCatch({village[j]<-as.character(output$results$components.village)},warning=function(e) NULL,error = function(e)NULL)
  tryCatch({hamlet[j]<-as.character(output$results$components.hamlet)},warning=function(e) NULL,error = function(e)NULL)
  
  
}


geo_data_fields$state=state
geo_data_fields$county=county
geo_data_fields$city_village_reversegeo=village
geo_data_fields$city_village_reversegeo[is.na(geo_data_fields$city_village_reversegeo)]=city[is.na(geo_data_fields$city_village_reversegeo)]
geo_data_fields$city_village_reversegeo[is.na(geo_data_fields$city_village_reversegeo)]=hamlet[is.na(geo_data_fields$city_village_reversegeo)]
geo_data_fields$city_village_reversegeo[is.na(geo_data_fields$city_village_reversegeo)]=str_remove_all(county[is.na(geo_data_fields$city_village_reversegeo)], " County")


##Add beginning and end of the growing season to geo_data_fields table

growingseason=read.table("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/original_data_files/growingseason_extendeddates.txt",sep = '\t',header = T)
growingseason$Date.Planted1=as.POSIXct(as.Date(growingseason$Date.Planted, "%m/%d/%Y"))
growingseason$Date.Harvested1=as.POSIXct(as.Date(growingseason$Date.Harvested, "%m/%d/%Y"))
growingseason$Date.Planted2=yday(growingseason$Date.Planted1)
growingseason$Date.Harvested2=yday(growingseason$Date.Harvested1)


## Binding geo_data_fields with info on Planted/Harvest dates

geo_data_fields$Year_Exp=paste0(geo_data_fields$Year,'_',geo_data_fields$Field.Location)
geo_data_fields=merge(geo_data_fields,growingseason[,c('Year_Exp',"Date.Planted2","Date.Harvested2")],by='Year_Exp',all.x=T)
colnames(geo_data_fields)[which(colnames(geo_data_fields)=='Date.Planted2')]='Date.Planted'
colnames(geo_data_fields)[which(colnames(geo_data_fields)=='Date.Harvested2')]='Date.Harvested'

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/")






##############################################################
################## SOIL INFORMATION ##########################
#######Info very limited for years 2014 + 2015 --> use of USDA Web SOil Survey data
###Load the soil data files to add the soil texture (at least)
###Prepare the data.frames with common column names
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/original_data_files")

soil2014=as.data.frame(cbind('Year'=2014,'Experiment'=as.character(field2014$Experiment),'Soil.test.type'=as.character(as.vector(field2014$Soil.test.type)),'Soil.texture'=as.character(as.vector(field2014$Soil.texture)),'Soil.pH'=as.character(as.vector(field2014$Soil.pH))))
soil2014$Year_Exp=paste(soil2014$Year,soil2014$Experiment,sep = '_')
colnames(soil2014)[4]='Texture'
soil2014$X..Sand=NA
soil2014$X..Silt=NA
soil2014$X..Clay=NA
soil2014$OM=NA
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


# Bind soil information across years
soil_1=rbind(soil2014[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2015[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2016[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2017[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')],soil2018[,c('Year_Exp','X..Sand','X..Silt','X..Clay','Texture','OM')])
soil_1=soil_1[-which(duplicated(soil_1$Year_Exp)),]
geo_data_fields=merge(geo_data_fields,soil_1,by='Year_Exp',all.x=T)


#Create column to indicate that the data were imputed from other years of experiments
geo_data_fields$imputed=NA
geo_data_fields[!is.na(geo_data_fields$X..Sand),'imputed']='NO'


######VERY CLOSE COORDINATES --> SAME SOIL COMPOSITION ASSUMED######
geo_data_fields[geo_data_fields$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay')]=c(11,30,59)

geo_data_fields[geo_data_fields$Year_Exp=='2015_ONH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_ONH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp=='2015_ONH1','imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp=='2018_DEH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_DEH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp=='2018_DEH1','imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_GAH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_GAH1','2015_GAH1','2018_GAH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2016_IAH1','2018_IAH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2017_IAH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2016_IAH1','2018_IAH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_IAH3'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2017_IAH3',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_IAH3'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_MOH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2018_MOH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_MOH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_MOH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_MOH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_MOH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_ILH1','2015_ILH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_ILH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_ILH1','2015_ILH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_MNH1','2015_MNH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_MNH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_MNH1','2015_MNH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2017_NCH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_NCH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_NCH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2016_NCH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2014_NCH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1','2018_TXH1- Dry','2018_TXH1- Early','2018_TXH1- Late'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2014_TXH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_TXH1','2016_TXH1','2017_TXH1','2018_TXH1- Dry','2018_TXH1- Early','2018_TXH1- Late'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_OHH1'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2017_OHH1',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_OHH1'),'imputed']='close_location_other_years'
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_NYH2','2014_NYH2'),c('X..Sand','X..Silt','X..Clay','Texture','OM')]=geo_data_fields[geo_data_fields$Year_Exp=='2017_NYH2',c('X..Sand','X..Silt','X..Clay','Texture','OM')]
geo_data_fields[geo_data_fields$Year_Exp%in%c('2015_NYH2','2014_NYH2'),'imputed']='close_location_other_years'


#Looking for missing soil information
# Note_ for 2014_NEH1: soil data based on estimation given the texture, because only location for which no web soil survey data were available.
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/original_data_files")

soil_imputed=as.data.frame(read.table('soil_imputed.txt',header=T,sep="\t"))[,c(1:8)]
soil_imputed$imputed='web_soil_survey'
geo_data_fields$Year_Exp=as.character(as.vector(geo_data_fields$Year_Exp))
geo_data_fields$Texture=as.character(as.vector(geo_data_fields$Texture))
soil_imputed$Year_Exp=as.character(as.vector(soil_imputed$Year_Exp))
soil_imputed$Texture=as.character(as.vector(soil_imputed$Texture))

geo_data_fields[which(is.na(geo_data_fields$X..Sand)),c('Year_Exp','X..Sand','X..Silt','X..Clay','OM','Texture','imputed')]=soil_imputed[match(geo_data_fields[which(is.na(geo_data_fields$X..Sand)),'Year_Exp'],soil_imputed$Year_Exp),c('Year_Exp','Sand','Silt','Clay','OrganicMatter','Texture','imputed')]





###################AGRONOMIC MANAGEMENT############################
###Add previous crop on the field (agronomic feature) and total N
agro_features=c('Year_Exp','Previous.crop')
field2014$Year_Exp=paste(field2014$Year,field2014$Experiment,sep = '_')
field2015$Year_Exp=paste(field2015$Year,field2015$Experiment,sep = '_')
field2016$Year_Exp=paste(field2016$Year,field2016$Experiment,sep = '_')
field2017$Year_Exp=paste(field2017$Year,field2017$Experiment,sep = '_')
field2018$Year_Exp=paste(field2018$Year,field2018$Experiment,sep = '_')
field2018[field2018$Year_Exp=='2018_MOH1- rep 1'|field2018$Year_Exp=='2018_MOH1- rep 2','Year_Exp']<-'2018_MOH1'
colnames(field2015)[c(10)]=c('Previous.crop')
colnames(field2016)[c(13)]=c('Previous.crop')
colnames(field2017)[c(16)]=c('Previous.crop')
colnames(field2018)[c(16)]=c('Previous.crop')

agro_management=rbind(field2014[,agro_features],field2015[,agro_features],field2016[,agro_features],field2017[,agro_features],field2018[,agro_features])
agro_management$Year_Exp=as.character(as.vector(agro_management$Year_Exp))
geo_data_fields=merge(geo_data_fields,agro_management,by='Year_Exp',all.x=T)
geo_data_fields$Previous.crop=as.character(as.vector(geo_data_fields$Previous.crop))

for (j in 1:nrow(geo_data_fields)) {
  if (geo_data_fields[j,'Previous.crop']%in%c('soybean','soybeans','Soybeans','Soybean')){geo_data_fields[j,'Previous.crop']<-'soybean'}
  if (geo_data_fields[j,'Previous.crop']%in%c('wheat/soybean double crop','Wheat and Double Crop Soybeans','Small Grains and Double Crop Soybeans')){geo_data_fields[j,'Previous.crop']<-'wheat/double crop soybean'}
  if (geo_data_fields[j,'Previous.crop']%in%c('Cotton')){geo_data_fields[j,'Previous.crop']<-'cotton'}
  if (geo_data_fields[j,'Previous.crop']%in%c('Corn','corn','Corn ')){geo_data_fields[j,'Previous.crop']<-'corn'}
  if (geo_data_fields[j,'Previous.crop']%in%c('Wheat','Winter Wheat','Fallow most of 2014 winter planted in fall of 2014 then sprayed with Glystar 24 floz/a on 5/3/15  and killed spring of 2015 spray ')){geo_data_fields[j,'Previous.crop']<-'wheat'}
  if (geo_data_fields[j,'Previous.crop']%in%c('Lima beans followed by rye cover crop')){geo_data_fields[j,'Previous.crop']<-'lima beans'}
  if (geo_data_fields[j,'Previous.crop']%in%c('Sorghum')){geo_data_fields[j,'Previous.crop']<-'sorghum'}
}
geo_data_fields[which(geo_data_fields$Previous.crop==''),'Previous.crop']<-NA



################### ELEVATION ############################

library(elevatr)
loc=geo_data_fields[,c('long','lat')]

elev<-get_elev_point(locations=loc,units='meters',prj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
elev_df=cbind(loc,elev@data)
geo_data_fields$elev=elev_df[,'elevation']
geo_data_fields=arrange(geo_data_fields,Year_Exp)

## Need to add manually for Canadian loc.
geo_data_fields[geo_data_fields$Year_Exp=='2014_ONH1',c('elev')]<-336
geo_data_fields[geo_data_fields$Year_Exp=='2014_ONH2',c('elev')]<-204
geo_data_fields[geo_data_fields$Year_Exp=='2015_ONH1',c('elev')]<-328
geo_data_fields[geo_data_fields$Year_Exp=='2015_ONH2',c('elev')]<-204
geo_data_fields[geo_data_fields$Year_Exp=='2016_ONH1',c('elev')]<-329
geo_data_fields[geo_data_fields$Year_Exp=='2016_ONH2',c('elev')]<-203
geo_data_fields[geo_data_fields$Year_Exp=='2017_ONH1',c('elev')]<-334
geo_data_fields[geo_data_fields$Year_Exp=='2017_ONH2',c('elev')]<-203
geo_data_fields[geo_data_fields$Year_Exp=='2018_ONH2',c('elev')]<-203

write.table(geo_data_fields,'/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/geo_data_fields.txt',col.names = T,row.names = F,quote = F,sep = '\t')

####################################################################################################################
#### Create daily data frame for each Year_Exp present in geo_data_fields between sowing date and harvest dates ####
####################################################################################################################


rep.row <- function(r, n){
  colwise(function(x) rep(x, n))(r)
}

dd=list()
s=1
for (v in geo_data_fields$Year_Exp) {
  dd[[s]]=as.data.frame(matrix(NA,ncol = 2,nrow = length(geo_data_fields[geo_data_fields$Year_Exp==v,"Date.Planted"]:geo_data_fields[geo_data_fields$Year_Exp==v,"Date.Harvested"])))
  dd[[s]][,2]=v
  dd[[s]][,1]=geo_data_fields[geo_data_fields$Year_Exp==v,"Date.Planted"]:geo_data_fields[geo_data_fields$Year_Exp==v,"Date.Harvested"]
  df=geo_data_fields[geo_data_fields$Year_Exp==v,2:ncol(geo_data_fields)]
  dd[[s]]=cbind(dd[[s]],df[rep(seq_len(nrow(df)),each=nrow(dd[[s]])),])
  s=s+1
}

df<-plyr::ldply(dd,data.frame)
colnames(df)[c(1,2)]<-c('Day.of.Year','Year_Exp')
daily_weather=df
rm(df,dd)

# Add months and day
daily_weather$month <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%m'))
daily_weather$day <- with(daily_weather, format(strptime(paste(Year, Day.of.Year), format = "%Y %j"), '%d'))

## Order daily data frame
daily_weather <-arrange(daily_weather,Year_Exp,Day.of.Year)
write.table(daily_weather,'/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/daily_weather_afterstep0.txt',col.names = T,row.names = F,sep = '\t',quote = F)



# Reduce weather data to Year_Exp remaining in geo_data_fields
weather_all$Year_Exp=paste0(weather_all$Year,'_',weather_all$Field.Location)
weather=weather_all[weather_all$Year_Exp%in%geo_data_fields$Year_Exp,]


## Add same weather_all data for 2016_NYH2 and 2016_NYH3, and for 2017_NYH2 and 2017_NYH3
to_add1<-weather_all[weather_all$Year_Exp=='2016_NYH2',]
to_add1$Year_Exp<-'2016_NYH3'
to_add1$Field.Location<-'NYH3'
weather_all<-rbind(weather_all,to_add1)

to_add2<-weather_all[weather_all$Year_Exp=='2017_NYH2',]
to_add2$Year_Exp<-'2017_NYH3'
to_add2$Field.Location<-'NYH3'
weather_all<-rbind(weather_all,to_add2)


## Order weather data frame
weather <-arrange(weather,Year_Exp,Day.of.Year,Time.Local.)






write.table(weather,'/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/weather_intermediate.txt',col.names = T,row.names = F,sep = '\t',quote = F)














