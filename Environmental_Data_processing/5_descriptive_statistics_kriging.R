rm(list = ls())


## Load or install packages

source('fahrenheit_to_celsius.R')

packages = c("xts","data.table","rgdal","readxl","opencage","intrval","dplyr","USAboundaries","ggplot2","rnoaa","stringr","ggmap","revgeo","lubridate")


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
# Load dataset with missing values #
# ------------------------------------------------------------------------------
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing"
)
daily_weather = read.table(
  '2_merged_dataset_before_interpolation.txt',
  header = T,
  sep = '\t',
  na.strings = c(NA, ''),
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)
daily_weather = plyr::arrange(daily_weather, Year_Exp, Day.of.Year)

row.names(daily_weather) = NULL
s = nrow(daily_weather)
set_locations=unique(daily_weather$Year_Exp)




# ------------------------------------------------------------------------------
# Solar radiation data
# ------------------------------------------------------------------------------

# 1) Reading the NASA files
setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/NASA"
)
NASA_FILES <- list()
for (i in list.files(pattern = ".txt")) {
  NASA_FILES[[i]] <- read.table(i, header = T,sep = '\t')
}

NASA_FILES <- do.call('rbind', NASA_FILES)
NASA_FILES$day <- lubridate::yday(NASA_FILES$date)
NASA_FILES[NASA_FILES$ALLSKY_SFC_SW_DWN<0,'ALLSKY_SFC_SW_DWN']<-NA
NASA_FILES = plyr::arrange(NASA_FILES, Year_Exp, day)


# 2) Correlation between field solar radiation values and NASA files
cor <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
names(cor)<-unique(NASA_FILES$Year_Exp)
na.count <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
names(na.count)<-unique(NASA_FILES$Year_Exp)
total_field_station <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))
total_NASA_data <-
  vector(mode = 'numeric', length = length(unique(NASA_FILES$Year_Exp)))


for (j in 1:length(unique(NASA_FILES$Year_Exp))) {
  var = unique(NASA_FILES$Year_Exp)[j]
  sub_data=NASA_FILES[NASA_FILES$Year_Exp == var, ]
  matched_values = sub_data[match(daily_weather[daily_weather$Year_Exp ==
                                                  var, 'Day.of.Year'], sub_data$day), 'ALLSKY_SFC_SW_DWN']
  matched_table = sub_data[match(daily_weather[daily_weather$Year_Exp ==
                                                  var, 'Day.of.Year'], sub_data$day), c('date','ALLSKY_SFC_SW_DWN')]
  write.table(matched_table,paste0('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/SOLAR_RADIATION/NASA_',var,'.txt'),sep = '\t',quote = F,col.names = T,row.names = F)
  tryCatch({
    cor[j] <- cor(matched_values,
                  daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'],
                  method = 'pearson',
                  use = 'complete.obs')
  }, error = function(e) 
    return(as.numeric(NA)))
  
  na.count[j] <-
    length(which(is.na(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2']))) /
    length(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'])
  
  total_field_station[j] <- sum(as.numeric(as.vector(daily_weather[daily_weather$Year_Exp == var, 'incoming_radiation_MJm2'])),na.rm=T)
  total_NASA_data[j] <- sum(matched_values,na.rm = T)
  
}

# 3) Generate data.frame with corr

setwd(
  '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/summary_NASA'
)
cor_table = as.data.frame(cbind(names(cor), na.count, cor,total_field_station,total_NASA_data))
colnames(cor_table) <-
  c('Year_Exp',
    'flagged.NA.percent.growing.season.field.data',
    'pearson.corr','total_field_station_(MJ/m2/day)','total_NASA_data_(MJ/m2/day)')

cor_table$total_NASA_data<-as.numeric(as.vector(cor_table$total_NASA_data))
cor_table$total_field_station<-as.numeric(as.vector(cor_table$total_field_station))

cor_table$difference <-
  cor_table$total_field_station - cor_table$total_NASA_data
cor_table$difference<-as.numeric(as.vector(cor_table$difference))
cor_table$Replacement<-NA
cor_table[abs(cor_table$difference)>300,'Replacement']<-'Total'
cor_table[abs(cor_table$difference)< 0.6,'Replacement']<-'Total'
cor_table[cor_table$Replacement%notin%'Total','Replacement']<-'Partial'




write.table(
  as.data.frame(cor_table),
  'cor.daily.solar.radiation.fieldstation.vs.NASA.txt',
  col.names = T,
  row.names = F,
  quote = F,
  sep = '\t'
)



# ------------------------------------------------------------------------------
# Metorological variables: 'TMAX', 'TMIN'
# ------------------------------------------------------------------------------

for (meteo_variable in c('TMAX', 'TMIN')) {
  print(paste0(meteo_variable,' : start'))
  # 1) Reading the imputed data using kriging and constructing summarized dataset from kriging results before imputation
  
  meteo_in_table = meteo_variable
  
  
  setwd(
    paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '/',
      sep = ''
    )
  )
  list_files1 = list.files(pattern = ".RDS")
  
  
  
  summary = as.data.frame(matrix(ncol = 13, nrow = length(list.files(pattern = ".RDS"))))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'model.spatiotemporal.variogram',
    'pearson.cor.5fold.CV.model',
    'rmse.5fold.CV.model',
    'Pearson.correlation.daily.interpolated.with.field.obs',
    'RMSE.interpolated.with.field.obs',
    'weather.network.used',
    'interpolation_method',
    'near_or_interpolated'
   )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  summary$interpolation_method = 'kriging'
  predictions <- list()
  
  n = 1
  for (i in list_files1) {
    print(i)
    
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
        meteo_variable,
        '/',
        sep = ''
      )
    )
    v = str_sub(i, 1, nchar(i) - 4)
    if(v%notin%set_locations){
      print(v)
      next}
    summary[n, 'Year_Exp'] <- v
    summary[n, 'flagged.NA.percent.growing.season.field.data'] <-
      length(which(is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable]))) /
      length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])
    
   
    if (length(which(!is.na(daily_weather[daily_weather$Year_Exp == v, meteo_variable])))>length(daily_weather[daily_weather$Year_Exp == v, meteo_variable])/2){
      summary[n, 'Pearson.correlation.daily.interpolated.with.field.obs'] <- cor(as.numeric(as.vector(readRDS(i)[[1]][, 2])),
                           daily_weather[daily_weather$Year_Exp == v, meteo_variable],
                           method = 'pearson',
                           use = 'complete.obs')
      summary[n,'RMSE.interpolated.with.field.obs'] <- sqrt( mean( (as.numeric(as.vector(readRDS(i)[[1]][, 2]))- daily_weather[daily_weather$Year_Exp == v, meteo_variable])^2 , na.rm = TRUE ) )
        
                                                                                  }
    else{ 
      summary[n, 'Pearson.correlation.daily.interpolated.with.field.obs'] <- NA
      summary[n,'RMSE.interpolated.with.field.obs'] <-NA
      }
      
      
      summary[n, 'pearson.cor.5fold.CV.model'] <-
        max(unlist(lapply(readRDS(i)[['5f.cv.kriging.cor']], mean)), na.rm = t)
      summary[n, 'rmse.5fold.CV.model'] <-
        min(unlist(lapply(readRDS(i)[['5f.cv.kriging.rmse']], mean)), na.rm = t)
      
      summary[n, 'model.spatiotemporal.variogram'] <-
        names(which.min(unlist(lapply(
          readRDS(i)[['5f.cv.kriging.rmse']], mean
        ))))
    
    
    values_interpolated=as.numeric(as.vector(readRDS(i)[[1]][, 2]))
    dates_v=as.vector(readRDS(i)[[1]][, 3])
    selected_data<-as.data.frame(cbind(v,values_interpolated,dates_v))
    #write.table(as.data.frame(cbind(v,values_interpolated,dates_v)),file=paste(meteo_variable,'_selected_model_',v,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    
      
    
    setwd(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files'
    )
    pat = paste(meteo_variable, '_', v, '.txt', sep = '')
    pat2 = paste(v, '_weather.txt', sep = '')
    
    if (length(list.files(pattern = pat, recursive = TRUE)) == 0 &
        length(list.files(pattern = pat2, recursive = TRUE)) > 0) {
      #if not a GHCND formatted file nor a GSOD file, it must be CANADA weather files.
      summary[n, 'weather.network.used'] <- 'WeatherCAN'
      s = list.files(pattern = pat2, recursive = T)
      summary[n, 'nb.weather.stations.used'] <-
        length(unique(read.table(
          s, header = T, sep = '\t'
        )[, 'station_name']))
      
    }
    
    else{
      s = list.files(pattern = pat, recursive = TRUE)
      summary[n, 'weather.network.used'] <- 'GHCND'
      summary[n, 'nb.weather.stations.used'] <-
        length(unique(read.table(
          s, header = T, sep = '\t'
        )[, 'id']))
    }
    setwd(
      paste(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
        meteo_variable,
        '/',
        sep = ""
      )
    )
    
    ## If a weather station is located <1.5 km nearby, use of that station instead of interpolated data
    
    
    setwd(
      paste0(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
        meteo_variable,
        '_nearest/'
      )
    )
    nearest<-readRDS(paste0(v,'_1station.RDS'))
    
    if(nearest$min_dist<2){
      near_data<-as.data.frame(nearest$predictions_YearExp[,c('Year_Exp',meteo_variable,'date')])
      summary$near_or_interpolated[n]<-'near'
      write.table(near_data,file=paste('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/',meteo_variable,'/selected_model_',v,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }
    else{
      colnames(selected_data)<-c('Year_Exp','TMAX','date')
      summary$near_or_interpolated[n]<-'interpolated'
      write.table(selected_data,file=paste('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/',meteo_variable,'/selected_model_',v,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      
    }
    
    n <- n + 1
  }
  
  summary$Replacement<-NA
  summary[which(summary$Pearson.correlation.daily.interpolated.with.field.obs<0.65),'Replacement']<-'Total'
  summary[which(summary$flagged.NA.percent.growing.season.field.data>0.5),'Replacement']<-'Total'
  summary[summary$Replacement%notin%'Total','Replacement']<-'Partial'
  
  
  
  
  write.table(
    summary,
    file = paste(
      "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
      meteo_variable,
      '_summary.txt',
      sep = ''
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep = '\t'
  )
}

