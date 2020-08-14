rm(list = ls())


## Load or install packages
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing")
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
# Meteorological variables: 'TMAX', 'TMIN'
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
  
  
  
  
  summary = as.data.frame(matrix(ncol = 13, nrow = length(set_locations)))
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
    'near_or_interpolated_selected'
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  summary$interpolation_method = 'kriging'
  predictions <- list()
  
  n = 1
  for (i in paste0(set_locations,'.RDS')) {
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
    
    ## If a weather station is located <2 km nearby and has complete data, use of that station instead of interpolated data
    
    
    setwd(
      paste0(
        "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
        meteo_variable,
        '_nearest/'
      )
    )
    nearest<-readRDS(paste0(v,'_1station.RDS'))
    
    if(nearest$min_dist<2&length(which(is.na(nearest$predictions_YearExp[,c(meteo_variable)])))==0){
      near_data<-as.data.frame(nearest$predictions_YearExp[,c('Year_Exp',meteo_variable,'date')])
      summary$near_or_interpolated_selected[n]<-'near'
      write.table(near_data,file=paste('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/',meteo_variable,'/selected_model_',v,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    }
    else{
      colnames(selected_data)<-c('Year_Exp',meteo_variable,'date')
      summary$near_or_interpolated_selected[n]<-'interpolated'
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




# ------------------------------------------------------------------------------
# Meteorological variables: 'PRCP'
# ------------------------------------------------------------------------------

setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/')

Winning_kriging_model = as.data.frame(matrix(NA, nrow = length(set_locations), ncol = 19))

colnames(Winning_kriging_model) = c(
  'Year_Exp',
  'Number.days.growing.season',
  'Total_field_after_processing(mm)',
  '%NA.field.season',
  'radius',
  'interpolation_method',
  'nb_stations_used_for_kriging',
  'Winning_kriging_model',
  'Total.winning.model(mm)',
  'Total.with.elevation(mm)',
  'RMSE.model.with.elevation',
  'Pearson.Cor.model.with.elevation',
  'Total.no.elevation(mm)',
  'RMSE.model.no.elevation',
  'Pearson.Cor.model.no.elevation',
  'weather.network.used',
  'Total.near.station(mm)',
  'Name.near.station',
  'comment'
)


Winning_kriging_model$interpolation_method<-'kriging'
Winning_kriging_model[,c(3,4,9,10,11,12,13,14,15,17)]<-apply(Winning_kriging_model[, c(3,4,9,10,11,12,13,14,15,17)], 2,as.numeric)

n = 1
for (s in set_locations) {
  
  print(s)
  
  
  Winning_kriging_model[n, 'Year_Exp'] <- s
  Winning_kriging_model[n, 'radius'] <-70
  if (s%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH1','2018_ONH2')){
    Winning_kriging_model[n, 'weather.network.used'] <-'WeatherCAN'}
  else{Winning_kriging_model[n, 'weather.network.used'] <-'GHCND'}
  
  
  Winning_kriging_model[n, 'Number.days.growing.season'] <-length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  Winning_kriging_model[n, 'Total_field_after_processing(mm)']<-sum(daily_weather[daily_weather$Year_Exp==s,'PRCP'],na.rm = T)
  Winning_kriging_model[n, '%NA.field.season']<-length(which(is.na(daily_weather[daily_weather$Year_Exp==s,'PRCP'])))/length(daily_weather[daily_weather$Year_Exp==s,'PRCP'])
  
  
  setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_nearest')
  
  
  vector_nearest_station<-readRDS(paste(s,'_1station.RDS',sep = ''))
  Winning_kriging_model[n, 'Total.near.station(mm)']<-sum(as.numeric(as.vector(vector_nearest_station[[1]][,2])),na.rm = T)
  Winning_kriging_model[n,  'Distance.near.station']<-as.numeric(as.character(vector_nearest_station[['min_dist']]))
  Winning_kriging_model[n, 'Name.near.station']<-vector_nearest_station[['min_dist_station_name']]
  
  
  
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL")
  
  if (!file.exists(paste('PRCP_elev',s,'.RDS',sep = ''))){
    print(paste0(s,': kriging could not finished due to size data -- we use nearest station'))
    Winning_kriging_model[n, 'Winning_kriging_model'] <- 'nearest_station'
    Winning_kriging_model[n, 'Total.winning.model(mm)']<-Winning_kriging_model[n, 'Total.near.station(mm)']
    Winning_kriging_model[n, 'comment']<-'kriging_never_finished'
    write.table(vector_nearest_station[[1]][,c('Year_Exp', 'PRCP', 'date')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    rm(selected_data)
    setwd(paste0("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files/",year))
    if (s%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH1','2018_ONH2')){
      dat=read.table(paste0(s,'_weather.txt'),header=T,sep='\t')
      Winning_kriging_model[n, 'nb_stations_used_for_kriging'] <-length(unique(dat[which(!is.na(dat$prcp)),'station_name']))}
    else{
      Winning_kriging_model[n, 'nb_stations_used_for_kriging']<-length(unique(read.table(paste0('PRCP_',s,'.txt'),header=T,sep='\t')$id))}
    
    n <- n + 1
    next
  }
  
  elev_model<-readRDS(paste('PRCP_elev',s,'.RDS',sep = ''))[[1]]
  
  elev_model<-cbind(elev_model,vector_nearest_station)
  
  if (!file.exists(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))|!file.exists(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))){
    Winning_kriging_model[n, 'comment']<-'missing LOOCV results'
    if (file.exists(paste('PRCPno_elev',s,'.RDS',sep = ''))){
    print(paste0(s,': no results regarding LOOCV here '))
    no_elev_model<-readRDS(paste('PRCPno_elev',s,'.RDS',sep = ''))[[1]]
    no_elev_model<-cbind(no_elev_model,vector_nearest_station)
    Winning_kriging_model[n, 'Total.no.elevation(mm)'] <-sum(no_elev_model[,LOOCV_min_RMSE_noelev])
    
    Winning_kriging_model[n, 'Winning_kriging_model'] <- paste('no_elevation_','summetric_Vgm3',sep='')
    Winning_kriging_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm3'],na.rm=T)
    selected_data<-no_elev_model[,c('predictions.table','summetric_Vgm3','as.character(dates)')]
    colnames(selected_data)<-c('Year_Exp','Kriging.data','dates')
    
    if(Winning_kriging_model[n,'Distance.near.station']<2&length(which(is.na(vector_nearest_station[[1]][,c('PRCP')])))==0){
      Winning_kriging_model[n, 'Winning_kriging_model'] <- 'nearest_station'
      Winning_kriging_model[n, 'Total.winning.model(mm)']<-Winning_kriging_model[n, 'Total.near.station(mm)']
      write.table(vector_nearest_station[[1]][,c('Year_Exp', 'PRCP', 'date')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    }
    else{write.table(selected_data,file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)}
    year=str_sub(s,1,4)
    rm(selected_data)
    setwd(paste0("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files/",year))
    if (s%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH1','2018_ONH2')){
      dat=read.table(paste0(s,'_weather.txt'),header=T,sep='\t')
      Winning_kriging_model[n, 'nb_stations_used_for_kriging'] <-length(unique(dat[which(!is.na(dat$prcp)),'station_name']))}
    else{
      Winning_kriging_model[n, 'nb_stations_used_for_kriging']<-length(unique(read.table(paste0('PRCP_',s,'.txt'),header=T,sep='\t')$id))}
    
    n <- n + 1
    next}
  } else{
    LOOCV<-readRDS(paste('elev_results_LOOCV_',s,'.RDS',sep = ''))
    
    
    LOOCV_min_RMSE_elev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    LOOCV_max_cor_elev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    
    Winning_kriging_model[n, 'Total.with.elevation(mm)'] <-sum(elev_model[,LOOCV_min_RMSE_elev])
    Winning_kriging_model[n, 'RMSE.model.with.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    Winning_kriging_model[n, 'Pearson.Cor.model.with.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    
    no_elev_model<-readRDS(paste('PRCPno_elev',s,'.RDS',sep = ''))[[1]]
    no_elev_model<-cbind(no_elev_model,vector_nearest_station)
    
    
    LOOCV<-readRDS(paste('no_elev_results_LOOCV_',s,'.RDS',sep = ''))
    names(LOOCV[['kriging_rmse']])[names(LOOCV[['kriging_rmse']])=='summetric_Vgm']<-'summetric_Vgm1'
    names(LOOCV[['kriging_cor']])[names(LOOCV[['kriging_cor']])=='summetric_Vgm']<-'summetric_Vgm1'
    LOOCV_min_RMSE_noelev<- names(which.min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    LOOCV_max_cor_noelev<- names(which.max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T)))))
    
    Winning_kriging_model[n, 'Total.no.elevation(mm)'] <-sum(no_elev_model[,LOOCV_min_RMSE_noelev])
    Winning_kriging_model[n, 'RMSE.model.no.elevation'] <-min(unlist(lapply(LOOCV[['kriging_rmse']][names(LOOCV[['kriging_rmse']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    Winning_kriging_model[n, 'Pearson.Cor.model.no.elevation'] <-max(unlist(lapply(LOOCV[['kriging_cor']][names(LOOCV[['kriging_cor']]) %in% c("productsum_Vgm","separable_Vgm") == FALSE]  , function(x)mean(x,na.rm=T))),na.rm=T)
    
    setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/PRCP_FINAL")
    
    if(Winning_kriging_model[n,'Distance.near.station']<2&length(which(is.na(vector_nearest_station[[1]][,c('PRCP')])))==0){
      Winning_kriging_model[n, 'Winning_kriging_model'] <- 'nearest_station'
      Winning_kriging_model[n, 'Total.winning.model(mm)']<-Winning_kriging_model[n, 'Total.near.station(mm)']
      write.table(vector_nearest_station[[1]][,c('Year_Exp', 'PRCP', 'date')],file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    }
    else if (Winning_kriging_model[n, 'RMSE.model.no.elevation']<Winning_kriging_model[n, 'RMSE.model.with.elevation']){
      if(max(no_elev_model[,LOOCV_min_RMSE_noelev])<230&var(no_elev_model[,LOOCV_min_RMSE_noelev])>10){
        Winning_kriging_model[n, 'Winning_kriging_model'] <- paste('no_elevation_',LOOCV_min_RMSE_noelev,sep='')
        Winning_kriging_model[n, 'Total.winning.model(mm)']<-Winning_kriging_model[n, 'Total.no.elevation(mm)']
        selected_data<-no_elev_model[,c('predictions.table',LOOCV_min_RMSE_noelev,'as.character(dates)')]
        colnames(selected_data)<-c('Year_Exp','Kriging.data','dates')
        write.table(selected_data,file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }else {Winning_kriging_model[n, 'Winning_kriging_model'] <-'issue_with_model'}
    }else {
      if(max(elev_model[,LOOCV_min_RMSE_elev])<230&var(elev_model[,LOOCV_min_RMSE_elev])>10){
        Winning_kriging_model[n, 'Winning_kriging_model'] <- paste('elevation_',LOOCV_min_RMSE_elev,sep='')
        Winning_kriging_model[n, 'Total.winning.model(mm)']<-Winning_kriging_model[n, 'Total.with.elevation(mm)']
        selected_data<-elev_model[,c('predictions.table',LOOCV_min_RMSE_elev,'as.character(dates)')]
        colnames(selected_data)<-c('Year_Exp','Kriging.data','dates')
        write.table(selected_data,file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      } else {Winning_kriging_model[n, 'Winning_kriging_model'] <-'issue_with_model'}}
    
    
    if (Winning_kriging_model[n, 'Winning_kriging_model']=='issue_with_model'){
      print(paste0(s,'issue model'))
      if('summetric_Vgm3'%in%colnames(no_elev_model)){
        Winning_kriging_model[n, 'Winning_kriging_model'] <- paste('elevation_','summetric_Vgm3',sep='')
        Winning_kriging_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm3'],na.rm=T)
        selected_data<-no_elev_model[,c('predictions.table','summetric_Vgm3','as.character(dates)')]
        colnames(selected_data)<-c('Year_Exp','Kriging.data','dates')
        
        write.table(selected_data,file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }
      else {
        Winning_kriging_model[n, 'Winning_kriging_model'] <- paste('elevation_','summetric_Vgm5',sep='')
        Winning_kriging_model[n, 'Total.winning.model(mm)']<-sum(no_elev_model[,'summetric_Vgm5'],na.rm=T)
        selected_data<-no_elev_model[,c('predictions.table','summetric_Vgm5','as.character(dates)')]
        colnames(selected_data)<-c('Year_Exp','Kriging.data','dates')
        
        write.table(selected_data,file=paste('selected_model_',s,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      }
      
    }
    
    year=str_sub(s,1,4)
    rm(selected_data)
    setwd(paste0("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/data_files/",year))
    if (s%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH1','2018_ONH2')){
      dat=read.table(paste0(s,'_weather.txt'),header=T,sep='\t')
      Winning_kriging_model[n, 'nb_stations_used_for_kriging'] <-length(unique(dat[which(!is.na(dat$prcp)),'station_name']))}
    else{
    Winning_kriging_model[n, 'nb_stations_used_for_kriging']<-length(unique(read.table(paste0('PRCP_',s,'.txt'),header=T,sep='\t')$id))}
    
    n <- n + 1
  }
  
  
  
}



Winning_kriging_model$Difference.selected.model.field.station=Winning_kriging_model$`Total.winning.model(mm)`-Winning_kriging_model$`Total_field_after_processing(mm)`
Winning_kriging_model$Difference.selected.model.nearest.station=Winning_kriging_model$`Total.winning.model(mm)`-Winning_kriging_model$`Total.near.station(mm)`

##Total or partial replacement by interpolated/nearest stations based on 1)% NA data from field 2)Difference between field and interpolated or nearest stations

Winning_kriging_model[which(Winning_kriging_model$`%NA.field.season`>0.4),'Total_partial_replacement']<-'Total'
Winning_kriging_model[which(Winning_kriging_model$`%NA.field.season`<0.4),'Total_partial_replacement']<-'Partial'

Winning_kriging_model[which(abs(Winning_kriging_model$Difference.selected.model.field.station)>100),'Total_partial_replacement']<-'Total'


setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation")
write.table(Winning_kriging_model,paste0('PRCP_summary.txt'),col.names = T,row.names = F,sep = '\t',quote = F)


# ------------------------------------------------------------------------------
# Meteorological variables: 'WDSP'
# ------------------------------------------------------------------------------

for (meteo_variable in c('WDSP')) {
  
  # 1) Reading the imputed data using IDW and constructing summarized dataset from IDW results before imputation
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/AWND_GSOD_WeatherCAN")
  
  summary = as.data.frame(matrix(NA, ncol = 13, nrow = length(set_locations)))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'growing_season_length_(days)',
    'pearson.cor.daily.field.vs.interpolated',
    'weather.network.used',
    'interpolation_method',
    'distance_1nearest_station_km',
    'name_1nearest_station',
    'pearson.cor.interpolated.vs.nearby.stations',
    'near_or_interpolated_selected'
    
  )
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  predictions <- list()
  
  n = 1
  for (i in set_locations) {
    print(i)
    v=paste0(i,'_IDW.RDS')
    
    if (file.exists(v)){
      
      summary[n, 'interpolation_method'] <- 'inverse_distance_weighting'
      
      if (i%in%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')) {
        summary[n, 'weather.network.used'] <- 'WeatherCAN'
        summary[n, "Year_Exp"] <- i
        summary[n, "flagged.NA.percent.growing.season.field.data"] <-
          sum(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])) /
          length(daily_weather[daily_weather$Year_Exp == i,'MEANWINDSPEED'])
        
        summary[n, "growing_season_length_(days)"] <-
          length(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])
        
        if (length(which(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])))>summary[n, "growing_season_length_(days)"]/2){
          summary[n, "pearson.cor.daily.field.vs.interpolated"] <- NA
        }
        
        else{summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
          cor(as.numeric(as.vector(readRDS(v)[[1]][, 'transformed_u2'])),
              daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'],
              method = 'pearson',
              use = 'complete.obs')}
        
        summary[n,'nb.weather.stations.used'] <- unique(readRDS(v)[[1]][,'nb_stations_used'])
        
        summary[n,'distance_1nearest_station_km'] <- readRDS(v)[['min_dist']]
        
        summary[n,'name_1nearest_station'] <- readRDS(v)[['min_dist_name']]
        
        summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(v)[[1]][, 'transformed_u2'])),
                                                                        as.numeric(as.vector(readRDS(v)[[1]][, 'nearest_station_values'])),
                                                                        method = 'pearson',
                                                                        use = 'complete.obs')
        if(readRDS(v)[['min_dist']]<2&length(which(is.na(readRDS(v)[[1]][,c('nearest_station_values')])))==0){
          print(v)
          df<-readRDS(v)[[1]][,c('dates', 'nearest_station_values')]
          df<-cbind(i,df)
          summary[n,'near_or_interpolated_selected'] <- 'near'
          write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
          
        }else{
          df<-readRDS(v)[[1]][,c('dates', 'transformed_u2')]
          df<-cbind(i,df)
          summary[n,'near_or_interpolated_selected'] <- 'interpolated'
          write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
        } 
        
        
      }  else{
        summary[n, "weather.network.used"] <- 'GSOD'
        summary[n, "Year_Exp"] <- i
        summary[n, "flagged.NA.percent.growing.season.field.data"] <-
          sum(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])) /
          length(daily_weather[daily_weather$Year_Exp == i,'MEANWINDSPEED'])
        
        summary[n, "growing_season_length_(days)"] <-
          length(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])
        
        
        if (length(which(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])))>summary[n, "growing_season_length_(days)"]/2){
          summary[n, "pearson.cor.daily.field.vs.interpolated"] <- NA
        }
        
        else{summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
          cor(as.numeric(as.vector(readRDS(v)[[1]][, 'transformed_u2'])),
              daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'],
              method = 'pearson',
              use = 'complete.obs')}
        
        summary[n,'nb.weather.stations.used'] <- unique(readRDS(v)[[1]][,'nb_stations_used'])
        
        summary[n,'distance_1nearest_station_km'] <- readRDS(v)[[5]]
        
        summary[n,'name_1nearest_station'] <- readRDS(v)[[4]]
        
        summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(v)[[1]][, 'transformed_u2'])),
                                                                        as.numeric(as.vector(readRDS(v)[[1]][, 'nearest_station_values'])),
                                                                        method = 'pearson',
                                                                        use = 'complete.obs')
        
        if(readRDS(v)[['min_dist']]<2&length(which(is.na(readRDS(v)[[1]][,c('nearest_station_values')])))==0){
          print(v)
          df<-readRDS(v)[[1]][,c('dates', 'nearest_station_values')]
          df<-cbind(i,df)
          colnames(df)=c('Year_Exp','dates','transformed_u2')
          summary[n,'near_or_interpolated_selected'] <- 'near'
          write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
          
        }else{
          df<-readRDS(v)[[1]][,c('dates', 'transformed_u2')]
          df<-cbind(i,df)
          summary[n,'near_or_interpolated_selected'] <- 'interpolated'
          write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
        }
      }
    } else {
      
      
      print(paste0(i,': only 1 station'))
      summary[n,'near_or_interpolated_selected'] <- 'near'
      summary[n, 'interpolation_method'] <- 'only_1_station_available'
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- i
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])) /
        length(daily_weather[daily_weather$Year_Exp == i,'MEANWINDSPEED'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'])
      
      if (!all(is.na(daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED']))){
        summary[n, "pearson.cor.daily.field.vs.interpolated"] <-
          cor(as.numeric(as.vector(readRDS(paste0(i,'_1station.RDS'))[[1]][, 'transformed_u2'])),
              daily_weather[daily_weather$Year_Exp == i, 'MEANWINDSPEED'],
              method = 'pearson',
              use = 'complete.obs')}
      else{summary[n, "pearson.cor.daily.field.vs.interpolated"] <-NA}
      
      summary[n,'nb.weather.stations.used'] <- 1
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(paste0(i,'_1station.RDS'))[[4]]
      
      
      summary[n,'name_1nearest_station'] <- readRDS(paste0(i,'_1station.RDS'))[[5]]
      
      
      summary[n,'pearson.cor.interpolated.vs.nearby.stations'] <- NA
      
      df<-readRDS(paste0(i,'_1station.RDS'))[[1]][,c('date', 'transformed_u2')]
      df<-cbind(i,df)
      
      
      write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      
      
      
    }
    
    
    
    n <- n + 1
    
    
  }
  
  ## Wind replaced by imputed values --> lots of missing data + measurements estimated at 10 and transformed to 2 meters, whereas weather stations measurements are probably at a lower height
  ## Values always lower from field stations compared to interpolated data, although the correlation looks good --> as said before, relative to wind height measurement.
  summary$Replacement<-'Total'
  
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



# ------------------------------------------------------------------------------
# Meteorological variables: 'HMEAN','HMAX','HMIN'
# ------------------------------------------------------------------------------


for (meteo_variable in c('HMEAN')) {
  
  # 1) Reading the imputed data using IDW and constructing summarized dataset from IDW results before imputation
  
  setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN")
  
  summary = as.data.frame(matrix(NA, ncol = 13, nrow = length(set_locations)-9))
  colnames(summary) = c(
    'Year_Exp',
    'meteo_variable',
    'flagged.NA.percent.growing.season.field.data',
    'nb.weather.stations.used',
    'radius',
    'weather.network.used',
    'growing_season_length_(days)',
    'RHMEAN.pearson.cor.daily.field.vs.interpolated',
    'interpolation_method',
    'distance_1nearest_station_km',
    'name_1nearest_station',
    'RHMEAN.pearson.cor.interpolated.vs.nearby.stations',
    'near_or_interpolated_selected'
  )
  
  
  summary$meteo_variable = eval(meteo_variable)
  summary$radius = 70
  
  predictions <- list()
  
  n = 1
  for (i in set_locations[set_locations%notin%c('2014_ONH1','2014_ONH2','2015_ONH1','2015_ONH2','2016_ONH1','2016_ONH2','2017_ONH1','2017_ONH2','2018_ONH2')]) {
    print(i)
    v=paste0(i,'_IDW.RDS')
    
    if (file.exists(v)){
      
      summary[n, 'interpolation_method'] <- 'inverse_distance_weighting'
      
      
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- i
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == i, 'HMEAN'])) /
        length(daily_weather[daily_weather$Year_Exp == i,'HMEAN'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == i, 'HMEAN'])
      
      
      if (length(which(is.na(daily_weather[daily_weather$Year_Exp == i, 'HMEAN'])))>summary[n, "growing_season_length_(days)"]/2){
        summary[n, "pearson.cor.daily.field.vs.interpolated"] <- NA
      }
      
      else{summary[n, "RHMEAN.pearson.cor.daily.field.vs.interpolated"] <-
        cor(as.numeric(as.vector(readRDS(v)[[1]][, 'var1.pred'])),
            daily_weather[daily_weather$Year_Exp == i, 'HMEAN'],
            method = 'pearson',
            use = 'complete.obs')}
      
      summary[n,'nb.weather.stations.used'] <- unique(readRDS(v)[[1]][,'nb_stations_used'])
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(v)[[5]]
      
      summary[n,'name_1nearest_station'] <- readRDS(v)[[4]]
      
      summary[n,'RHMEAN.pearson.cor.interpolated.vs.nearby.stations'] <- cor(as.numeric(as.vector(readRDS(v)[[1]][, 'var1.pred'])),
                                                                             as.numeric(as.vector(readRDS(v)[[1]][, 'nearest_station_values'])),
                                                                             method = 'pearson',
                                                                             use = 'complete.obs')
      
      if(readRDS(v)[['min_dist']]<2&length(which(is.na(readRDS(v)[[1]][,c('nearest_station_values')])))==0){
        print(v)
        df<-readRDS(v)[[1]][,c('dates', 'nearest_station_values')]
        df<-cbind(i,df)
        colnames(df)=c('Year_Exp','dates','HMEAN')
        summary[n,'near_or_interpolated_selected'] <- 'near'
        write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
        
      }else{
        df<-readRDS(v)[[1]][,c('dates', 'var1.pred')]
        df<-cbind(i,df)
        summary[n,'near_or_interpolated_selected'] <- 'interpolated'
        write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
        
      }
    } else {
      
      
      print(paste0(i,': only 1 station'))
      summary[n,'near_or_interpolated_selected'] <- 'near'
      summary[n, 'interpolation_method'] <- 'only_1_station_available'
      summary[n, "weather.network.used"] <- 'GSOD'
      summary[n, "Year_Exp"] <- i
      summary[n, "flagged.NA.percent.growing.season.field.data"] <-
        sum(is.na(daily_weather[daily_weather$Year_Exp == i, 'HMEAN'])) /
        length(daily_weather[daily_weather$Year_Exp == i,'HMEAN'])
      
      summary[n, "growing_season_length_(days)"] <-
        length(daily_weather[daily_weather$Year_Exp == i, 'HMEAN'])
      
      if (!all(is.na(daily_weather[daily_weather$Year_Exp == i, 'HMEAN']))){
        summary[n, "RHMEAN.pearson.cor.daily.field.vs.interpolated"] <-
          cor(as.numeric(as.vector(readRDS(paste0(i,'_1station.RDS'))[[1]][, 'HMEAN'])),
              daily_weather[daily_weather$Year_Exp == i, 'HMEAN'],
              method = 'pearson',
              use = 'complete.obs')}
      else{summary[n, "RHMEAN.pearson.cor.daily.field.vs.interpolated"] <-NA}
      
      summary[n,'nb.weather.stations.used'] <- 1
      
      summary[n,'distance_1nearest_station_km'] <- readRDS(paste0(i,'_1station.RDS'))[[4]]
      
      
      summary[n,'name_1nearest_station'] <- readRDS(paste0(i,'_1station.RDS'))[[5]]
      
      
      summary[n,'RHMEAN.pearson.cor.interpolated.vs.nearby.stations'] <- NA
      
      df<-readRDS(paste0(i,'_1station.RDS'))[[1]][,c('date', 'HMEAN')]
      df<-cbind(i,df)
      
      
      
      write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
      
      
      
    }
    
    
    
    n <- n + 1
    
    
  }
  
  summary$Replacement<-NA
  summary[which(summary$RHMEAN.pearson.cor.daily.field.vs.interpolated<0.6),'Replacement']<-'Total'
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



setwd(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/HMEAN_GSOD_WeatherCAN"
)

summary = as.data.frame(matrix(NA, ncol = 15, nrow = 9))
colnames(summary) = c(
  'Year_Exp',
  'meteo_variable',
  'flagged.NA.percent.growing.season.field.data',
  'nb.weather.stations.used',
  'radius',
  'weather.network.used',
  'growing_season_length_(days)',
  'RHMAX.pearson.cor.daily.field.vs.interpolated',
  'RHMIN.pearson.cor.daily.field.vs.interpolated',
  'interpolation_method',
  'distance_1nearest_station_km',
  'name_1nearest_station',
  'RHMAX.pearson.cor.interpolated.vs.nearby.stations',
  'RHMIN.pearson.cor.interpolated.vs.nearby.stations',
  'near_or_interpolated_selected'
)


summary$radius = 70
summary$meteo_variable = 'HMIN;HMAX'




n = 1
for (i in set_locations[set_locations %in% c(
  '2014_ONH1',
  '2014_ONH2',
  '2015_ONH1',
  '2015_ONH2',
  '2016_ONH1',
  '2016_ONH2',
  '2017_ONH1',
  '2017_ONH2',
  '2018_ONH2'
)]) {
  print(i)
  v=paste0(i,'_IDW.RDS')
  
  
  summary[n, 'interpolation_method'] <-
    'inverse_distance_weighting'
  
  
  summary[n, 'weather.network.used'] <- 'WeatherCAN'
  summary[n, "Year_Exp"] <- i
  
  summary[n, "flagged.NA.percent.growing.season.field.data"] <-
    sum(is.na(daily_weather[daily_weather$Year_Exp == i, c('HMAX')])) /
   length(daily_weather[daily_weather$Year_Exp == i, c('HMAX')])
  summary[n, "growing_season_length_(days)"] <-
    length(daily_weather[daily_weather$Year_Exp == i, c('HMAX')])
  
  
  
  if (i %in% c('2017_ONH1', '2017_ONH2')) {
    summary[n,  'RHMAX.pearson.cor.daily.field.vs.interpolated'] <-
      NA
    summary[n,  'RHMIN.pearson.cor.daily.field.vs.interpolated'] <-
      NA
  }
  
  else{
    
    summary[n,  'RHMAX.pearson.cor.daily.field.vs.interpolated'] <-
      cor(as.numeric(as.vector(readRDS(v)[[1]][, 'rhmax'])),
          daily_weather[daily_weather$Year_Exp == i, 'HMAX'],
          method = 'pearson',
          use = 'complete.obs')
    summary[n,  'RHMIN.pearson.cor.daily.field.vs.interpolated'] <-
      cor(as.numeric(as.vector(readRDS(v)[[1]][, 'rhmin'])),
          daily_weather[daily_weather$Year_Exp == i, 'HMIN'],
          method = 'pearson',
          use = 'complete.obs')
  }
  
  summary[n,  'RHMAX.pearson.cor.interpolated.vs.nearby.stations'] <-
    cor(
      as.numeric(as.vector(readRDS(v)[[1]][, 'rhmax'])),
      as.numeric(as.vector(readRDS(v)[[1]][, 8][, 1])),
      method = 'pearson',
      use = 'complete.obs'
    )
  summary[n,  'RHMIN.pearson.cor.interpolated.vs.nearby.stations'] <-
    cor(
      as.numeric(as.vector(readRDS(v)[[1]][, 'rhmin'])),
      as.numeric(as.vector(readRDS(v)[[1]][, 8][, 2])),
      method = 'pearson',
      use = 'complete.obs'
    )
  summary[n, 'nb.weather.stations.used'] <-
    unique(readRDS(v)[[1]][, 'nb_stations_used'])
  
  summary[n, 'distance_1nearest_station_km'] <-
    readRDS(v)[['min_dist']]
  
  summary[n, 'name_1nearest_station'] <-
    readRDS(v)[['min_dist_name']]
  
  
  if(readRDS(v)[['min_dist']]<2&length(which(is.na(readRDS(v)[[1]][,8])))==0){
    print(v)
    df<-readRDS(v)[[1]][,c(3,8,9)]
    df<-cbind(i,df)
    summary[n,'near_or_interpolated_selected'] <- 'near'
    write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
    
  }else{
    df<-readRDS(v)[[1]][,c('dates', 'rhmin','rhmax')]
    df<-cbind(i,df)
    summary[n,'near_or_interpolated_selected'] <- 'interpolated'
    write.table(df,file=paste('selected_model_',i,'.txt',sep=''),col.names = T,row.names = F,sep = '\t',quote=F)
  }
  
  
  
  n <- n + 1
  
  
  
  
  
}
summary$Replacement<-NA
summary[which(summary$RHMAX.pearson.cor.daily.field.vs.interpolated<0.6),'Replacement']<-'Total'
summary[which(summary$flagged.NA.percent.growing.season.field.data>0.5),'Replacement']<-'Total'
summary[summary$Replacement%notin%'Total','Replacement']<-'Partial'
summary[summary$Year_Exp%in%c('2017_ONH1', '2017_ONH2'),'Replacement']<-'Total'



write.table(
  summary,
  file = paste(
    "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/INTERPOLATED_DATA/imputation/",
    'HMIN_HMAX',
    '_summary.txt',
    sep = ''
  ),
  col.names = T,
  row.names = F,
  quote = F,
  sep = '\t'
)
