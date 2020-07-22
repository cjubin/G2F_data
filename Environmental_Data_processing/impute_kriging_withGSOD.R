#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable_GSOD Character. GSOD element names (= variable measured at the ISD STATION): TEMP, MAX, MIN, DEWP, STP, WDSP, PRCP
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested' and @variable_to_impute
#' @param variable_to_impute. Character.
#' @param name_in_table. Character


impute_kriging_withGSOD <- function(Year_Exp,radius=70,meteo_variable_GSOD,daily_weather=daily_weather,variable_to_impute,name_in_table) {
  source('fahrenheit_to_celsius.R')
  
  #options(noaakey = "ueWgGjcckAdRLEXbpNtePVgbRWXmiQBG")
  
  print(Year_Exp)
  #Retrieve information about the experiment
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  date_start=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Planted'])),origin=paste(year-1,'12','31',sep = '-'))
  date_end=as.Date(as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Date.Harvested'])),origin=paste(year-1,'12','31',sep = '-'))
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  #Values recorded for this variable at the field station
  if (!is.null(name_in_table)&name_in_table%in%colnames(daily_weather)) {
    field_values = daily_weather[daily_weather$Year_Exp == Year_Exp, name_in_table]
  }
  
  #Finding the closest stations in a certain radius and select those for which coverage period 2013-2019 is sure
  stations_close=as.data.frame(rnoaa::isd_stations_search(lat =latitude, lon = longitude, radius = radius))
  stations_close$idstation=paste(stations_close$usaf,stations_close$wban,sep='')
  stations_close<-arrange(stations_close,distance)
  stations_close$yearstart<-substr(stations_close$begin,0,nchar(stations_close$begin)-4)
  stations_close$yearend<-substr(stations_close$end,0,nchar(stations_close$end)-4)
  stations_close<-filter(stations_close,yearstart<2013&yearend>2019)
  
  #files_stations: path to find weather files for the stations used in kriging (radius 80 km) (if data for theses ISD stations is available)
  
  files_stations <-
    paste0('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GSOD/gsod_files',
           "/",
	   year,
	   '/',
           stations_close$idstation,
           ".csv")
  #List of the csv files in the temporary folder
  GSOD_list <-
    list.files(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GSOD/gsod_files',
      pattern = "*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
  
  #Overlap between all initial selected stations and those for which data could be eventually downloaded on the computer.
  GSOD_list2 <-
    subset(GSOD_list, GSOD_list %in% files_stations)
  
  data_frames=lapply(GSOD_list2,function(x)read.csv(x))
  all_data<-plyr::compact(data_frames)
  all_data=as.data.frame(do.call('rbind',all_data))
  
  d = cbind(
    'station' = all_data$STATION,
    'longitude' = all_data$LONGITUDE,
    'latitude' = all_data$LATITUDE,
     all_data[, colnames(all_data) %in% meteo_variable_GSOD],
    'dates' = as.character(as.vector(all_data$DATE))
  )
  d=as.data.frame(d) 
  
  if (length(meteo_variable_GSOD)==1){colnames(d)[4]=meteo_variable_GSOD}
  
  d$longitude=as.numeric(as.vector(d$longitude))  
  d$latitude=as.numeric(as.vector(d$latitude)) 
  
  if (any(meteo_variable_GSOD%in%c('DEWP'))) {
    index=which(d$DEWP>999)
    if(length(index)>0){d<-d[-index,]}
    
  }
  
  if (any(meteo_variable_GSOD%in%c('TEMP', 'DEWP', 'MAX', 'MIN'))) {
    if (length(meteo_variable_GSOD)>1){ d[,eval(meteo_variable_GSOD)] <- fahrenheit_to_celsius(d[,eval(meteo_variable_GSOD)])}
    if (length(meteo_variable_GSOD)==1){ d[,eval(meteo_variable_GSOD)] <- fahrenheit_to_celsius(as.numeric(as.vector(d[,eval(meteo_variable_GSOD)])))}
    
  }
  
  if (any(meteo_variable_GSOD%in%c('DEWP'))){
    d$HMEAN<-NA
    d[,'HMEAN'] <- 100*(exp((17.625*d$DEWP)/(243.04+d$DEWP))/exp((17.625*d$TEMP)/(243.04+d$TEMP)))
  }
  
  d$dates=as.Date(d$dates)
  
  for (i in colnames(d)[colnames(d)%notin%c('station','dates','latitude','longitude')]) {
    d[,i]<-as.numeric(as.vector(d[,i]))
  }
  
  print('Data from GSOD stations prepared')
  
  ########################
  ####Ordinary kriging####
  
  sub=d[,which(colnames(d)%in%c('station','longitude','latitude',variable_to_impute,'dates'))]
  
  sp::coordinates(sub)=c('longitude','latitude')
  proj4string(sub) = "+proj=longlat +datum=WGS84"
  #projection(sub)=CRS("+init=epsg:4326")
  
  #Transform into Mercator Projection
  tempmin.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 
  
  
  tempminSP <- SpatialPoints(tempmin.UTM@coords,CRS("+init=epsg:3395"))
  
  ######
  
  #tempminDF <- data.frame(values=tempmin.UTM$HMEAN )
  tempminDF <- data.frame(d[,variable_to_impute])
  colnames(tempminDF)<-'values'
  ####
  
  
  tempminTM <- as.POSIXct(date(tempmin.UTM$dates))
  #combine the 3 objects
  timeDF <- STIDF(tempminSP,tempminTM,data=tempminDF) 
  
  
  #variogram
  print('Computation variogram starts:')
  var <- variogramST(values~1,data=timeDF,tunit="days",assumeRegular=F,na.omit=T) 
  pdf(paste('GSOD/imputation/',variable_to_impute,'/variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
  print(plot(var,map=F))
  dev.off()

  
  ######Different models under assessment: we need to fit a model to our variogram.
  ######5 variograms models are possible: separable, product sum, metric, sumMetric, simpleSum metric
  
  #0 Separable metric model
  separable <- vgmST("separable", space = vgm(-60,"Sph", 500, 1),time = vgm(35,"Sph", 500, 1), sill=0.56) 
  #Automatic fit
  tryCatch({
    separable_Vgm=fit.StVariogram(var,separable)
    separable_mse<-attr(separable_Vgm, "MSE")}
    , error = function(e)
      NULL
  ) 
 
  #1 product sum
  productsum <-vgmST("productSum",space = vgm(1, "Exp", 150, 0.5),time = vgm(1, "Exp", 5, 0.5),k = 50) 
  #Automatic fit
  tryCatch({
    productsum_Vgm=fit.StVariogram(var,productsum)
    productsum_mse<-attr(productsum, "MSE")}
    , error = function(e)
      NULL
  )
  
  
  #2 metric model
  metric <-vgmST("metric", joint = vgm(50,"Mat", 500, 0), stAni=200) 
  #Automatic fit
  tryCatch({
  metric_Vgm=fit.StVariogram(var,metric)
  metric_mse<-attr(metric_Vgm, "MSE")}
  , error = function(e)
    NULL
  )
  
  #3 Sum metric model
  summetric <-vgmST("sumMetric", space = vgm(psill=5,"Sph", range=500, nugget=0),time = vgm(psill=500,"Sph", range=500, nugget=0), joint = vgm(50,"Mat", range=500, nugget=10), stAni=1) 
  #Automatic fit
  tryCatch({
  summetric_Vgm=fit.StVariogram(var,summetric)
  summetric_mse<-attr(summetric_Vgm, "MSE")}
  , error = function(e)
    NULL
  )

  
  #4 Simple sum metric model
  simplesummetric <-vgmST("simpleSumMetric",space = vgm(5,"Sph", 500, 0),time = vgm(500,"Sph", 500, 0), joint = vgm(1,"Sph", 500, 0), nugget=1, stAni=500) 
  #Automatic fit
  tryCatch({
  simplesummetric_Vgm=fit.StVariogram(var,simplesummetric)
  simplesummetric_mse<-attr(simplesummetric_Vgm, "MSE")}
  , error = function(e)
  NULL
  )
  
  
  
  ##Cross-validate the method 
  library(dismo)
  nfolds<-5
  k<-kfold(timeDF,nfolds)
  
  list_models <- list()
  if (exists(paste(quote(separable_Vgm)))){list_models<-append(list_models,list(separable_Vgm))
  names(list_models)[[length(list_models)]]<-'separable_Vgm'}
  if (exists(paste(quote(productsum_Vgm)))){list_models<-append(list_models,list(productsum_Vgm))
  names(list_models)[[length(list_models)]]<-'productsum_Vgm'}
  if (exists(paste(quote(metric_Vgm)))){list_models<-append(list_models,list(metric_Vgm))
  names(list_models)[[length(list_models)]]<-'metric_Vgm'}
  if (exists(paste(quote(summetric_Vgm)))){list_models<-append(list_models,list(summetric_Vgm))
  names(list_models)[[length(list_models)]]<-'summetric_Vgm'}
  if (exists(paste(quote(simplesummetric_Vgm)))){list_models<-append(list_models,list(simplesummetric_Vgm))
  names(list_models)[[length(list_models)]]<-'simplesummetric_Vgm'}
  
  kriging_rmse<-vector(mode = 'list',length = length(list_models))
  kriging_cor<-vector(mode = 'list',length = length(list_models))
  names(kriging_cor)<-names(list_models)
  names(kriging_rmse)<-names(list_models)
  
  
  for (i in 1:nfolds) {
    train<- timeDF[k != i, ]
    test <- timeDF[k == i, ]
    values=test@data$values
    
    for (model in 1:length(list_models)) {
      tryCatch({
        p1 <-
          krigeST(
            values ~ 1,
            data = train,
            model = list_models[[model]],
            newdata = test
          )
        predictions_test = p1@data$var1.pred
        kriging_cor[[model]][i] <- cor(predictions_test, values)
        kriging_rmse[[model]][i] <-
          sqrt(mean((predictions_test - values) ^ 2))
        
      }, error = function(e)
        NULL)
      
    }
  }
  
  
  
  ##Choose the appropriate model based on cross-validation results (minimum average RMSE)
  fitted.stvgm=get(names(which.min(sapply(kriging_rmse, function(x)mean(x,na.rm=TRUE)))))
  
  
  ##Plot results
  par(mfrow=c(2,1))
  pdf(paste('GSOD/imputation/',variable_to_impute,'/fitted_variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
  print(plot(var,fitted.stvgm,map=F))
  dev.off()

  
  #Prediction grid: growing season for the field experiment described by Year_Exp
  field=vector(mode = 'numeric',length = 2)
  field<-c(longitude,latitude)
  
  field=as.data.frame(matrix(field,ncol = 2))
  colnames(field)<-c('longitude','latitude')
  sp::coordinates(field)=c('longitude','latitude')
  proj4string(field) = "+proj=longlat +datum=WGS84"
  field_grid<- spTransform(field,CRS("+init=epsg:3395")) 
  tm.grid<-seq(as.POSIXct(date_start,tz="CET"),as.POSIXct(date_end,tz="CET"),by='days')
  
  grid.ST <- STF(field_grid,tm.grid) 
  
  
  pred<-krigeST(values~1,data=timeDF,modelList =fitted.stvgm,newdata = grid.ST )
  predicted.values=pred@data$var1.pred
  if (!is.null(variable_to_impute)) {
    cors=cor(predicted.values,field_values,use = 'complete.obs')}
  
  
  dates=seq(as.Date(date_start,tz="CET"),as.Date(date_end,tz="CET"),by='days')
  predictions.table=cbind(Year_Exp,predicted.values,as.character(dates))
  colnames(predictions.table)=c('Year_Exp',paste(variable_to_impute),'dates')

  to_save=list(predictions.table,cors,kriging_cor,kriging_rmse)
  names(to_save)<-c('predictions_YearExp','cors_YearExp','5f.cv.kriging.cor','5f.cv.kriging.rmse')
  saveRDS(to_save,file=paste('GSOD/imputation/',variable_to_impute,'/',Year_Exp,'.RDS',sep=''))

  #return(to_save)
  
  
  
  
  
  
}





