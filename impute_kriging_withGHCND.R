#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGHCND} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable_GHCND Character. GHCND element names (= variable measured at the GHCND station): TMAX, TMIN, PRCP
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested' and @variable_to_impute
#' @param variable_to_impute. Character.
#' @param name_in_table. Character


impute_kriging_withGHCND <- function(Year_Exp,radius=50,meteo_variable_GHCND,daily_weather=daily_weather,variable_to_impute,name_in_table) {
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
  
  d=read.table(file=paste('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GHCND/ghcnd_files/',year,'/',variable_to_impute,'_',Year_Exp,'.txt',sep = ''),header=T)
  
  colnames(d)[3]='dates'
  ind=which(colnames(d)%in%c('Year_Exp','variable'))
  d=d[,-ind]
  
  d$date=as.Date(as.character(d$date))
  
  ind2 <- which(!is.na(d[,eval(meteo_variable_GHCND)]))
  d <- d[ind2,]
  
  print('Data from GHCND stations prepared')
  
  
  ########################
  ########################
  
  
  if(length(unique(d$id))<3){
    d$day<-lubridate::yday(d$dates)
    predicted.values<-d[which(d$day==yday(date_start)):which(d$day==yday(date_end)),meteo_variable_GHCND]
    
    
    if (!is.null(variable_to_impute)) {
      cors=cor(predicted.values,field_values,use = 'complete.obs')}
    
    
    dates=seq(as.Date(date_start,tz="CET"),as.Date(date_end,tz="CET"),by='days')
    predictions.table=cbind(Year_Exp,predicted.values,as.character(dates))
    colnames(predictions.table)=c('Year_Exp',paste(variable_to_impute),'dates')
    
    to_save=list(predictions.table,cors)
    names(to_save)<-c('predictions_YearExp','cors_YearExp')
    
    
    }
  ########################
  ####Ordinary kriging####
  if(length(unique(d$id))>=3){
  sub=d[,which(colnames(d)%in%c('station','longitude','latitude',meteo_variable_GHCND,'dates'))]
  
  sp::coordinates(sub)=c('longitude','latitude')
  proj4string(sub) = "+proj=longlat +datum=WGS84"
  #projection(sub)=CRS("+init=epsg:4326")
  
  #Transform into Mercator Projection
  tempmin.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 
  
  
  tempminSP <- SpatialPoints(tempmin.UTM@coords,CRS("+init=epsg:3395"))
  
  ######
  
  #tempminDF <- data.frame(values=tempmin.UTM$HMEAN )
  tempminDF <- data.frame(d[,meteo_variable_GHCND])
  colnames(tempminDF)<-'values'
  ####
  
  
  tempminTM <- as.POSIXct(date(tempmin.UTM$dates))
  #combine the 3 objects
  timeDF <- STIDF(tempminSP,tempminTM,data=tempminDF) 
  
 
  
  #variogram
  print('Computation variogram starts:')
  var <- variogramST(values~1,data=timeDF,tunit="days",assumeRegular=F,na.omit=T) 
  pdf(paste('GHCND/imputation/',variable_to_impute,'/variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
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
  pdf(paste('GHCND/imputation/',variable_to_impute,'/fitted_variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
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
  }
  
  saveRDS(to_save,file=paste('GHCND/imputation/',variable_to_impute,'/',Year_Exp,'.RDS',sep=''))
  
  #return(to_save)
  
  
  
  
  
  
}





