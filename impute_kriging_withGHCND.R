#' Imputation of a meteorological variable based on a field location with geographical coordinates
#'
#' \code{impute_kriging_withGHCND} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param radius Numeric. Distance from the field location to consider to interpolate.
#' @param meteo_variable_GHCND Character. GHCND element names (= variable measured at the GHCND station): TMAX, TMIN, PRCP
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, 'long', 'lat', 'Date.Planted', 'Date.Harvested', and optionally meteo_variable_in_table with its value as column name. Ex: 'TMIN'.
#' @param meteo_variable_in_table. Character with the column name in the table daily_weather of the meteorological variable of interest which has to be imputed from surrounding stations. It does not have to be the same as @meteo_variable_GHCND



impute_kriging_withGHCND <- function(Year_Exp,radius=50,meteo_variable_GHCND,daily_weather=daily_weather,meteo_variable_in_table=NULL) {
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
  if (!is.null(meteo_variable_in_table)) {
    field_values = daily_weather[daily_weather$Year_Exp == Year_Exp, meteo_variable_in_table]
  }
  
  #Download data alreafdy formatted for GHCND files
  
  d=read.table(file=paste('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/GHCND/ghcnd_files/',meteo_variable_GHCND,'_',Year_Exp,'.txt',sep = ''))
  colnames(d)[2]='values'
  colnames(d)[3]='dates'
  d=d[,which(colnames(d)%in%c('Year_Exp','variable'))]
  
  d$date=as.Date(d$date)
  d=arrange(d,dates)
  print('Data from GHCND stations prepared')
  
  ########################
  ####Ordinary kriging####
  
  sub=d
  sp::coordinates(sub)=c('longitude','latitude')
  proj4string(sub) = "+proj=longlat +datum=WGS84"
  #projection(sub)=CRS("+init=epsg:4326")
  
  #Transform into Mercator Projection
  tempmin.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 
  
  
  tempminSP <- SpatialPoints(tempmin.UTM@coords,CRS("+init=epsg:3395"))
  tempminDF <- data.frame(values=tempmin.UTM$values) 
  tempminTM <- as.POSIXct(date(tempmin.UTM$dates))
  #combine the 3 objects
  timeDF <- STIDF(tempminSP,tempminTM,data=tempminDF) 
  
  
  #variogram
  print('Computation variogram starts:')
  var <- variogramST(values~1,data=timeDF,tunit="days",assumeRegular=F,na.omit=T) 
  pdf(paste('GHCND/imputation/',meteo_variable_in_table,'/variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
  print(plot(var,map=F))
  dev.off()
  
  
  ######Different models under assessment: we need to fit a model to our variogram.
  ######5 variograms models are possible: separable, product sum, metric, sumMetric, simpleSum metric
  
  #1 Separable
  separable <-vgmST("separable", space = vgm(-60,"Sph", 500, 1),time = vgm(35,"Sph", 500, 1), sill=0.56) 
  #Automatic fit
  separable_Vgm=fit.StVariogram(var,separable)
  separable_mse<-attr(separable_Vgm, "MSE")
  
  #2 product sum
  productsum <-vgmST("productSum",space = vgm(1, "Exp", 150, 0.5),time = vgm(1, "Exp", 5, 0.5),k = 50) 
  #Automatic fit
  productsum_Vgm=fit.StVariogram(var,productsum)
  productsum_mse<-attr(productsum, "MSE")
  
  #3 metric model
  metric <-vgmST("metric", joint = vgm(50,"Mat", 500, 0), stAni=200) 
  #Automatic fit
  metric_Vgm=fit.StVariogram(var,metric)
  metric_mse<-attr(metric_Vgm, "MSE")
  
  #4 Sum metric model
  summetric <-vgmST("sumMetric", space = vgm(psill=5,"Sph", range=500, nugget=0),time = vgm(psill=500,"Sph", range=500, nugget=0), joint = vgm(50,"Mat", range=500, nugget=10), stAni=1) 
  #Automatic fit
  summetric_Vgm=fit.StVariogram(var,summetric)
  summetric_mse<-attr(summetric_Vgm, "MSE")
  
  
  #5 Simple sum metric model
  simplesummetric <-vgmST("simpleSumMetric",space = vgm(5,"Sph", 500, 0),time = vgm(500,"Sph", 500, 0), joint = vgm(1,"Sph", 500, 0), nugget=1, stAni=500) 
  #Automatic fit
  simplesummetric_Vgm=fit.StVariogram(var,simplesummetric)
  simplesummetric_mse<-attr(simplesummetric_Vgm, "MSE")
  
  
  
  ##Cross-validate the method 
  library(dismo)
  nfolds<-5
  k<-kfold(timeDF,nfolds)
  
  list_models <- list(separable_Vgm,productsum_Vgm,metric_Vgm,summetric_Vgm,simplesummetric_Vgm)
  
  kriging_rmse<-vector(mode = 'list',length = 5)
  kriging_cor<-vector(mode = 'list',length = 5)
  
  
  for (i in 1:nfolds) {
    train<- timeDF[k != i, ]
    test <- timeDF[k == i, ]
    
    
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
        kriging_cor[[model]][i] <- cor(predictions_test, test$values)
        kriging_rmse[[model]][i] <-
          sqrt(mean((predictions_test - test$values) ^ 2))
        
      }, error = function(e)
        NULL)
      
    }
  }
  
  names(kriging_cor)<-c('separable_Vgm','productsum_Vgm','metric_Vgm','summetric_Vgm','simplesummetric_Vgm')
  names(kriging_rmse)<-c('separable_Vgm','productsum_Vgm','metric_Vgm','summetric_Vgm','simplesummetric_Vgm')
  
  
  ##Choose the appropriate model based on cross-validation results (minimum average RMSE)
  fitted.stvgm=get(names(which.min(sapply(kriging_rmse, function(x)mean(x,na.rm=TRUE)))))
  
  
  ##Plot results
  par(mfrow=c(2,1))
  pdf(paste('GHCND/imputation/',meteo_variable_in_table,'/fitted_variogram',Year_Exp,'.pdf',sep='') ,width = 8,height = 8)
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
  if (!is.null(meteo_variable_in_table)) {
    cors=cor(predicted.values,field_values,use = 'complete.obs')}
  
  
  dates=seq(as.Date(date_start,tz="CET"),as.Date(date_end,tz="CET"),by='days')
  predictions.table=cbind(Year_Exp,predicted.values,as.character(dates))
  colnames(predictions.table)=c('Year_Exp',paste(meteo_variable_GHCND),'dates')
  
  to_save=list(predictions.table,cors,kriging_cor,kriging_rmse)
  names(to_save)<-c('predictions_YearExp','cors_YearExp','5f.cv.kriging.cor','5f.cv.kriging.rmse')
  saveRDS(to_save,file=paste('GHCND/imputation/',meteo_variable_in_table,'/',Year_Exp,'.RDS',sep=''))
  
  return(to_save)
  
  
  
  
  
  
}





