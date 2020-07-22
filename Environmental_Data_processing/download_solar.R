#' Download from NCEI the GSOD data files corresponding to the Year_Exp for the complete year (for further modeling of temporal correlation) and for stations present nearby the Year_Exp.
#' Recommend to download locally files (not on server)
#' \code{impute_kriging_withGSOD} interpolates values for a specific meteorological variable from ISD stations given a time frame for a specific location
#' @param Year_Exp Character. Experiment (associated iwth a specific field location) in the G2F dataset which needs to be imputed.
#' @param daily_weather. Data.frame containing at least the following columns: a column 'Year_Exp' containing the specific element used in @Year_Exp, a column 'long', and column 'lat'.



download_solar <- function(Year_Exp,daily_weather=daily_weather) {
  library(nasapower)
  library(data.table)
  source('safeguarding.R')
  
  
  # ------------------------------------------------------------------------------
  #Retrieve information about the experiment
  # ------------------------------------------------------------------------------
  year=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'Year']))])
  longitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'long']))])
  latitude=as.numeric(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat'])[!is.na(unique(daily_weather[daily_weather$Year_Exp==Year_Exp,'lat']))])
  
  
  dat=nasapower::get_power(community = 'AG',lonlat = c(longitude,latitude),dates = c(paste(year,'-01-01',sep = ''),paste(year,'-12-31',sep = '')),pars='ALLSKY_SFC_SW_DWN',temporal_average = 'DAILY')
 
  
  dat2=cbind(dat$ALLSKY_SFC_SW_DWN,as.character(dat$YYYYMMDD),dat$LON,dat$LAT,Year_Exp,variable='ALLSKY_SFC_SW_DWN')
  colnames(dat2)<-c('ALLSKY_SFC_SW_DWN','date','LON','LAT','Year_Exp','variable')
  write.table(dat2,file=paste('C:/Users/cathyjubin/Documents/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/NASA/',Year_Exp,'.txt',sep = ''),col.names = T,row.names = F,sep = '\t',quote = F)
  
  return(all_data)
  
 
  
}




