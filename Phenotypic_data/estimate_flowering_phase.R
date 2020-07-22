#' Find approximate dates for flowering phase, i.e. from
#' initiation of silks to the date for which abortion does not occur anymore,
#' based on the weather data (only thermal time)
#'
#' \code{estimate_flowering_phase} returns a data.frame with approximate start and end of flowering time for each hybrid (i.e. pedigree column),
#' based on data from other experiments with these hybrids and using weather data at the location to predict.
#'
#' @param data_to_predict data.frame with 3 or 4 columns:
#' - 3 columns: both begin and end of flowering phase to predict. First column 'Pedigree', second column the 'Year_Exp' ID, third column the 'Planting.Date '(numeric, in Day.of.Year) 
#' - 4 columns: only end of flowering phase to predict (silking was scored). Fourth column for BLUEs of silking dates: 'Date.Silking'
#' @param data_estimation_model data.frame with 6 columns containing first column 'Pedigree', second column the 'Year_Exp' ID, 
#' third column the 'Planting.Date' (numeric, in Day.of.Year), 
#' fourth column the 'Silking.Date' (as Day.of.Year, numeric), 
#' fifth column the 'GDU' reached at Silk DAP, 
#' and sixth column 'GS.duration',the duration of the growing season. 
#' 
#' @param GDD_place_to_predict data.frame with colnames 'Date' and 'GDU'. 
#' The data.frame containsin first column all dates from Planting to Harvest dates corresponding to the place to predict.
#' Second column contains the cooresponding cumulative GDUs (real weather data from the concerned place).
#' Estimation of the end of the flowering phase --> 280 GDUs after beginning of silking (R1 + R2 corn reproductive growth stages)
#'  
#'  
#'  
#'  Note: other alternative is to take 14 days after emergence of silks as end of the flowering phase
#'  
#' @return 



estimate_flowering_phase<-function(data_to_predict,data_estimation_model=NULL,GDD_place_to_predict){
  
  ## Estimation of silk emergence for each hybrid (as date)
  if (length(data_estimation_model)>1){
  data_estimation_model$GDU=as.numeric(as.vector(data_estimation_model$GDU))
  unique_pedigree=unique(data_to_predict$Pedigree)
  average_GDD_to_flower=vector()
  silk_emergence_date=vector()
  flowering_phase_end=vector()
  
  for (ped in unique_pedigree) {
    average_GDD_to_flower[ped]=mean(data_estimation_model[data_estimation_model$Pedigree==ped,'GDU'])
    silk_emergence_date[ped]=GDD_place_to_predict[which(GDD_place_to_predict$GDU>average_GDD_to_flower[ped])[1],'Date']
    n<-GDD_place_to_predict[which(GDD_place_to_predict$GDU>average_GDD_to_flower[ped])[1],'GDU']
    flowering_phase_end[ped]=GDD_place_to_predict[which(GDD_place_to_predict$GDU>n+280),'Date']
  }
  
  df<-cbind(names(average_GDD_to_flower),average_GDD_to_flower,silk_emergence_date,flowering_phase_end)
  colnames(df)<-c('Pedigree','GDD_to_silk','Date.Silking','Date.end.lag.phase')
  df<-as.data.frame(df)
  df$Date.pre.flowering<-as.numeric(as.vector(df$Date.Silking))-10
  return(df)
  }
  if(length(data_estimation_model)==0){
    data_estimation_model$GDU=as.numeric(as.vector(data_estimation_model$GDU))
    unique_pedigree=unique(data_to_predict$Pedigree)
    
    flowering_phase_end=vector()
    
    for (ped in unique_pedigree) {
      d<-data_to_predict[data_to_predict$Pedigree==ped,'Date.Silking']
      n<-GDD_place_to_predict[GDD_place_to_predict$Date==d,'GDU']
      flowering_phase_end[ped]=GDD_place_to_predict[which(GDD_place_to_predict$GDU>n+280),'Date']
    }
    df<-cbind(names(flowering_phase_end),flowering_phase_end)
    colnames(df)<-c('Pedigree','Date.end.lag.phase')
    df<-as.data.frame(df)
    df$Date.pre.flowering<-as.numeric(as.vector(df$Date.Silking))-10
    return(df)
  }
  
  
}
