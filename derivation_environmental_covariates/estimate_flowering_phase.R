#' Find approximate dates for flowering phase, i.e. from
#' initiation of silks to the date for which abortion does not occur anymore,
#' based on the weather data (only thermal time)
#' Especially needed when silking was not scored.
#'
#' \code{estimate_flowering_phase} returns a data.frame with approximate start and end of flowering time for each hybrid (i.e. pedigree column),
#' based on data from other experiments with these hybrids and using weather data at the location to predict.
#'
#' @param data_to_predict data.frame with 3 or 4 columns:
#' - 3 columns: both begin and end of flowering phase to predict. First column 'pedigree', second column the 'Year_Exp' ID, third column the 'Planting.Date '(numeric, in Day.of.Year) 
#' - 4 columns: only end of flowering phase to predict (silking was scored). Fourth column for BLUEs of silking dates: 'Date.Silking'
#' @param data_estimation_model data.frame with 6 columns containing first column 'Pedigree', second column the 'Year_Exp' ID, 
#' third column the 'Planting.Date' (numeric, in Day.of.Year), 
#' fourth column the 'Silking.Date' (as Day.of.Year, numeric), 
#' fifth column the 'GDD' reached at Silk DAP, 
#'
#' 
#' @param GDD_place_to_predict data.frame with colnames 'Date' and 'GDD'. 
#' The data.frame contains in first column all dates from Planting to Harvest dates corresponding to the places to predict.
#' Second column contains the corresponding cumulative GDDs (real weather data from the concerned place).
#' 
#' 
#' Estimation of the end of the flowering phase --> 280 GDDs after beginning of silking (R1 + R2 corn reproductive growth stages)
#'  
#'  
#'  
#'  Note: other alternative is to take 14 days after emergence of silks as end of the flowering phase
#'  
#' @return 



estimate_flowering_phase<-function(i,data_to_predict='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals/Growth_stages_intervals/phenos_to_predict.txt',data_estimation_model='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals/data_estimation_model.txt',GDD_place_to_predict='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals/GDD_place_to_predict.txt'){
  data_estimation_model=read.table(data_estimation_model,header = T,sep = '\t')
  
  data_to_predict=read.table(data_to_predict,header = T,sep = '\t')
  data_to_predict$Year_Exp=as.character(data_to_predict$Year_Exp)
  data_to_predict=data_to_predict[which(data_to_predict$Year_Exp==i),]
  
  GDD_place_to_predict=read.table(GDD_place_to_predict,header = T,sep = '\t')
  GDD_place_to_predict=GDD_place_to_predict[which(GDD_place_to_predict$Year_Exp==i),]
  GDD_place_to_predict$cumGDD= cumsum(GDD_place_to_predict$GDD)
  
  
  ## Estimation of silk emergence for each hybrid (as date)
  if (nrow(data_estimation_model)>1){
    data_estimation_model$GDD=as.numeric(as.vector(data_estimation_model$GDD))
    unique_pedigree=unique(data_to_predict$pedigree)
    average_GDD_to_flower=vector()
    silk_emergence_date=vector()
    flowering_phase_end=vector()
    
    for (ped in unique_pedigree) {
      average_GDD_to_flower[ped]=mean(data_estimation_model[data_estimation_model$Pedigree==ped,'GDD'])
      silk_emergence_date[ped]=GDD_place_to_predict[which(GDD_place_to_predict$cumGDD>average_GDD_to_flower[ped])[1],'Day.of.Year']
      #n<-GDD_place_to_predict[which(GDD_place_to_predict$cumGDD>average_GDD_to_flower[ped])[1],'cumGDD']
      #flowering_phase_end[ped]=GDD_place_to_predict[which(GDD_place_to_predict$cumGDD>n+280),'Date']
    }
    
    df<-cbind(names(average_GDD_to_flower),average_GDD_to_flower,silk_emergence_date)
    colnames(df)<-c('Pedigree','GDD_to_silk','Date.Silking')
    df<-as.data.frame(df)
    return(df)
  }
  
  if(length(data_estimation_model)==0){
    data_estimation_model$GDD.at.silking=as.numeric(as.vector(data_estimation_model$GDD.at.silking))
    unique_pedigree=unique(data_to_predict$Pedigree)
    
    flowering_phase_end=vector()
    
    for (ped in unique_pedigree) {
      d<-data_to_predict[data_to_predict$Pedigree==ped,'Date.Silking']
      n<-GDD_place_to_predict[GDD_place_to_predict$Date==d,'cumGDD']
      flowering_phase_end[ped]=GDD_place_to_predict[which(GDD_place_to_predict$cumGDD>n+280),'Date']
    }
    df<-cbind(names(flowering_phase_end),flowering_phase_end)
    colnames(df)<-c('Pedigree','Date.end.lag.phase')
    df<-as.data.frame(df)
    df$Date.pre.flowering<-as.numeric(as.vector(df$Date.Silking))-10
    return(df)
  }
  
  
}
phenos_to_predict=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals/Growth_stages_intervals/phenos_to_predict.txt',header = T,sep = '\t')
table(phenos_to_predict$Year_Exp)
phenos_to_predict[phenos_to_predict$Year_Exp%in%c('2014_IAH2','2014_IAH3'),]
GDD_place_to_predict=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals/GDD_place_to_predict.txt',header = T,sep = '\t')
year_exps=unique(GDD_place_to_predict$Year_Exp)
estimate_flowering_phase(i=unique(year_exps)[1])

