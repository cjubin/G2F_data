#' Derive environmental covariates from daily weather data based on estimated growth stages (using flowering tame, planting and harvest date)
#'
#' \code{derive_W_growth_stages} 
#'
#' @param 
#' @param 
#' 
#' 
#' @return 
#' @examples
#'  
#'
#' @importFrom magrittr %>%
#' @export
#' 
#' 

`%notin%` <- Negate(`%in%`)

derive_W_growth_stages=function(version_days=T,use_ref_season_length=T,length_ref=145,version_GDD=F,weather_file='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/replaced_daily_weather.txt',pheno_file='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/phenotypes_after_processing_only_with_GBS.txt'){
 
  if (version_days){
  weather=read.table(weather_file,header = T,sep = '\t')
  
  phenos=read.table(pheno_file,header = T,sep = '\t')
  phenos$Year_Exp=as.character(phenos$Year_Exp)
  weather$Year_Exp=as.character(as.vector(weather$Year_Exp))
  
  #Round days needed to determine temporal intervals corresponding to growth stages. 
  phenos$pollendap=ceiling(phenos$pollendap)
  phenos$silkdap=ceiling(phenos$silkdap)
  phenos$location=as.character(phenos$location)
  
  ##Create a column referring to the weather data index : same weather data to use but with different planting/harvest dates
  
  
  phenos[phenos$location=='IAH1a'&phenos$year==2014,'Year_Exp']='2014_IAH1'
  phenos[phenos$location=='IAH1b'&phenos$year==2014,'Year_Exp']='2014_IAH1'
  phenos[phenos$location=='IAH1c'&phenos$year==2014,'Year_Exp']='2014_IAH1'
  
  
  #############
  
  add1<-weather[which(weather$Year_Exp=='2016_ILH1'),]
  add1$Year_Exp='2016_ILH1.b'
  add1$Date.Planted=117
  weather[weather$Year_Exp=='2016_ILH1','Year_Exp']<-'2016_ILH1.a'
  weather[which(weather$Year_Exp=='2016_ILH1.a'),'Date.Planted']<-127
  
  add2<-weather[which(weather$Year_Exp=='2017_TXH1'),]
  add2$Year_Exp='2017_TXH1-Dry'
  add2$Date.Planted=62
  add2$Date.Harvested=206
  add3<-weather[which(weather$Year_Exp=='2017_TXH1'),]
  add3$Year_Exp='2017_TXH1-Early'
  add3$Date.Planted=62
  add3$Date.Harvested=212
  
  weather[weather$Year_Exp=='2017_TXH1','Year_Exp']<-'2017_TXH1-Late'
  weather[which(weather$Year_Exp=='2016_ILH1.a'),'Date.Planted']<-96
  weather[which(weather$Year_Exp=='2016_ILH1.a'),'Date.Harvested']<-222
  
  
  add4<-weather[which(weather$Year_Exp=='2018_KSH1'),]
  add4$Year_Exp='2018_KSH1.drought'
  weather[weather$Year_Exp=='2018_KSH1','Year_Exp']<-'2018_KSH1.irrigated'
  
  weather=rbind(weather,add1,add2,add3,add4)
  
  

  weather[weather$Year_Exp=='2018_TXH1- Dry','Year_Exp']<-'2018_TXH1-Dry'
  weather[weather$Year_Exp=='2018_TXH1- Early','Year_Exp']<-'2018_TXH1-Early'
  weather[weather$Year_Exp=='2018_TXH1- Late','Year_Exp']<-'2018_TXH1-Late'
  
  weather = plyr::arrange(weather, Year_Exp, Day.of.Year)
  phenos = plyr::arrange(phenos, Year_Exp,pedigree)
  #############
  
  phenos=phenos[!is.na(phenos$silkdap),]
  
  #############
  
  phenos$P.V=NA
  phenos$P.F=NA
  phenos$P.G=NA
  
  phenos$FreqP5.V=NA
  phenos$FreqP5.F=NA
  phenos$FreqP5.G=NA
  
  phenos$Wind.V=NA
  phenos$Wind.F=NA
  phenos$Wind.G=NA
  
  phenos$MeanT.V=NA
  phenos$MeanT.F=NA
  phenos$MeanT.G=NA
  
  phenos$MinT.V=NA
  phenos$MinT.F=NA
  phenos$MinT.G=NA
  
  phenos$MaxT.V=NA
  phenos$MaxT.F=NA
  phenos$MaxT.G=NA
  
  phenos$GDD.V=NA
  phenos$GDD.F=NA
  phenos$GDD.G=NA
  
  
  phenos$FreqMaxT30.V=NA
  phenos$FreqMaxT30.F=NA
  phenos$FreqMaxT30.G=NA
  
  phenos$FreqMaxT35.V=NA
  phenos$FreqMaxT35.F=NA
  phenos$FreqMaxT35.G=NA
  
  
  phenos$EvTot.V=NA
  phenos$EvTot.F=NA
  phenos$EvTot.G=NA
  
  phenos$Max.consecutive.dry.days.V=NA
  phenos$Max.consecutive.dry.days.F=NA
  phenos$Max.consecutive.dry.days.G=NA
  
  phenos$Photoperiod.hrs.Tot.V=NA
  phenos$Photoperiod.hrs.Tot.F=NA
  phenos$Photoperiod.hrs.Tot.G=NA  
  
  phenos$Sdrad.V=NA
  phenos$Sdrad.F=NA
  phenos$Sdrad.G=NA
  
  phenos$RatioP.ET0.V=NA
  phenos$RatioP.ET0.F=NA
  phenos$RatioP.ET0.G=NA
  
  
  for (i in unique(phenos$Year_Exp)){
    
    
    
    for (s in unique(as.numeric(as.vector(phenos[phenos$Year_Exp==i,'silkdap'])))){
      
      if(is.na(s)){next}
      start=unique(weather[weather$Year_Exp==i,'Date.Planted'])
      end=unique(weather[weather$Year_Exp==i,'Date.Harvested'])
      length_gs=end-start
      print(length_gs)
      
      o=start+s
      c=length_gs/length_ref
      # Flowering time period: 3 weeks in total
      FD1=o-7
      FD2=o+14
      if (use_ref_season_length){PM=round(FD2+60*c)}
      
      
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'wind.speed.ms.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMAX'])
      #phenos[phenos$Year_Exp==i,]$St32.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'TMAX'])
      #phenos[phenos$Year_Exp==i,]$St35.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1&weather$TMAX>35,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'et0'])
      #phenos[phenos$Year_Exp==i,]$Nb.dry.days.V=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'photothermal_time'])
      #print(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'PRCP2']/weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start+7)&weather$Day.of.Year<FD1,'solar_radiation_NASA'])
      
      
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.ms.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St32.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St35.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2&weather$TMAX>35,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Nb.dry.days.F=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'photothermal_time'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'solar_radiation_NASA'])
      
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'wind.speed.ms.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St32.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St35.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM&weather$TMAX>35,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Nb.dry.days.G=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'photothermal_time'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'solar_radiation_NASA'])
      
      
    }
    
    
  }
  }
  
  return(phenos)
}


pheno=derive_W_growth_stages()

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals")
write.table(pheno,'phenos_with_EC_covariates.txt',col.names = T,row.names=F,sep="\t",quote = F)

