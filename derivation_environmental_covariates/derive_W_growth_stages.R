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

derive_W_growth_stages=function(version_days=T,use_ref_season_length=F,length_ref=145,version_GDD=F,weather_file='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/replaced_daily_weather.txt',pheno_file='/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/phenotypes_after_processing_only_with_GBS.txt'){
 
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
  
  
  
  phenos = plyr::arrange(phenos, Year_Exp,pedigree)
  
  phenos$Planting.Date=weather[match(phenos$Year_Exp,weather$Year_Exp),'Date.Planted']
  phenos$SandProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Sand']
  phenos$SiltProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Silt']
  phenos$ClayProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Clay']
  phenos$OM.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'OM']
  #phenos$city
  phenos$counties=weather[match(phenos$Year_Exp,weather$Year_Exp),'county']
  phenos$state=weather[match(phenos$Year_Exp,weather$Year_Exp),'state']
  phenos_to_pred= phenos[which(is.na(phenos$silkdap)),c('pedigree','Year_Exp','Planting.Date')]
  
  
  write.table(phenos_to_pred,'phenos_to_predict.txt',col.names = T,row.names = F,sep = '\t',quote = F)
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
  phenos$GDD.at.silking=NA
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
  
  #phenos$Max.consecutive.dry.days.V=NA
  #phenos$Max.consecutive.dry.days.F=NA
  #phenos$Max.consecutive.dry.days.G=NA
  
  phenos$Photoperiod.hrs.Tot.V=NA
  phenos$Photoperiod.hrs.Tot.F=NA
  phenos$Photoperiod.hrs.Tot.G=NA  
  
  phenos$Sdrad.V=NA
  phenos$Sdrad.F=NA
  phenos$Sdrad.G=NA
  
  #phenos$RatioP.ET0.V=NA
  #phenos$RatioP.ET0.F=NA
  #phenos$RatioP.ET0.G=NA
  
  phenos$SumDiffP_ETP.V=NA
  phenos$SumDiffP_ETP.F=NA
  phenos$SumDiffP_ETP.G=NA
  
  phenos$growing_season_length=NA
  phenos$coeff_ref_length_season=NA
  
  
  for (i in unique(phenos$Year_Exp)){
    
    
    
    for (s in unique(as.numeric(as.vector(phenos[phenos$Year_Exp==i,'silkdap'])))){
      
      if(is.na(s)){next}
      
      
      start=unique(weather[weather$Year_Exp==i,'Date.Planted'])
      end=unique(weather[weather$Year_Exp==i,'Date.Harvested'])
      length_gs=end-start
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,'growing_season_length']<-length_gs
      
      o=start+s
      c=length_gs/length_ref
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,'coeff_ref_length_season']<-c
      print(c)
      # Flowering time period: 3 weeks in total
      FD1=o-7
      FD2=o+14
      # Grain filling time lasts about 60 days for 140 days growing season reference length
      if (use_ref_season_length){PM=FD2+ceiling(60*c)}
      else{PM=FD2+60}
      
      # 3 time periods considered:
      # 1. From planting time to 1 week before start of silking (start to FD1)
      # 2. From 1 week before start of silking to 2 weeks after start of silking (FD1 to FD2)
      # 3. Based on a reference season length (but set as option in the function): about 60 days after end of flowering period (FD2 to PM)
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      #phenos[phenos$Year_Exp==i,]$St32.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      #phenos[phenos$Year_Exp==i,]$St35.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1&weather$TMAX>35,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'et0'],na.rm = T)
      #phenos[phenos$Year_Exp==i,]$Nb.dry.days.V=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'photothermal_time'])
      #print(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2']/weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'et0'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$SumDiffP_ETP.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'daily_diffPrec_ETP.V'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'solar_radiation_NASA'],na.rm = T)
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.at.silking=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<o,'GDD'])
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'et0'],na.rm=T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St32.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St35.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2&weather$TMAX>35,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Nb.dry.days.F=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'photothermal_time'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$SumDiffP_ETP.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'daily_diffPrec_ETP.F'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'solar_radiation_NASA'],na.rm = T)
      
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Wind.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$EvTot.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'et0'],na.rm=T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St32.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM&weather$TMAX>32,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St35.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM&weather$TMAX>35,'TMAX'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Nb.dry.days.G=sum(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'diff.PminusETP']<0))
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'photothermal_time'])
      #phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$RatioP.ET0.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'PRCP2'])/sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'et0'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$SumDiffP_ETP.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'daily_diffPrec_ETP.G'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<PM,'solar_radiation_NASA'],na.rm = T)
      
      
    }
    
    
  }
  }
  
  return(phenos)
}


pheno=derive_W_growth_stages()

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Growth_stages_intervals")
write.table(pheno,'phenos_with_EC_covariates.txt',col.names = T,row.names=F,sep="\t",quote = F)

