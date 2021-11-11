`%notin%` <- Negate(`%in%`)

derive_W_growth_stages=function(version_days=T,use_ref_season_length=F,length_ref=145,version_GDD=F,weather_file='replaced_daily_weather.txt',pheno_file='phenotypes_after_processing2.txt'){
 
  if (version_days){
  weather=read.table(weather_file,header = T,sep = '\t')
  
  phenos=read.table(pheno_file,header = T,sep = '\t')
  phenos$Year_Exp=as.character(phenos$Year_Exp)
  weather$Year_Exp=as.character(as.vector(weather$Year_Exp))
  
  #Round days needed to determine temporal intervals corresponding to growth stages. 
  phenos$pollendap=ceiling(phenos$pollendap)
  phenos$silkdap=ceiling(phenos$silkdap)
  phenos$location=as.character(phenos$location)
  
  
  
  #############
  
  
  
  phenos = plyr::arrange(phenos, Year_Exp,pedigree)
  
  phenos$Planting.Date=weather[match(phenos$Year_Exp,weather$Year_Exp),'Date.Planted']
  phenos$Harvest.Date=weather[match(phenos$Year_Exp,weather$Year_Exp),'Date.Harvested']
  phenos$SandProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Sand']
  phenos$SiltProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Silt']
  phenos$ClayProp.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'X..Clay']
  phenos$Longitude=weather[match(phenos$Year_Exp,weather$Year_Exp),'long']
  phenos$Latitude=weather[match(phenos$Year_Exp,weather$Year_Exp),'lat']


  phenos$OM.SC=weather[match(phenos$Year_Exp,weather$Year_Exp),'OM']
  #phenos$city
  phenos$counties=weather[match(phenos$Year_Exp,weather$Year_Exp),'county']
  phenos$state=weather[match(phenos$Year_Exp,weather$Year_Exp),'state']
  phenos_to_pred= phenos[which(is.na(phenos$silkdap)),c('pedigree','Year_Exp','Planting.Date')]
  
  
  write.table(phenos_to_pred,'phenos_to_predict.txt',col.names = T,row.names = F,sep = '\t',quote = F)
  #############
  
  #phenos=phenos[!is.na(phenos$silkdap),]
  
  #############
  
  phenos$P.V=NA
  phenos$P.F=NA
  phenos$P.G=NA
  
  phenos$FreqP5.V=NA
  phenos$FreqP5.F=NA
  phenos$FreqP5.G=NA
  
  phenos$Mean.Wind.V=NA
  phenos$Mean.Wind.F=NA
  phenos$Mean.Wind.G=NA
  
  phenos$MaxWind.V=NA
  phenos$MaxWind.F=NA
  phenos$MaxWind.G=NA
  
  phenos$FreqDaysWindSpeedSup12.V=NA
  phenos$FreqDaysWindSpeedSup12.F=NA
  phenos$FreqDaysWindSpeedSup12.G=NA
  
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
  
  phenos$St30.V=NA
  phenos$St30.F=NA
  phenos$St30.G=NA
  
  phenos$CumSumET0.V=NA
  phenos$CumSumET0.F=NA
  phenos$CumSumET0.G=NA

  
  #phenos$Max.consecutive.dry.days.V=NA
  #phenos$Max.consecutive.dry.days.F=NA
  #phenos$Max.consecutive.dry.days.G=NA
  
  phenos$Photothermal.time.Tot.V=NA
  phenos$Photothermal.time.Tot.F=NA
  phenos$Photothermal.time.Tot.G=NA  
  
  phenos$Sdrad.V=NA
  phenos$Sdrad.F=NA
  phenos$Sdrad.G=NA
  
  phenos$CumDailyWaterBalance.V=NA
  phenos$CumDailyWaterBalance.F=NA
  phenos$CumDailyWaterBalance.G=NA
    
  phenos$growing_season_length=NA
  phenos$coeff_ref_length_season=NA
  weather$water_balance=weather$PRCP2-weather$et0
  
  weather$etp=NA
  
  for (i in unique(phenos$Year_Exp)){
  print(i)
    
    
    for (s in unique(as.numeric(as.vector(phenos[phenos$Year_Exp==i,'silkdap'])))){
      print(s)
      if(is.na(s)){next}
      
      
      start=unique(phenos[phenos$Year_Exp==i&phenos$silkdap==s,'Planting.Date'])
      #print(start)
      end=unique(phenos[phenos$Year_Exp==i&phenos$silkdap==s,'Harvest.Date'])
      #print(end)
      length_gs=end-start
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,'growing_season_length']<-length_gs
      
      o=start+s
      c=length_gs/length_ref
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,'coeff_ref_length_season']<-c
      
      
      # Flowering time period: 3 weeks in total
      FD1=o-7
      FD2=o+14

    
      print(weather[weather$Year_Exp==i,'et0'])

      
      # Grain filling time lasts about 65 days for 140 days growing season reference length
      if (use_ref_season_length){
        PM=FD2+ceiling(65*c)
        if(PM>end){print('the harvest was done before +65 days after FW.')
          PM=end}}
      else{PM=FD2+65}
      
      # For other covariates: 3 time periods considered:
      # 1. From planting time to 1 week before start of silking (start to FD1)
      # 2. From 1 week before start of silking to 2 weeks after start of silking (FD1 to FD2)
      # 3. Based on a reference season length (but set as option in the function): about 65 days after end of flowering period (FD2 to PM)
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Mean.Wind.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxWind.V=max(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.V=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      phenos[phenos$Year_Exp==i,]$St30.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1&weather$TMAX>30,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumSumET0.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'et0'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumDailyWaterBalance.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'water_balance'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photothermal.time.Tot.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'photothermal_time'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.V=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'solar_radiation_NASA'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqDaysWindSpeedSup12.V=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'wind.speed.m.s']>12))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<FD1,'wind.speed.m.s'])
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.at.silking=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=(start)&weather$Day.of.Year<o,'GDD'])
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Mean.Wind.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxWind.F=max(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.F=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumSumET0.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'et0'],na.rm=T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumDailyWaterBalance.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'water_balance'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St30.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2&weather$TMAX>30,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photothermal.time.Tot.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'photothermal_time'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.F=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'solar_radiation_NASA'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqDaysWindSpeedSup12.F=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.m.s']>12))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD1&weather$Day.of.Year<FD2,'wind.speed.m.s'])
      
      
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$P.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqP5.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'PRCP2']>5))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'PRCP2'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Mean.Wind.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxWind.G=max(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'wind.speed.m.s'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MeanT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMEAN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MinT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMIN'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$MaxT.G=mean(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$GDD.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'GDD'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumSumET0.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'et0'],na.rm=T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$CumDailyWaterBalance.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'water_balance'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT30.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMAX']>30))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$St30.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM&weather$TMAX>30,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqMaxT35.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMAX']>35))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'TMAX'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Photothermal.time.Tot.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'photothermal_time'])
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$Sdrad.G=sum(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'solar_radiation_NASA'],na.rm = T)
      phenos[phenos$Year_Exp==i&phenos$silkdap==s,]$FreqDaysWindSpeedSup12.G=length(which(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'wind.speed.m.s']>12))/length(weather[weather$Year_Exp==i&weather$Day.of.Year>=FD2&weather$Day.of.Year<=PM,'wind.speed.m.s'])
      
    }
    
    
  }
  }
  
  return(phenos)
}


pheno=derive_W_growth_stages()


write.table(pheno,'phenos_with_EC_covariates_without_univ_factor.txt',col.names = T,row.names=F,sep="\t",quote = F)


## Add the university factor, similar to locations but more useful because it gives a hint at the geographical localization in the US.
## If different locations were used for one university, then an integer is added at the end of the factor. This integer refers to the same location i fused several times throughout years of experiment.
factoruniversity=read.table('factoruniversity.txt',header = T,sep='\t')
pheno=merge(pheno,factoruniversity[,c(1,4,5,6)],by='Year_Exp',all.x=T)


write.table(pheno,'phenos_with_EC_covariates.txt',col.names = T,row.names=F,sep="\t",quote = F)


