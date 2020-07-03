
derive_W_growth_stages=function(file='C:/Users/cathyjubin/Documents/allweather_with_irrigation.txt'){
  #setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/weather_growth_stages")
  weather=read.table(file,header = T)
  
  phenos=read.table('C:/Users/cathyjubin/Documents/hybrids_no_exp_removed.txt',header = T,sep = '\t')
  
  phenos=phenos[phenos$year%in%c(2014:2018),]
  
  phenos$pollendap=round(phenos$pollendap)
  phenos$silkdap=round(phenos$silkdap)
  
  
  ##Restrict in weather data to first day of planting and last day harvesting (since sometimes the weather station was remained few days/week after...)
  weather=do.call(rbind, mget(ls(pattern = 'weather')))
  weather$loc_year=paste0(weather$location,'_',weather$year)
  
  phenos$location=as.character(phenos$location)
  
  phenos[phenos$location=='IAH1a'&phenos$year==2014,'location']='IAH1'
  phenos[phenos$location=='IAH1b'&phenos$year==2014,'location']='IAH1'
  phenos[phenos$location=='IAH1c'&phenos$year==2014,'location']='IAH1'
  phenos$loc_year=paste0(phenos$location,"_",phenos$year)
  
  ##Remove SDH1 from phenos
  phenos=phenos[-which(phenos$loc_year=='SDH1_2015'),]
  
  phenos$lon=weather[match(phenos$loc_year,weather$loc_year),'lon.x']
  phenos$lat=weather[match(phenos$loc_year,weather$loc_year),'lat.x']
  phenos$Date.Planted=weather[match(phenos$loc_year,weather$loc_year),'Date.Planted1.x']
  phenos$Date.Harvested=weather[match(phenos$loc_year,weather$loc_year),'Date.Harvested1.x']
  
  ##Remove locations with no weather data
  phenos=phenos[-which(is.na(phenos$lat)),]
  phenos$length.growing.season=phenos$Date.Harvested-phenos$Date.Planted
  
  
  #############
  
  
  
  ##
  phenos$P.EV=NA
  phenos$P.F=NA
  phenos$P.G=NA
  
  phenos$FreqP5.EV=NA
  phenos$FreqP5.F=NA
  phenos$FreqP5.G=NA
  
  phenos$Wind.EV=NA
  phenos$Wind.F=NA
  phenos$Wind.G=NA
  
  phenos$MeanT.EV=NA
  phenos$MeanT.F=NA
  phenos$MeanT.G=NA
  
  phenos$MinT.EV=NA
  phenos$MinT.F=NA
  phenos$MinT.G=NA
  
  phenos$MaxT.EV=NA
  phenos$MaxT.F=NA
  phenos$MaxT.G=NA
  
  phenos$GDD.EV=NA
  phenos$GDD.F=NA
  phenos$GDD.G=NA
  
  
  phenos$FreqMaxT30.EV=NA
  phenos$FreqMaxT30.F=NA
  phenos$FreqMaxT30.G=NA
  
  phenos$FreqMaxT35.EV=NA
  phenos$FreqMaxT35.F=NA
  phenos$FreqMaxT35.G=NA
  
  
  phenos$EvTot.EV=NA
  phenos$EvTot.F=NA
  phenos$EvTot.G=NA
  
  #phenos$Nb.dry.days.EV=NA
  #phenos$Nb.dry.days.F=NA
  #phenos$Nb.dry.days.G=NA
  
  phenos$Photoperiod.hrs.Tot.EV=NA
  phenos$Photoperiod.hrs.Tot.F=NA
  phenos$Photoperiod.hrs.Tot.G=NA  
  
  phenos$Sdrad.EV=NA
  phenos$Sdrad.F=NA
  phenos$Sdrad.G=NA
  
  phenos$Mean.Daily.RatioP.ET0.EV=NA
  phenos$Mean.Daily.RatioP.ET0.F=NA
  phenos$Mean.Daily.RatioP.ET0.G=NA
  
  
  
  for (i in unique(phenos$loc_year)){
    
    
    
    for (s in unique(as.numeric(as.vector(phenos[phenos$loc_year==i,'silkdap'])))){
      
      
      start=unique(weather[weather$loc_year==i,'Date.Planted1.x'])
      end=unique(weather[weather$loc_year==i,'Date.Harvested1.x'])
      length_gs=end-start
      print(length_gs)
      length_ref=130
      o=start+s
      c=length_gs/length_ref
      FD1=round(o-5*c) 
      FD2=round(o+10*c)
      PM=round(FD2+55*c)
      
      
      
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$P.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqP5.EV=sum(which(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'rainfall.mm']>5))/length(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Wind.EV=mean(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'wind.speed.ms.s'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MeanT.EV=mean(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'mean.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MinT.EV=mean(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'min.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MaxT.EV=mean(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$GDD.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'GDD'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT30.EV=sum(which(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'max.temperature']>30))/length(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'max.temperature'])
      #phenos[phenos$loc_year==i,]$St32.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1&weather$max.temperature>32,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT35.EV=sum(which(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'max.temperature']>35))/length(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'max.temperature'])
      #phenos[phenos$loc_year==i,]$St35.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1&weather$max.temperature>35,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$EvTot.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'ref.evapotranspiration'])
      #phenos[phenos$loc_year==i,]$Nb.dry.days.EV=sum(which(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'diff.PminusETP']<0))
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'photothermal.time'])
      #print(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'rainfall.mm']/weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Mean.Daily.RatioP.ET0.EV=mean(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'rainfall.mm']/weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Sdrad.EV=sum(weather[weather$loc_year==i&weather$day_year>=(start+7)&weather$day_year<FD1,'sum.solar.radiation.MJ.m2'])
      
      
      
      
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$P.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqP5.F=sum(which(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'rainfall.mm']>5))/length(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Wind.F=mean(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'wind.speed.ms.s'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MeanT.F=mean(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'mean.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MinT.F=mean(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'min.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MaxT.F=mean(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$GDD.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'GDD'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$EvTot.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT30.F=sum(which(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'max.temperature']>30))/length(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$St32.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2&weather$max.temperature>32,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT35.F=sum(which(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'max.temperature']>35))/length(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$St35.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2&weather$max.temperature>35,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$Nb.dry.days.F=sum(which(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'diff.PminusETP']<0))
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'photothermal.time'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Mean.Daily.RatioP.ET0.F=mean(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'rainfall.mm']/weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Sdrad.F=sum(weather[weather$loc_year==i&weather$day_year>=FD1&weather$day_year<FD2,'sum.solar.radiation.MJ.m2'])
      
      
      
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$P.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqP5.G=sum(which(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'rainfall.mm']>5))/length(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'rainfall.mm'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Wind.G=mean(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'wind.speed.ms.s'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MeanT.G=mean(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'mean.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MinT.G=mean(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'min.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$MaxT.G=mean(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$GDD.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'GDD'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$EvTot.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT30.G=sum(which(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'max.temperature']>30))/length(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$St32.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM&weather$max.temperature>32,'max.temperature'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$FreqMaxT35.G=sum(which(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'max.temperature']>35))/length(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$St35.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM&weather$max.temperature>35,'max.temperature'])
      #phenos[phenos$loc_year==i&phenos$silkdap==s,]$Nb.dry.days.G=sum(which(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'diff.PminusETP']<0))
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Photoperiod.hrs.Tot.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'photothermal.time'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Mean.Daily.RatioP.ET0.G=mean(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'rainfall.mm']/weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'ref.evapotranspiration'])
      phenos[phenos$loc_year==i&phenos$silkdap==s,]$Sdrad.G=sum(weather[weather$loc_year==i&weather$day_year>=FD2&weather$day_year<PM,'sum.solar.radiation.MJ.m2'])
      
      
    }
    
    
  }
  
  return(phenos)
}


pheno=derive_W_growth_stages()
#setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/Weather/weather_growth_stages")

write.table(pheno,'phenos_with_EC_covariates_oneEV.txt',col.names = T,row.names=F,sep="\t",quote = F)
