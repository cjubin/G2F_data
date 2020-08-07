rm(list = ls())
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/datasets")
library(lme4)
library(lsmeans)
library(plyr)


source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/compute_redres.R")
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/pearsonres.R")
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/plot_redres.R")
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/rawres.R")
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/redres.R")
source("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/outliers_functions/stdres.R")


`%notin%` <- Negate(`%in%`)
dat <- read.csv("g2f_2014_hybrid_data_clean.csv", header=T, stringsAsFactors=F)
dat$Year=2014
dat$Plot.Discarded..enter..yes..or..blank..=as.character(as.vector(dat$Plot.Discarded..enter..yes..or..blank..))
dat=dat[-which(dat$Plot.Discarded..enter..yes..or..blank..=='Yes'),]



## Peculiarities regarding sowing date?

for (s in unique(dat$Field.Location)) {
  print('Date.Planted')
  print(s)
  print(unique(dat[dat$Field.Location==s,'Date.Planted']))
  if (length(unique(dat[dat$Field.Location==s,'Date.Planted']))>1){
    #print('Date.Planted')
    #print(s)
    #print(unique(dat[dat$Field.Location==s,'Date.Planted']))
  }
  if (length(unique(dat[dat$Field.Location==s,'Date.Harvested']))>1){
    print('Date.Harvested')
    print(s)
    print(unique(dat[dat$Field.Location==s,'Date.Harvested']))
    #print(table(dat[dat$Field.Location==s,'Date.Harvested']))
  }
}

## Removed pheno. observations based on collaborators comments
critical_comments<-c("harvest error","combine ??","Sprayer damange left row","put in filler due to planting error","plant issue","TWT machine failure","planting error","Sprayer damange right row", "TWT machine failure; low stand" ,"Sprayer damange Left row","put in filler at planting due to planting error","filled in due to planting error","Sprayer damange Right row","Half of the row was empty", "planter issue","Jabbed for short seed","jabbed for short seed","wet edge of field; plot was yellow and/or border germinated poorly" )
matches <- grep(paste(critical_comments,collapse="|"), 
                dat$Comments, value=TRUE)
dat<-dat[-which(dat$Comments%in%matches),]


traitCols <- c('Year',"Field.Location", "Pedigree","Replicate",'Block', "Pollen.DAP..days.", "Silk.DAP..days.", "Plant.Height..cm.", "Ear.Height..cm.", "Grain.Yield..bu.A.")
set <- dat$Trial
dat <- dat[,traitCols]
colnames(dat) <- c("year",'location', "pedigree", "rep",'block', "pollendap", "silkdap", "pltht", "earht", "yld_bu_ac")

write.csv(dat,file = 'pheno2014_after_cleaning.csv',row.names = F)



dat$pedigree=as.factor(as.character(dat$pedigree))
dat$rep=as.factor(as.character(dat$rep))
dat$block=as.factor(as.character(dat$block))


# ANALYSIS OF 5 TRAITS



means1=list()


j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$yld_bu_ac),]) == 0) {
    next
  }
  if (length(unique(dat[dat$location == v &
                        !is.na(dat$yld_bu_ac), 'rep'])) > 1){
    print(paste0(v,' : at least 2 reps'))
    mod1 = lmer(
    yld_bu_ac ~ pedigree + (1 |
                              rep),
    data = dat[dat$location == v &
                 !is.na(dat$yld_bu_ac),])
  
  
  
  sc_resids <- compute_redres(mod1, type = "std_cond")
  ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
  means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
  means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
  
  j = 1 + j}
  else{
    print(paste0(v,' : only 1 rep'))
    
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$yld_bu_ac), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                            !is.na(dat$yld_bu_ac),])){
     
      
      
        print('case 1')
        mod1 = lm(
          yld_bu_ac ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$yld_bu_ac),]
        )
        means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
        
        
        
       
        
      }
      else{
        print('case 2')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$yld_bu_ac), c(year, field.location, 'pedigree', 'yld_bu_ac')]
      }
      j = 1 + j
    } 
  }
  

gy1=ldply(means1, rbind)
colnames(gy1)=c('year','location','pedigree','yld_bu_ac')




## Plant height 




means1=list()


j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$pltht),]) == 0) {
    next
  }
  if (length(unique(dat[dat$location == v &
                        !is.na(dat$pltht), 'rep'])) > 1){
    print(paste0(v,' : at least 2 reps'))
    mod1 = lmer(
      pltht ~ pedigree + (1 |
                                rep),
      data = dat[dat$location == v &
                   !is.na(dat$pltht),])
    
    
    
    sc_resids <- compute_redres(mod1, type = "std_cond")
    ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
    means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
    means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
    
    j = 1 + j}
  else{
    print(paste0(v,' : only 1 rep'))
    
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$pltht), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                          !is.na(dat$pltht),])){
     
        print('case 1')
        mod1 = lm(
          pltht ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$pltht),]
        )
        means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      }
      
      
      
      
      
    
    else{
      print('case 2')
      means1[[j]] = dat[dat$location == v &
                          !is.na(dat$pltht), c(year, field.location, 'pedigree', 'pltht')]
    }
    j = 1 + j
  } 
}


pltht1=ldply(means1, rbind)
colnames(pltht1)=c('year','location','pedigree','pltht')





## Ear height 




means1=list()


j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$earht),]) == 0) {
    next
  }
  if (length(unique(dat[dat$location == v &
                        !is.na(dat$earht), 'rep'])) > 1){
    print(paste0(v,' : at least 2 reps'))
    mod1 = lmer(
      earht ~ pedigree + (1 |
                            rep),
      data = dat[dat$location == v &
                   !is.na(dat$earht),])
    
    
    
    sc_resids <- compute_redres(mod1, type = "std_cond")
    ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
    means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
    means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
    
    j = 1 + j}
  else{
    print(paste0(v,' : only 1 rep'))
    
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$earht), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                      !is.na(dat$earht),])){
    
        print('case 1')
        mod1 = lm(
          earht ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$earht),]
        )
        means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      
      
      
      
      
      
    }
    else{
      print('case 2')
      means1[[j]] = dat[dat$location == v &
                          !is.na(dat$earht), c(year, field.location, 'pedigree', 'earht')]
    }
    j = 1 + j
  } 
}


earht1=ldply(means1, rbind)
colnames(earht1)=c('year','location','pedigree','earht')

## Silk dap




means1=list()


j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$silkdap),]) == 0) {
    next
  }
  if (length(unique(dat[dat$location == v &
                        !is.na(dat$silkdap), 'rep'])) > 1){
    print(paste0(v,' : at least 2 reps'))
    mod1 = lmer(
      silkdap ~ pedigree + (1 |
                            rep),
      data = dat[dat$location == v &
                   !is.na(dat$silkdap),])
    
    
    
    sc_resids <- compute_redres(mod1, type = "std_cond")
    ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
    means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
    means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
    
    j = 1 + j}
  else{
    print(paste0(v,' : only 1 rep'))
    
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$silkdap), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                      !is.na(dat$silkdap),])){
      
        print('case 1')
        mod1 = lm(
          silkdap ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$silkdap),]
        )
        means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      
      
      
      
      
      
    }
    else{
      print('case 2')
      means1[[j]] = dat[dat$location == v &
                          !is.na(dat$silkdap), c(year, field.location, 'pedigree', 'silkdap')]
    }
    j = 1 + j
  } 
}


silkdap1=ldply(means1, rbind)
colnames(silkdap1)=c('year','location','pedigree','silkdap')

## Pollen dap




means1=list()


j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$pollendap),]) == 0) {
    next
  }
  if (length(unique(dat[dat$location == v &
                        !is.na(dat$pollendap), 'rep'])) > 1){
    print(paste0(v,' : at least 2 reps'))
    mod1 = lmer(
      pollendap ~ pedigree + (1 |
                                rep),
      data = dat[dat$location == v &
                   !is.na(dat$pollendap),])
    
    
    
    sc_resids <- compute_redres(mod1, type = "std_cond")
    ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
    means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
    means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
    
    j = 1 + j}
  else{
    print(paste0(v,' : only 1 rep'))
    
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$pollendap), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                          !is.na(dat$pollendap),])){
     
        print('case 1')
        mod1 = lm(
          pollendap ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$pollendap),]
        )
        means1[[j]] = cbind('2014', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      
      
      
      
      
      
    }
    else{
      print('case 2')
      means1[[j]] = dat[dat$location == v &
                          !is.na(dat$pollendap), c(year, field.location, 'pedigree', 'pollendap')]
    }
    j = 1 + j
  } 
}


pollendap1=ldply(means1, rbind)
colnames(pollendap1)=c('year','location','pedigree','pollendap')



#Merge multiple data frames with the 5 traits (all rows full)


data1=merge(gy1,pltht1,by=c('pedigree','location'),all.x = T,all.y = T)
data2=merge(data1,earht1,by=c('pedigree','location'),all.x = T,all.y = T)
data3=merge(data2,silkdap1,by=c('pedigree','location'),all.x = T,all.y = T)
data4=merge(data3,pollendap1,by=c('pedigree','location'),all.x = T,all.y = T)
df_final<-data4[,c('pedigree','year','location','silkdap','pollendap','yld_bu_ac','pltht','earht')]
write.table(df_final,file = 'pheno2014.txt',col.names = T,row.names = F,sep = '\t',quote = F)
