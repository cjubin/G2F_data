rm(list = ls())
library(lme4)
library(lsmeans)
library(plyr)


source("compute_redres.R")
source("pearsonres.R")
source("plot_redres.R")
source("rawres.R")
source("redres.R")
source("stdres.R")


`%notin%` <- Negate(`%in%`)
dat <- read.csv("g2f_2015_hybrid_data_clean.csv", header=T, stringsAsFactors=F)
dat$Year=2015
dat$Plot.Discarded..enter..yes..or..blank..=as.character(as.vector(dat$Plot.Discarded..enter..yes..or..blank..))
dat=dat[-which(dat$Plot.Discarded..enter..yes..or..blank..=='Yes'),]



## Peculiarities regarding sowing date?

for (s in unique(dat$Field.Location)) {
  if (length(unique(dat[dat$Field.Location==s,'Date.Planted']))>1){
    print('Date.Planted')
    print(s)
    print(unique(dat[dat$Field.Location==s,'Date.Planted']))
  }
  if (length(unique(dat[dat$Field.Location==s,'Date.Harvested']))>1){
    print('Date.Harvested')
    print(s)
    print(unique(dat[dat$Field.Location==s,'Date.Harvested']))
    #print(table(dat[dat$Field.Location==s,'Date.Harvested']))
  }
}

## Removed pheno. observations based on collaborators comments
critical_comments<-c("gap","inbred/hybrid mixture","large gaps","gaps", "0.75 inbred; 0.25 hybrid","half plot","most plot missing","combine ??","Half of the row was empty" )
matches <- grep(paste(critical_comments,collapse="|"), 
                dat$Comments, value=TRUE)
dat<-dat[-which(dat$Comments%in%matches),]


traitCols <- c('Year',"Field.Location", "Pedigree","Replicate",'Block', "Pollen.DAP..days.", "Silk.DAP..days.", "Plant.Height..cm.", "Ear.Height..cm.", "Grain.Yield..bu.A.")
set <- dat$Trial
dat <- dat[,traitCols]
colnames(dat) <- c("year",'location', "pedigree", "rep",'block', "pollendap", "silkdap", "pltht", "earht", "yld_bu_ac")

write.table(dat,file = 'pheno2015_after_cleaning.txt',col.names = T,row.names = F,sep = '\t',quote = F)



dat$pedigree=as.factor(as.character(dat$pedigree))
dat$rep=as.factor(as.character(dat$rep))
dat$block=as.factor(as.character(dat$block))


# Pedigrees were randomly assigned to reps within locations



means1=list()




j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$yld_bu_ac),]) == 0) {
    next
  }
  else{
    print(v)
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$yld_bu_ac), 'rep'])) == 1) {
      print('Only 1 rep')
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$yld_bu_ac), 'block'])) > 1) {
        print('case 1')
        mod1 = lmer(
          yld_bu_ac ~ pedigree + (1 |
                                    block),
          data = dat[dat$location == v &
                       !is.na(dat$yld_bu_ac),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        
        
        
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      
      else if (length(unique(dat[dat$location == v &
                                 !is.na(dat$yld_bu_ac), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                                 !is.na(dat$yld_bu_ac),])) {
        print('case 2')
        mod1 = lm(
          yld_bu_ac ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$yld_bu_ac),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        
        
        
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      else{
        print('case 3')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$yld_bu_ac), c(year, v, 'pedigree', 'yld_bu_ac')]
      }
      j = 1 + j
    } else{
      mod1 = lmer(
        yld_bu_ac ~ pedigree + (1 |
                                  rep),
        data = dat[dat$location == v &
                     !is.na(dat$yld_bu_ac),],
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
      )
      
      
      sc_resids <- compute_redres(mod1, type = "std_cond")
      ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
      means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
      
      j = 1 + j
    }
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
  else{
    print(v)
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$pltht), 'rep'])) == 1) {
      print('Only 1 rep')
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$pltht), 'block'])) > 1) {
        print('case 1')
        mod1 = lmer(
          pltht ~ pedigree + (1 |
                                block),
          data = dat[dat$location == v &
                       !is.na(dat$pltht),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      
      else if (length(unique(dat[dat$location == v &
                                 !is.na(dat$pltht), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                             !is.na(dat$pltht),])) {
        print('case 2')
        mod1 = lm(
          pltht ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$pltht),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      else{
        print('case 3')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$pltht), c(year, field.location, 'pedigree', 'pltht')]
      }
      j = 1 + j
    } else{
      mod1 = lmer(
        pltht ~ pedigree + (1 |
                              rep),
        data = dat[dat$location == v &
                     !is.na(dat$pltht),],
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
      )
      sc_resids <- compute_redres(mod1, type = "std_cond")
      ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
      means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
      
      j = 1 + j
    }
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
  else{
    print(v)
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$earht), 'rep'])) == 1) {
      print('Only 1 rep')
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$earht), 'block'])) > 1) {
        print('case 1')
        mod1 = lmer(
          earht ~ pedigree + (1 |
                                block),
          data = dat[dat$location == v &
                       !is.na(dat$earht),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      
      else if (length(unique(dat[dat$location == v &
                                 !is.na(dat$earht), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                             !is.na(dat$earht),])) {
        print('case 2')
        mod1 = lm(
          earht ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$earht),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      else{
        print('case 3')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$earht), c(year, field.location, 'pedigree', 'earht')]
      }
      j = 1 + j
    } else{
      mod1 = lmer(
        earht ~ pedigree + (1 |
                              rep),
        data = dat[dat$location == v &
                     !is.na(dat$earht),],
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
      )
      sc_resids <- compute_redres(mod1, type = "std_cond")
      ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
      means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
      
      j = 1 + j
    }
  }
}  

earht1=ldply(means1, rbind)
colnames(earht1)=c('year','location','pedigree','earht')


## Pollen dap

means1=list()




j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$pollendap),]) == 0) {
    next
  }
  else{
    print(v)
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$pollendap), 'rep'])) == 1) {
      print('Only 1 rep')
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$pollendap), 'block'])) > 1) {
        print('case 1')
        mod1 = lmer(
          pollendap ~ pedigree + (1 |
                                    block),
          data = dat[dat$location == v &
                       !is.na(dat$pollendap),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      
      else if (length(unique(dat[dat$location == v &
                                 !is.na(dat$pollendap), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                                 !is.na(dat$pollendap),])) {
        print('case 2')
        mod1 = lm(
          pollendap ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$pollendap),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      else{
        print('case 3')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$pollendap), c(year, field.location, 'pedigree', 'pollendap')]
      }
      j = 1 + j
    } else{
      mod1 = lmer(
        pollendap ~ pedigree + (1 |
                                  rep),
        data = dat[dat$location == v &
                     !is.na(dat$pollendap),],
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
      )
      sc_resids <- compute_redres(mod1, type = "std_cond")
      ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
      means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
      
      j = 1 + j
    }
  }
}  

pollendap1=ldply(means1, rbind)
colnames(pollendap1)=c('year','location','pedigree','pollendap')


# Silk dap



means1=list()




j = 1
for (v in unique(dat$location)) {
  if (nrow(dat[dat$location == v & !is.na(dat$silkdap),]) == 0) {
    next
  }
  else{
    print(v)
    if (length(unique(dat[dat$location == v &
                          !is.na(dat$silkdap), 'rep'])) == 1) {
      print('Only 1 rep')
      if (length(unique(dat[dat$location == v &
                            !is.na(dat$silkdap), 'block'])) > 1) {
        print('case 1')
        mod1 = lmer(
          silkdap ~ pedigree + (1 |
                                  block),
          data = dat[dat$location == v &
                       !is.na(dat$silkdap),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      
      else if (length(unique(dat[dat$location == v &
                                 !is.na(dat$silkdap), 'pedigree'])) < nrow(dat[dat$location == v &
                                                                               !is.na(dat$silkdap),])) {
        print('case 2')
        mod1 = lm(
          silkdap ~ pedigree,
          data = dat[dat$location == v &
                       !is.na(dat$silkdap),],
          control = lmerControl(check.conv.singular = .makeCC(
            action = "ignore",  tol = 1e-4
          ))
        )
        means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
        
      }
      else{
        print('case 3')
        means1[[j]] = dat[dat$location == v &
                            !is.na(dat$silkdap), c(year, field.location, 'pedigree', 'silkdap')]
      }
      j = 1 + j
    } else{
      mod1 = lmer(
        silkdap ~ pedigree + (1 |
                                rep),
        data = dat[dat$location == v &
                     !is.na(dat$silkdap),],
        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
      )
      sc_resids <- compute_redres(mod1, type = "std_cond")
      ind<-sc_resids[abs(sc_resids$res)>3,'pedigree']
      means1[[j]] = cbind('2015', v, as.data.frame(lsmeans(mod1, "pedigree"))[, c(1:2)])
      means1[[j]]<-means1[[j]][which(means1[[j]]$pedigree%notin%ind),]
      
      j = 1 + j
    }
  }
}  

silkdap1=ldply(means1, rbind)
colnames(silkdap1)=c('year','location','pedigree','silkdap')





#Merge multiple data frames with the 5 traits (all rows full)


data1=merge(gy1,pltht1,by=c('pedigree','location'),all.x = T,all.y = T)
data2=merge(data1,earht1,by=c('pedigree','location'),all.x = T,all.y = T)
data3=merge(data2,silkdap1,by=c('pedigree','location'),all.x = T,all.y = T)
data4=merge(data3,pollendap1,by=c('pedigree','location'),all.x = T,all.y = T)
df_final<-data4[,c('pedigree','location','silkdap','pollendap','yld_bu_ac','pltht','earht')]
df_final<-cbind(df_final[,1:2],year=2015,df_final[,3:ncol(df_final)])
write.table(df_final,file = 'pheno2015.txt',col.names = T,row.names = F,sep = '\t',quote = F)
