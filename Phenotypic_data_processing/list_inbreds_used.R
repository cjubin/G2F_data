setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/datasets')
library(stringr)
library(stringi)

ped1=list()
ped2=list()
hybrids=list()

## Retrieving all the inbred lines names used in the field experimnets 2014-2018
i=1
for (j in list.files(pattern = 'hybrid_data')){
  v=str_sub(j,5,15)
  assign(v,read.csv(j, header=T, stringsAsFactors=F,fill = T))
  hybrids[[i]]<-unique(get(v)[,'Pedigree'])
  ped1[[i]]<-stri_extract_first(get(v)[,'Pedigree'], regex = '[^/]+')
  ped2[[i]] <-stri_extract_last(get(v)[,'Pedigree'], regex = '[^/]+')
  i<-i+1
}

all_hybrids=sort(unique(c(unlist(hybrids))))
all_inbreds=sort(unique(c(unlist(ped1),unlist(ped2))))
write.table(file='all_hybrids_pheno_trials.txt',as.data.frame(all_hybrids),sep = '\t',quote = F,row.names = T,col.names = T)

write.table(file='all_inbreds_pheno_trials.txt',as.data.frame(all_inbreds),sep = '\t',quote = F,row.names = T,col.names = T)
