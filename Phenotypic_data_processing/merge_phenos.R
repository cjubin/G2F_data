rm(list = ls())
library(data.table)
library(stringr)
library(stringi)
phenos2014=read.table('pheno2014.txt',header = T,sep = '\t',)
pheno2014_after_cleaning=read.csv('pheno2014_after_cleaning.csv',header = T)
phenos2014$year<-2014
phenos2015=read.table('pheno2015.txt',header = T,sep = '\t')
pheno2015_after_cleaning=read.table('pheno2015_after_cleaning.txt',header = T,sep = '\t')
phenos2015$year<-2015
phenos2016=read.table('pheno2016.txt',header = T,sep = '\t')
pheno2016_after_cleaning=read.table('pheno2016_after_cleaning.txt',header = T,sep = '\t')
phenos2016$year<-2016
phenos2017=read.table('pheno2017.txt',header = T,sep = '\t')
pheno2017_after_cleaning=read.table('pheno2017_after_cleaning.txt',header = T,sep = '\t')
phenos2017$year<-2017
phenos2018=read.table('pheno2018.txt',header = T,sep = '\t')
pheno2018_after_cleaning=read.table('pheno2018_after_cleaning.txt',header = T,sep = '\t')
phenos2018$year<-2018


require(plyr)
unique(ddply(phenos2014, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])
unique(ddply(pheno2014_after_cleaning, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])

unique(ddply(phenos2015, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])
unique(ddply(pheno2015_after_cleaning, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])

unique(ddply(phenos2016, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])
unique(ddply(pheno2016_after_cleaning, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])

unique(ddply(phenos2017, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])
unique(ddply(pheno2017_after_cleaning, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])

unique(ddply(phenos2018, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])
unique(ddply(pheno2018_after_cleaning, .(location), mutate, count = length(unique(pedigree)))[,c('location','count')])

all_phenos<-rbind(phenos2014,phenos2015,phenos2016,phenos2017,phenos2018)
all_phenos$Year_Exp<-paste0(all_phenos$year,'_',all_phenos$location)


### First step: remove phenotypes which do not have grain yield data ###

all_phenos_1<-all_phenos[which(!is.na(all_phenos$yld_bu_ac)),]


### Second step: remove hybrids of environments which should be discarded: ###

# NYH1 2015, 2016 and 2017 --> disease trials REMOVED
# TXH2 2015, 2016, 2017 and 2018 REMOVED
# ARH1+ARH2 2016; GAH1 2016 + ARH1 2018 removed --> many dates with irrigation present but no records on the amount REMOVED
# MNH1 2017: lots of critical comments regarding the whole trial experiment : REMOVED
# Notes in g2f_field_metadata for 2017_MNH1: "Northwest portion of field was underwater for part of the early growing season, affecting phenotypic data.",Yeild in font of field much lower than the back of the field due to water pooling early in the season.,"Stand counts were taken after harvest and before tilling. Some plots were impossible to count, therefor no data are provided.""
# Note: NEH1 indicated in metadata as irrigated but no further info; only 2 dates with irrigation for TXH1_2014, no amount tracked --> those environments were however kept

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("NYH1")&all_phenos_1$year%in%c(2015,2016,2017)),]

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("GAH1")&all_phenos_1$year%in%c(2016)),]

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("MNH1")&all_phenos_1$year%in%c(2017)),]

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("TXH2")&all_phenos_1$year%in%c(2015,2016,2017,2018)),]

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("ARH1","ARH2")&all_phenos_1$year%in%c(2016)),]

all_phenos_1<-all_phenos_1[-which(all_phenos_1$location%in%c("ARH1")&all_phenos_1$year%in%c(2018)),]



### Third step: remove hybrids for which no genotypic data of the parents are available ###

# Load list of genotypes retained after filtering (with samples names)
samples_kept<-readRDS('name_kept_lines_procedure.RDS')
key_pheno_geno<-read.table('phenos_to_genos.txt',header=T,sep = '\t')
key_pheno_geno<-key_pheno_geno[key_pheno_geno$GBSname%in%samples_kept,]

# Add some phenotypes to the key in Pedigree --> some phenotype names slightly differ but refer to the same genotype (cf script match_geno_pheno.R)
supp1<-key_pheno_geno[grep('PHW65_MOG',key_pheno_geno$Pedigree),]
supp1$Pedigree<- gsub("PHW65_MOG", "PHW65_MoG",supp1$Pedigree)
supp2<-key_pheno_geno[grep('MO44_PHW65',key_pheno_geno$Pedigree),]
supp2$Pedigree<- gsub("MO44_PHW65", 'Mo44_PHW65',supp2$Pedigree)
key_pheno_geno<-rbind(key_pheno_geno,supp1,supp2)


all_phenos_1$parent1<-stri_extract_first(all_phenos_1$pedigree, regex = '[^/]+')
all_phenos_1$parent2<-stri_extract_last(all_phenos_1$pedigree, regex = '[^/]+')
all_phenos_1$genotyped<-NA
all_phenos_1$parent1.GBS.sample<-NA
all_phenos_1$parent2.GBS.sample<-NA

options(warn = 0)
for (i in 1:nrow(all_phenos_1)) {
  if (all_phenos_1$parent1[i]%in%key_pheno_geno$Pedigree&all_phenos_1$parent2[i]%in%key_pheno_geno$Pedigree){
    print(i)
    all_phenos_1[i,'genotyped']<-'YES'
    all_phenos_1$parent1.GBS.sample[i]<-unique(as.character(key_pheno_geno[key_pheno_geno$Pedigree==all_phenos_1$parent1[i],'GBSname']))
    all_phenos_1$parent2.GBS.sample[i]<-unique(as.character(key_pheno_geno[key_pheno_geno$Pedigree==all_phenos_1$parent2[i],'GBSname']))
    
    
    
  }else{all_phenos_1[i,'genotyped']<-'NO'}
  
}


### Final step: construct table to obtain synthetic genotypes from parental inbreds using Pedigree ### 
write.table(all_phenos_1[all_phenos_1$genotyped=='YES',],file='phenotypes_after_processing_all.txt',col.names = T,row.names = F,sep = '\t',quote = F)
write.table(all_phenos_1[all_phenos_1$genotyped=='YES',],file='phenotypes_after_processing_only_with_GBS.txt',col.names = T,row.names = F,sep = '\t',quote = F)



write.table(unique(all_phenos_1[all_phenos_1$genotyped=='YES',c("pedigree","parent1","parent2","parent1.GBS.sample","parent2.GBS.sample")]),file='table_hybrids_to_generate.txt',col.names = T,row.names = F,sep = '\t',quote = F)

