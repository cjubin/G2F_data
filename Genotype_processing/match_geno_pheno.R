`%notin%` <- Negate(`%in%`)
library(dplyr)
library(data.table)
library(gtools)
library(plyr)


## Identify the inbred lines used in hybrids combinations in the genotypic file
## Some names differ slightly.


# Load names of inbred lines present in pheno files.
all_inbreds_phenos=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/PHENOTYPES_PROCESSING/datasets/all_inbreds_pheno_trials.txt',header = T,sep = '\t',stringsAsFactors = F,quote = NULL)


# Load genotypes names together with GBS sample ID
genotypes_names=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/g2f_2017_gbs_hybrid_codes.txt',header=T,fill = T,sep = '\t')
df1<-as.data.frame(genotypes_names[,c('Female.Pedigree','Female.GBS')])
df2<-as.data.frame(genotypes_names[,c('Male.Pedigree','Male.GBS')])
colnames(df1)<-colnames(df2)<-c('Pedigree','GBSname')
genotypes_names=smartbind(df1,df2)
genotypes_names=unique(genotypes_names)


phenos_to_identify_in_genos=all_inbreds_phenos$all_inbreds[which(all_inbreds_phenos$all_inbreds%notin%genotypes_names$Pedigree)]


# First, need to eliminate set of genotypes from 2018 genotypes with no GBS data + unknown phenotyped hybrids not identified in genotypic files.

toMatch<- c('CHECK','DKC','NILAS')
matches <- unique (grep(paste(toMatch,collapse="|"), 
                        phenos_to_identify_in_genos, value=TRUE))
phenos_remaining_not_identified<-phenos_to_identify_in_genos[-which(phenos_to_identify_in_genos%in%matches)]


all_gbs_samples_table=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/taxa_list_g2f.txt',header = T,sep = '\t')



# Read phenotypes identified in original Zea mays GBS file (pedigree --> GBS sample)

phenos_remaining_not_identified <- gsub("PHW65_MOG", "PHW65_MoG",phenos_remaining_not_identified )
phenos_remaining_not_identified <- gsub("MO44_PHW65", "Mo44_PHW65",phenos_remaining_not_identified )
phenos_remaining_not_identified <- gsub("DK78010", "78010",phenos_remaining_not_identified )

phenos_remaining_not_identified <- phenos_remaining_not_identified[which(phenos_remaining_not_identified%notin%all_gbs_samples_table$DNASample)]


# Final list: Names of phenotypes with genotypic data
phenos_with_SNPsdata=all_inbreds_phenos$all_inbreds[all_inbreds_phenos$all_inbreds%notin%phenos_remaining_not_identified]

GBS_samples1<-as.character(as.vector(all_gbs_samples_table[match(phenos_with_SNPsdata,all_gbs_samples_table$DNASample),'Taxa']))
GBS_samples2<-as.character(as.vector(genotypes_names[match(phenos_with_SNPsdata,genotypes_names$Pedigree),'GBSname']))

phenos_to_genos<-as.data.frame(cbind(phenos_with_SNPsdata,GBS_samples1,GBS_samples2))

phenos_to_genos$GBS_samples3=NA
for (j in 1:nrow(phenos_to_genos)) {
  if (!is.na(phenos_to_genos$GBS_samples1)){
    phenos_to_genos[j,'GBS_samples3']<-as.character(phenos_to_genos$GBS_samples1[j])
  } else{phenos_to_genos[j,'GBS_samples3']<-as.character(phenos_to_genos$GBS_samples2[j])}
}

phenos_to_genos<-phenos_to_genos[,c('phenos_with_SNPsdata','GBS_samples3')]
phenos_to_genos<-as.matrix(phenos_to_genos)
phenos_to_genos[phenos_to_genos=='']<-NA
phenos_to_genos<-phenos_to_genos[complete.cases(phenos_to_genos[,c('phenos_with_SNPsdata','GBS_samples3')]),]
colnames(phenos_to_genos)<-c('Pedigree','GBSname')
phenos_to_genos<-rbind(phenos_to_genos,read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/additional_6inbreds.txt',header = T))
phenos_to_genos$Pedigree<-as.character(as.vector(phenos_to_genos$Pedigree))
phenos_to_genos <- phenos_to_genos[order(phenos_to_genos$Pedigree),]
write.table(phenos_to_genos,file='phenos_to_genos.txt',col.names=T,row.names=F,sep='\t',quote=F)

##

gbs_2017_samples=readRDS('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/g2f_GBS_2017.rds')
gbs_2017_samples=as.data.frame(gbs_2017_samples)
gbs_6samples=fread('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/gbs_6samples.txt.txt',header = T,sep='\t')
gbs_6samples=as.data.frame(gbs_6samples)
colnames(gbs_2017_samples)[1]='Taxa'
colnames(gbs_6samples)[1]='Taxa'

g2f_2017_markers=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/pos_list_g2f.txt',header = T)
g2f_Zea_markers=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/position_snps_list_zea.txt',header = T)
colnames(gbs_6samples)[2:ncol(gbs_6samples)]<-as.character(as.vector(g2f_Zea_markers$Name))


colnames(gbs_2017_samples)[2:ncol(gbs_2017_samples)]<-as.character(as.vector(g2f_2017_markers$Name))
gbs_final1<-gbs_2017_samples[which(gbs_2017_samples$Taxa%in%phenos_to_genos$GBSname),which(colnames(gbs_2017_samples)%in%g2f_2017_markers$Name)]
names_taxa_1<-gbs_2017_samples$Taxa[which(gbs_2017_samples$Taxa%in%phenos_to_genos$GBSname)]
gbs_final2<-gbs_6samples[,which(colnames(gbs_6samples)%in%g2f_2017_markers$Name)]
names_taxa_2<-gbs_6samples$Taxa
gbs_final<-bind_rows(gbs_final1,gbs_final2)
names_taxa<-c(names_taxa_1,names_taxa_2)
gbs_final<-cbind(as.character(as.vector(names_taxa)),gbs_final)
colnames(gbs_final)[1]<-'line'
saveRDS(gbs_final,'gbs_final.rds')
write.table(gbs_final,'gbs_final.txt',col.names=T,row.names=F,sep='\t',quote=F)



