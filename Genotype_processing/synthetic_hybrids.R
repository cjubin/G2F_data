rm(list = ls())
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING")

hybrids_to_generate<-read.table('table_hybrids_to_generate.txt',header = T,sep = '\t')

geno<-readRDS('geno_after_imputation.rds')
geno<-as.matrix(geno)



geno_hybrids<-matrix(NA,nrow = nrow(hybrids_to_generate),ncol = ncol(geno)-1)

# x is the row number in the hybrids_to_generate table

compute_hybrids<-function(x){
  geno_parent_1<-as.numeric(geno[geno[,1]==as.character(hybrids_to_generate[x,'parent1.GBS.sample']),-1])
  geno_parent_2<-as.numeric(geno[geno[,1]==as.character(hybrids_to_generate[x,'parent2.GBS.sample']),-1])
  geno_hybrid<-rowMeans(cbind(geno_parent_1, geno_parent_2), na.rm=TRUE)
  geno_hybrid
}

res<-lapply(1:nrow(hybrids_to_generate),compute_hybrids)
res<-do.call('rbind',res)
#Note: 2,105 unique hybrids (sometimes parent2 and parent 1 are exchanged in pedigree but actually same hybrid)
colnames(res)<-colnames(geno)[2:ncol(geno)]
res2<-unique(res)

## MAF filter

m <- ncol(res2)     ## number of snps
n <- nrow(res2)     ## number of individuals

## assign all non {0,1,2} to NA

res2 <- as.matrix(res2)

## calc_n
n0 <- apply(res2==0,2,sum,na.rm=T)
n1 <- apply(res2==1,2,sum,na.rm=T)
n2 <- apply(res2==2,2,sum,na.rm=T)

n <- n0 + n1 + n2

## calculate allele frequencies
p <- ((2*n0)+n1)/(2*n)
q <- 1 - p
maf <- pmin(p, q)
mgf <- apply(cbind(n0,n1,n2),1,min) / n



kept_markers=names(maf[which(maf>0.01)])
res<-res[,kept_markers]

geno_hybrids<-cbind(hybrids_to_generate,res)

saveRDS(geno_hybrids,'geno_hybrids.rds')
write.table(geno_hybrids,'geno_hybrids.txt',col.names = T,row.names = F,sep = '\t',quote = F)
