## Read imputed genotype matrix
setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING")
imputed_geno=fread('/home/uni08/jubin1/LinkImpute/g2f/procedure1_imputed_G2F.txt',header=F)
markers<-read.table('markers_before_imputation_procedure.txt')
lines<-readRDS('name_kept_lines_procedure.RDS')
colnames(imputed_geno)<-as.character(as.vector(markers$V1))


geno<-cbind(lines,imputed_geno)
colnames(geno)[1]<-'line'


saveRDS(geno,'geno_after_imputation.rds')
