##### Code author: C. Jubin #######
##### Pre-processing of G2F genotype data #######
args = commandArgs(trailingOnly = TRUE)
x1 = args[1]
x2 = args[2]
x3 = args[3]
x4 = args[4]
x5 = args[5]
x6 = args[6]
c_gwdg = as.integer(Sys.getenv('SLURM_CPUS_ON_NODE'))

library(data.table)
library(doParallel)
library(purrr)
library(stringr)
source(
  '/home/uni08/jubin1/Data/GenomesToFields/QTL_ML/spaeml_nam/data/Rfunctions_maf/recode_snps.R'
)
source(
  '/home/uni08/jubin1/Data/GenomesToFields/QTL_ML/spaeml_nam/data/Rfunctions_maf/compute_maf.R'
)
`%notin%` <- Negate(`%in%`)

heterozygotes = c('R', 'Y', 'S', 'W', 'K', 'M')
tri_allelic = c('B', 'D', 'H', 'V')
homozygotes = c('A', 'C', 'G', 'T')
any_base = c('N', '.', '-', '0')

df=readRDS('replaced_val.RDS')
print(dim(df))
print(colnames(df))
print(df[,1])
map <- read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING/datasets/pos_list_g2f.txt',header=T)
name_snps=map$Name
print(map$Position[500])
print(colnames(df)[501])
print(map$Position[25000])
print(colnames(df)[25001])

#colnames(df)[2:ncol(df)]<-as.character(as.vector(name_snps)) 



df=as.matrix(df)
print(row.names(df))

row.names(df)=NULL

# ------------------------------------------------------------------------------
# Step1: Missing call rate per G2F line
# ------------------------------------------------------------------------------

library(doParallel)
print('Lines processing start:')

lines_processing = function(geno_table = df[, 2:ncol(df)],
                            threshold_na = as.numeric(x1),
                            threshold_het = as.numeric(x2)) {
  
  
  
  
  missing_per_line <- rowSums(is.na(geno_table)) / ncol(geno_table)
  
  heterozygosity_per_line <-
    rowSums(geno_table=='R' |geno_table=='Y'|geno_table=='S'|geno_table=='W'|geno_table=='K'|geno_table=='M' ,na.rm = T) / ncol(geno_table)
  
  to_remove = unique(c(which(
    missing_per_line > threshold_na  ), 
    which(heterozygosity_per_line > threshold_het)
  ))  
  df_final = df[which(1:nrow(df) %notin% to_remove),]
  
  
  l = list(missing_per_line, heterozygosity_per_line, df_final)
  return(l)
}


l = lines_processing()

cat('Step 2 is achieved:lines processing done')

df_final = l[[3]]
saveRDS(
  df_final[,1],
  paste(
    'name_kept_lines_',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)

setwd("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/GENOTYPE_PROCESSING")

saveRDS(
  list(l[[1]],l[[2]]),
  paste(
    'Lines',
    as.numeric(x1),
    'NA_',
    as.numeric(x2),
    'Het',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)
rm(df)

# ------------------------------------------------------------------------------
# Step2: Removing SNPs which are monomorphic, and missing (>10% of the inbred lines) or heterozygous (>5% of the inbred lines)
# ------------------------------------------------------------------------------

snps_processing = function(geno_table = df_final[, 2:ncol(df_final)],
                           threshold_na = as.numeric(x3),
                           threshold_het = as.numeric(x4)) {
  monomorphic <-
    apply(geno_table, 2, function(x) {
      length(unique(x[!is.na(x)]))
    })
  missing_per_snp <-
    colSums(is.na(geno_table)) / nrow(geno_table)
  
  
  heterozygosity_per_snp <-
    colSums(geno_table=='R' |geno_table=='Y'|geno_table=='S'|geno_table=='W'|geno_table=='K'|geno_table=='M' ,na.rm = T) / nrow(geno_table)
  
  print(names(which(monomorphic == 1)))
  print(names(which(missing_per_snp > threshold_na)))
  print(names(which(heterozygosity_per_snp > threshold_het)))
  
  to_remove = unique(c(names(which(monomorphic == 1)), names(which(
    missing_per_snp > threshold_na
  )), names(
    which(heterozygosity_per_snp > threshold_het)
  )))
  
  print(length(to_remove))

  to_keep=colnames(geno_table)[colnames(geno_table)%notin%to_remove]
  
  to_keep=subset(geno_table,select = to_keep)    
  
  l = list(monomorphic,
           missing_per_snp,
           heterozygosity_per_snp,
           to_keep)
  return(l)
}



dat = snps_processing()
cat('Step 3 is achieved.')

df_final_2 = dat[[4]]

saveRDS(
  colnames(df_final_2),
  paste(
    'name_kept_markers_',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)

saveRDS(
  df_final_2,
  paste(
    'geno_matrix_after_filtering_beforerenaming_',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)


saveRDS(
  list(dat[[1]],dat[[2]],dat[[3]]),
  paste(
    'SNPs',
    as.numeric(x3),
    'NA_',
    as.numeric(x4),
    'Het',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)

# ------------------------------------------------------------------------------
# Step3: For simplicity, all the heterozygous genotypes were assigned as missing and subsequently imputed. Matrix coded numerically for further imputation with LinkImpute.
# ------------------------------------------------------------------------------


df_final_2[df_final_2 == 'A'] <- 'AA'
df_final_2[df_final_2 == 'G'] <- 'GG'
df_final_2[df_final_2 == 'C'] <- 'CC'
df_final_2[df_final_2 == 'T'] <- 'TT'
df_final_2[df_final_2 == 'R'] <- NA
df_final_2[df_final_2 == 'Y'] <- NA
df_final_2[df_final_2 == 'S'] <- NA
df_final_2[df_final_2 == 'W'] <- NA
df_final_2[df_final_2 == 'K'] <- NA
df_final_2[df_final_2 == 'M'] <- NA


monomorphic <-
  apply(df_final_2, 2, function(x) {
    length(unique(x[!is.na(x)]))
  })
to_remove2 = unique(c(names(which(monomorphic == 1))))
print('Second removal of SNPs monomorphic after assigning het as NA')
print(length(to_remove2))
df_final_2 = df_final_2[, colnames(df_final_2) %notin% to_remove2]
print(dim(df_final_2))


write.table(
  df_final_2,
  file = paste(
    'after_filtering_',
    as.character(x6),
    '.RDS',
    sep = ''
  )
  ,
  col.names = F,
  sep = ' ',
  row.names = F
)

write.table(
  colnames(df_final_2),
  file = paste(
    'markers_before_imputation_',
    as.character(x6),
    '.txt',
    sep = ''
  )
  ,
  col.names = F,
  sep = ' ',
  row.names = F
)


print(dim(df_final_2))



geno_maf= compute_maf(
  x = df_final_2,
  output = c("marker_maf"),
  missing_value = NA,
  maf_threshold = as.numeric(x5)
)
saveRDS(
  geno_maf,
  paste(
    'maf_marker',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)

geno_lst = compute_maf(
  x = df_final_2,
  output = c("geno_list"),
  missing_value = NA,
  maf_threshold = as.numeric(x5)
)

major <- geno_lst[["major_genotype"]]

minor <- geno_lst[["minor_genotype"]]





saveRDS(
  list(major,minor),
  paste(
    'major_minor_',
    as.character(x6),
    '.RDS',
    sep = ''
  )
)


# Recode the major allele as '2', the minor allele as '0' and heterozygotes as '1'.
coded_m = recode_snps(
  df_final_2,
  major = major,
  minor = minor,
  major_coding = 0,
  minor_coding = 2,
  het_coding = 1,
  missing_value = "??",
  na_coding = NA_real_
)

## calc_n
n0 <- apply(coded_m==0,2,sum,na.rm=T)
n1 <- apply(coded_m==1,2,sum,na.rm=T)
n2 <- apply(coded_m==2,2,sum,na.rm=T)

n <- n0 + n1 + n2

## calculate allele frequencies
p <- ((2*n0)+n1)/(2*n)
q <- 1 - p
maf <- pmin(p, q)
mgf <- apply(cbind(n0,n1,n2),1,min) / n
to_remove=names(which(maf<as.numeric(x5)))
print(length(to_remove))
coded_m=coded_m[,colnames(coded_m)%notin%to_remove]
print(dim(coded_m))

coded_m[is.na(coded_m)] <- -1
coded_m[coded_m==1] <- -1

print(dim(coded_m))
print(sum(coded_m==-1))

write.table(
  coded_m,
  file = paste(
    '/home/uni08/jubin1/LinkImpute/',
    as.character(x6),
    'G2F_preprocessed',
    '.txt',
    sep = ''
  ),
  col.names = F,
  sep = ' ',
  row.names = F
)

