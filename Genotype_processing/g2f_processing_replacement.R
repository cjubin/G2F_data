##### Code author: C. Jubin #######
##### Pre-processing of NAM genotype data #######
args = commandArgs(trailingOnly = TRUE)
x1 = args[1]
x2 = args[2]
x3 = args[3]
x4 = args[4]
x5 = args[5]
x6 = as.character(args[6])
c_gwdg = as.integer(Sys.getenv('SLURM_CPUS_ON_NODE'))

library(data.table)
library(doParallel)
library(purrr)
source(
  '/home/uni08/jubin1/Data/GenomesToFields/QTL_ML/spaeml_nam/data/Rfunctions_maf/recode_snps.R'
)
source(
  '/home/uni08/jubin1/Data/GenomesToFields/QTL_ML/spaeml_nam/data/Rfunctions_maf/compute_maf.R'
)
`%notin%` <- Negate(`%in%`)

# ------------------------------------------------------------------------------
# Step1: Load table with the 5874 taxa belonging to G2F populations
# ------------------------------------------------------------------------------




df=readRDS('gbs_final.rds')
print(df[,1])
print(dim(df))
## IUPAC nucleotide code

heterozygotes = c('R', 'Y', 'S', 'W', 'K', 'M')
tri_allelic = c('B', 'D', 'H', 'V')
homozygotes = c('A', 'C', 'G', 'T')
any_base = c('N', '.', '-', '0')



# ------------------------------------------------------------------------------
# Step2: Removing G2F lines with more than 10% NA and 5% heterozygous rates (per line)
# ------------------------------------------------------------------------------

## Convert unknown genotypic values to NA, and tri-allelic values to NA
df = as.matrix(df)
pos_snps = colnames(df) 

df <- apply(df, 2, function(x) {
  x[x %in% any_base] <- NA
  x
})
print('First replacement done')
df<- apply(df, 2, function(x) {
  x[x %in% tri_allelic] <- NA
  x
})
print('Second replacement done')

saveRDS(df,'replaced_val.RDS')
quit()
