#!/bin/bash
#SBATCH -p medium 
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --mem=90G
module load intel/mkl/64/2017/2.174
echo $PATH
module load openmpi/intel/64/1.10.7
echo $PATH
Rscript 4_run_kriging_all_variables_e2.R 2
Rscript 4_run_kriging_all_variables_e2.R 3
Rscript 4_run_kriging_all_variables_e2.R 4
Rscript 4_run_kriging_all_variables_e2.R 5
Rscript 4_run_kriging_all_variables_e2.R 6
Rscript 4_run_kriging_all_variables_e2.R 7
Rscript 4_run_kriging_all_variables_e2.R 8
Rscript 4_run_kriging_all_variables_e2.R 9
Rscript 4_run_kriging_all_variables_e2.R 10
Rscript 4_run_kriging_all_variables_e2.R 11
Rscript 4_run_kriging_all_variables_e2.R 12
Rscript 4_run_kriging_all_variables_e2.R 13
Rscript 4_run_kriging_all_variables_e2.R 14
Rscript 4_run_kriging_all_variables_e2.R 15
Rscript 4_run_kriging_all_variables_e2.R 16
Rscript 4_run_kriging_all_variables_e2.R 17
Rscript 4_run_kriging_all_variables_e2.R 18
Rscript 4_run_kriging_all_variables_e2.R 19
Rscript 4_run_kriging_all_variables_e2.R 20
Rscript 4_run_kriging_all_variables_e2.R 21
