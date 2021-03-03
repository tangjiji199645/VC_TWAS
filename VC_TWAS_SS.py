###################################################################
## for more details, please see https://github.com/yanglab-emory/TIGAR
###################################################################
# Import packages needed
import numpy as np
import pandas as pd
import SKAT

###############################################################
#load weight
weight = pd.read_csv(weight_path,sep='\t')
weight_mat = weight["BETA"].values
#load reference covariance matrix, assume the order of SNP in weight,reference covariance matrix and gwas result are the same
MCOV = pd.read_csv(MCOV_path,sep='\t')
MCOV = MCOV.values
#get the diag of reference covariance matrix
D= np.diag(MCOV)
#load reference GWAS result, assume the order of SNP in weight, reference covariance matrix and gwas result are the same
GWAS_result = pd.read_csv(GWAS_path,sep='\t')
GWAS_beta = GWAS_result["BETA"].values
#sample size of summary stat
sample_size=50000
#the value of y'y, phenotype of summary stat, if unknown, please use SKAT.SKAT_summary
y_est=600000
#known y'y
p_val = SKAT.SKAT_summary_withy(y_est, GWAS_beta, weight_mat, args.sample_size, MCOV, D)
#unknown y'y
GWAS_var = pow(GWAS_result["SE"].values,2)
p_val = SKAT.SKAT_summary_withy(GWAS_var, GWAS_beta, weight_mat, args.sample_size, MCOV, D)

