#########################################################################################
import numpy as np
import pandas as pd
import SKAT

#########################################################################################
###For more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR
#########################################################################################

#########################################################################################
###Load gwas, weights, and ld, the order of SNP are the same in these files
gwas_path = "/sample_data/TWAS/GWAS.txt"
ld_path = "/sample_data/TWAS/ld.txt"
weight_path = "/sample_data/TWAS/DPR_weights.txt"

gwas = pd.read_csv(gwas_path,sep='\t')
ld = pd.read_csv(ld_path,sep='\t')
weight = pd.read_csv(weight_path,sep='\t')

V = ld.values
D = np.diag(V)

beta = gwas["BETA"].values
beta_var=pow(gwas["SE"].values,2)

p_val=SKAT.SKAT_summary(beta_var, beta, weight['ES_sum'], 128, V, D)
