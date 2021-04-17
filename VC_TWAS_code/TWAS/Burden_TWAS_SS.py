#########################################################################################
import numpy as np
import pandas as pd
import SKAT
from scipy.stats import chi2
#########################################################################################
###For more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR
#########################################################################################

def get_V_cor(V_cov):
    V_cov = V_cov.copy()
    v = np.sqrt(np.diag(V_cov))
    outer_v = np.outer(v, v)
    V_cor = V_cov / outer_v
    V_cor[V_cov == 0] = 0
    return V_cor

def get_z_denom(V, w):
    return np.sqrt(np.linalg.multi_dot([w, V, w]))

def get_spred_zscore(V_cov, w, Z_gwas, snp_sd):
    Z_twas = snp_sd.dot(w * Z_gwas) / get_z_denom(V_cov, w)
    return Z_twas, get_pval(Z_twas)
    
def get_fusion_zscore(V_cov, w, Z_gwas, snp_sd=None):
    V_cor = get_V_cor(V_cov)
    Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
    return Z_twas, get_pval(Z_twas)
 
def get_pval(z): return np.format_float_scientific(1-chi2.cdf(z**2, 1), precision=15, exp_digits=0)

#########################################################################################
###Load GWAS result, LD matrix and estimated weights, the order of SNP are the same in these file

gwas_path = "/sample_data/TWAS/GWAS.txt"
ld_path = "/sample_data/TWAS/ld.txt"
weight_path = "/sample_data/TWAS/DPR_weights.txt"


gwas = pd.read_csv(gwas_path,sep='\t')
gwas['Zscore'] = gwas['BETA'] / gwas['SE']

ld = pd.read_csv(ld_path,sep='\t')
V = ld.values
weight = pd.read_csv(weight_path,sep='\t')


get_zscore_args = [V, gwas['BETA'].values, gwas['Zscore'].values, gwas['SE'].values]

##Burden test, FUSION
Z_FUSION, pval_FUSION = get_fusion_zscore(*get_zscore_args)

##Burden test, SPrediXcan
Z_SPred, pval_SPred = get_spred_zscore(*get_zscore_args)

