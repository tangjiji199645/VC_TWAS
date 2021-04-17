##########################################################################################
import numpy as np
import pandas as pd
import SKAT
#########################################################################################
############################################################
##SNP in LD, GWAS, DPR weights are in the same order
############################################################
LD_path = "/sample_data/simulation/ld.txt"
gwas_path = "/sample_data/simulation/model2_GWAS.txt"
DPR_weight_path = "/sample_data/simulation/DPR_weights.param.txt"

###load LD matrix
ld = pd.read_csv(LD_path, sep="\t", header=None)
###load GWAS result
gwas=pd.read_csv(gwas_path, sep="\t")


###load DPR weights
dpr_weight = pd.read_csv(DPR_weight_path, sep = "\t")
###the final estimated effect size from DPR is b + beta 
dpr_weight['ES'] = dpr_weight['b'] + dpr_weight['beta']

###set up threshold. In simulation, take the mean of median of final estimated effect size in each simulation
threshold = np.median(abs(dpr_weight['ES']))
dpr_weight_filtered  = dpr_weight[abs(dpr_weight.ES) > threshold]
dpr_weight_filtered_index = dpr_weight_filtered.index

###get GWAS result for filtered DPR weights
gwas_f=gwas[abs(dpr_weight.ES) > threshold]

###get Covariance matrix and diag Matrix of Covariance matrix
V = ld.values
D = np.diag(V)

V_f = ld.loc[dpr_weight_filtered_index, dpr_weight_filtered_index].values
D_f = np.diag(V_f)

###get gwas result
beta_estimate_all = gwas["BETA"].values
beta_var_all=pow(gwas["SE"].values,2)

beta_estimate_f = gwas_f["BETA"].values
beta_var_f=pow(gwas_f["SE"].values,2)


###VC_TWAS using summary stat
###SKAT_summary(beta_var, beta_estimate, weight, sample_size, COV, D)
p_val_all=SKAT.SKAT_summary(beta_var_all, beta_estimate_all, dpr_weight['ES'], 800, V, D)
p_val_f=SKAT.SKAT_summary(beta_var_f, beta_estimate_f, dpr_weight_filtered['ES'], 800, V_f, D_f)

