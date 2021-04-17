#########################################################################################
import numpy as np
import pandas as pd
import random
import math
import statistics 
import SKAT
import scipy.stats 
import statsmodels.api as sm
#########################################################################################
###Burden test
def Burden_test(geno_mat, weights, phenotype):
    
    C = geno_mat.dot(weights)
    C_mean = C.values - np.mean(C)
    Y_mean = phenotype.values - np.mean(phenotype)
    u = C_mean.T.dot(Y_mean)
    varu = np.var(Y_mean)*C_mean.T.dot(C_mean)
    Q = pow(u,2) / varu
    pvalue = scipy.stats.distributions.chi2.sf(x = Q, df = 1) 
    
    return pvalue

###Prepare data for summary stat
def summary_stat_data(geno_mat, phenotype):
    geno_center = geno_mat - np.mean(geno_mat)
    Y_mean = phenotype.values - np.mean(phenotype)
    nsnp = np.shape(geno_center)[1]
    beta = np.zeros(nsnp)
    std = np.zeros(nsnp)
    
    for i in range(nsnp):
        fit = sm.OLS(Y_mean,geno_center.iloc[:,i]).fit()
        beta[i] = fit.params
        std[i] = fit.bse
    
    return beta, std

#########################################################################################
###set up phenotype heritability
pheno_h2 = 0.2
stdev = np.sqrt(1 - pheno_h2)

test_dosage_path = "/sample_data/simulation/test_dosage.txt"
true_expression_path = "/sample_data/simulation/expression_test.txt"
DPR_weight_path = "/sample_data/simulation/DPR_weights.param.txt"
PrediXcan_weight_path = "/sample_data/simulation/prediXcan_weights.txt"
gwas_SS_path = "/sample_data/simulation/model1_GWAS.txt"

###load test_dosage
test_dosage = pd.read_csv(test_dosage_path, header = None,sep = "\t")
test_dosage_mat = test_dosage.loc[:,5:].T

###load DPR weights
dpr_weight = pd.read_csv(DPR_weight_path, sep = "\t")
###the final estimated effect size from DPR is b + beta 
dpr_weight['ES'] = dpr_weight['b'] + dpr_weight['beta']


###set up threshold. In simulation, take the mean of median of final estimated effect size in each simulation
threshold = np.median(abs(dpr_weight['ES']))
dpr_weight_filtered  = dpr_weight[abs(dpr_weight.ES) > threshold]
test_dosage_mat_filtered = test_dosage_mat.loc[:,abs(dpr_weight.ES) > threshold]


###load PrediXcan weights
PrediXcan_weight = pd.read_csv(PrediXcan_weight_path, header = None, sep = "\t")
###get PrediXcan dosage
test_dosage_mat_PrediXcan = test_dosage_mat.loc[:,test_dosage.loc[:,1].isin(PrediXcan_weight.loc[:,1])]


###load true expression
test_expression = pd.read_csv(true_expression_path, header = None,sep = "\t")

###simulated phenotype from gene expression level
var_exp = np.var(test_expression)
gamma = np.sqrt(pheno_h2 / var_exp)
error = np.random.normal(0, stdev, len(test_expression)).reshape(len(test_expression), 1)
###exp_simulated = gamma * train_expression(true) + error
exp_simulated = gamma * test_expression + error
exp_simulated = exp_simulated.iloc[:,0] #phenotype format pd.DataFrame, shape is(n,)

###VC_TWAS approach
pvalue_SKAT_dpr=SKAT.SKAT(test_dosage_mat,exp_simulated, None, dpr_weight['ES'], 'C')
pvalue_SKAT_dpr_f=SKAT.SKAT(test_dosage_mat_filtered,exp_simulated, None, dpr_weight_filtered['ES'], 'C')
pvalue_SKAT_PrediXcan=SKAT.SKAT(test_dosage_mat_PrediXcan,exp_simulated, None, PrediXcan_weight.loc[:,5], 'C')

##Burden TWAS approach
pvalue_burden_dpr = Burden_test(test_dosage_mat, dpr_weight['ES'], exp_simulated)
pvalue_burden_dpr_f = Burden_test(test_dosage_mat_filtered, dpr_weight_filtered['ES'], exp_simulated)
pvalue_burden_PrediXcan = Burden_test(test_dosage_mat_PrediXcan, PrediXcan_weight.loc[:,5].values, exp_simulated)


##pepare data for summary stat simulation
beta, std = summary_stat_data(test_dosage_mat, exp_simulated)
SNP_info = test_dosage.loc[:,1:4]
SNP_info['BETA'] = beta
SNP_info['SE'] = std
SNP_info.columns = ['rsID','POS','REF','ALT','BETA','SE']


SNP_info.to_csv(
    gwas_SS_path, index=None, sep='\t', mode='w')


