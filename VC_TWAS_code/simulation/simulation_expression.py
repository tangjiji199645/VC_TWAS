#########################################################################################
import numpy as np
import pandas as pd
import random
import math
import statistics 
#########################################################################################
###For more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR
#########################################################################################

#########################################################################################
###Load train dosage and test dosage, the order of SNP are the same in train dosage and test dosage
###First 5 rows is chrom, rsID, POS, REF, ALT
train_dosage_path = "/sample_data/simulation/train_dosage.txt"
test_dosage_path = "/sample_data/simulation/test_dosage.txt"
true_beta_path = "/sample_data/simulation/true_beta.txt"
true_expression_test_path = "/sample_data/simulation/expression_test.txt"
true_expression_train_path = "/sample_data/simulation/expression_train.txt"

train_dosage = pd.read_csv(train_dosage_path, header = None, sep ="\t")
n_train = len(train_dosage.columns) - 5

test_dosage = pd.read_csv(test_dosage_path, header = None,sep = "\t")
n_test = len(test_dosage.columns) - 5

#Merge train and test genotype together
train_test_merge = pd.concat([train_dosage.loc[:,5:], test_dosage.loc[:,5:]], axis=1)

##Set up gene expression heritability and causal probability
heritability = 0.2
causal_prop = 0.01
stdev = math.sqrt(1 - heritability)


##Simulate number of causal snp 
n_snp = round(causal_prop * len(train_dosage))
causal_snp = random.sample(range(1000),n_snp)
beta_mat = np.random.normal(0, 1, n_snp)

##Generate gene expression level based on true causal and heritability
exp_inital = train_test_merge.loc[causal_snp].values.T.dot(beta_mat)
gamma = math.sqrt(heritability / statistics.variance(exp_inital))
beta_gamma_mat = gamma * beta_mat
error_term = np.random.normal(0, stdev, len(train_test_merge.columns))
expression = gamma * exp_inital + error_term

##Combined true beta data
true_beta_info = train_dosage.loc[causal_snp,1:4]
true_beta_info['beta'] = beta_gamma_mat
true_beta_info.columns = ['rsID','POS','REF','ALT','BETA']

##True gene expression for train and test 
expression_train = expression[0:n_train]
expression_test = expression[n_train:]

##Save output
true_beta_info.to_csv(
    "/sample_data/simulation/true_beta.txt", index=None, sep='\t', mode='w')

np.savetxt(true_expression_train_path, expression_train)
np.savetxt(true_expression_test_path, expression_test)


