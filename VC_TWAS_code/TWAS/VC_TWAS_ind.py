#########################################################################################
import numpy as np
import pandas as pd
import SKAT

#########################################################################################
###For more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR
#########################################################################################

#########################################################################################
###Load dosage, weights, and phenotype data, the order of SNP are the same in these files
###Phenotype continous 
genotype_path = "/sample_data/TWAS/genotype.txt"
phenotype_path ="/sample_data/TWAS/phenotype.txt"
###Phenotype dichotomous (0,1)
phenotype_path_D ="/sample_data/TWAS/phenotype_D.txt"
DPR_path ="/sample_data/TWAS/DPR_weights.txt"


genotype = pd.read_csv(genotype_path,header=None, sep="\t")
genotype_mat = genotype.iloc[:,4:].T

phenotype = pd.read_csv(phenotype_path, sep="\t")
weights = pd.read_csv(DPR_path, sep="\t")
weight_mat = weights['ES_sum']
pheno= phenotype['phenotype']
cov = phenotype[['PC1','PC2','SEX','AGE']]

###Continous phenotype 
pval_C = SKAT.SKAT(genotype_mat,pheno, cov, weight_mat, 'C')

###Dichotomous phenotype 
phenotype_D = pd.read_csv(phenotype_path_D, sep="\t")
pheno_D = phenotype_D['phenotype']
cov_D = phenotype_D[['PC1','PC2','SEX','AGE']]

pval_D= SKAT.SKAT(genotype_mat, pheno_D, cov_D, weight_mat, 'D')

