
###################################################################
## for more details, please see https://github.com/yanglab-emory/TIGAR
###################################################################
# Import packages needed
import numpy as np
import pandas as pd
import SKAT
###############################################################
#load genotype data (assume the order of SNP in genotype and weight are the same, the order of sample ID in genotype and phenotype are the same)
dosage = pd.read_csv(geno_path,sep='\t',header=None)
#dosage input should be numpy format
dosage_mat =dosage.values
#load weight
weight_mat = pd.read_csv(weight_path,sep='\t')
#dosage input should be numpy format
weight =weight.values

#load phenotype (phenotype should be pd format,phenotype_name is name of phenotype, cov_name is name of cov)
phenotype = pd.read_csv(pheno_path,sep='\t')
phenotype_y=phenotype[['phenotype_name']]
phenotype_cov=phenotype[['cov_name']]


p_value=SKAT.SKAT(dosage_mat,phenotype_y,phenotype_cov,weight_mat,"C") #phenotype is continuous, set "C"
p_value=SKAT.SKAT(dosage_mat,phenotype_y,phenotype_cov,weight_mat,"D") #phenotype is dichotomous, set "D"


