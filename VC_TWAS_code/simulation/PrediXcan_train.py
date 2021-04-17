#########################################################################################
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import ElasticNetCV
import statsmodels.api as sm

#########################################################################################
###for more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR/
#########################################################################################

#########################################################################################
###Load train dosage and true gene expression, the order of sample ID are the same in these files
###First 5 rows is chrom, rsID, POS, REF, ALT
train_dosage_path = "/sample_data/simulation/train_dosage.txt"
train_expression_path = "/sample_data/simulation/expression_train.txt"
PrediXcan_weight_path = "/sample_data/simulation/prediXcan_weights.txt"

train_dosage = pd.read_csv(train_dosage_path, header = None,sep ="\t")

train_exp = pd.read_csv(train_expression_path, header = None,sep ="\t")

train_dosage_mat= train_dosage.iloc[:,5:].T

# k : k-fold cross validation,default 5
#Alpha : ratio for L1 and L2 penalty in elastic net regression,default 0.5
#         Alpha=0: Lasso Regression
#         Alpha=1: Ridge Regression
#         0 < Alpha < 1: Elastic Net Regression

def elastic_net(train_dosage, train_exp, test_dosage = None, test_exp = None, k = 5, Alpha = 0.5):

    if test_dosage is None:
        test_dosage = train_dosage
        test_exp = train_exp

    reg = ElasticNetCV(
        l1_ratio = Alpha,
        fit_intercept = False,
        alphas = np.arange(0,1.01,0.01),
        selection = 'random',
        cv = k).fit(train_dosage, train_exp)

    Lambda = reg.alpha_
    cvm = np.min(reg.mse_path_)
    beta = reg.coef_

    predY = reg.predict(test_dosage)

    lm = sm.OLS(test_exp, sm.add_constant(predY)).fit()

    Rsquared = lm.rsquared

    Pvalue = lm.f_pvalue

    return beta, Rsquared, Pvalue, Alpha, Lambda, cvm

##get PrediXcan weights
weights = train_dosage.iloc[:,0:5]
weights.columns = ['CHROM', 'rsID', 'POS', 'REF', 'ALT']
weights['ES'], R2, Pvalue, Alpha, Lambda, cvm = elastic_net(train_dosage_mat, train_exp)
weights = weights[weights['ES']!=0]
weights.to_csv(
    PrediXcan_weight_path, header=False, index=None, sep='\t', mode='w')
