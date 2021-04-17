#########################################################################################
import pandas as pd

#########################################################################################
###for more details, please see our TIGAR tool https://github.com/yanglab-emory/TIGAR/
#########################################################################################

#########################################################################################
train_dosage_path = "/sample_data/simulation/train_dosage.txt"
dpr_snp_info_path = "/sample_data/simulation/snp_annot.txt"
dpr_weight_format_path = "/sample_data/simulation/bimbam.txt"

train_dosage = pd.read_csv(train_dosage_path, header = None, sep ="\t")

##prepare train genotype data into DPR format
##SNP anno info: rsID, POS, CHR (can't change the order)
snp_anno = train_dosage.loc[:,[1, 2, 0]]

snp_anno.to_csv(
        dpr_snp_info_path,
        header=False,
        index=None,
        sep='\t',
        mode='w',
        float_format='%f')

##dosage: rsID, POS, CHR (can't change the order) same in SNP anno
bimbam =  pd.concat([train_dosage.loc[:,[1, 2, 0]], train_dosage.loc[:,5:]], axis=1)

bimbam.to_csv(
        dpr_weight_format_path,
        header=False,
        index=None,
        sep='\t',
        mode='w',
        float_format='%f')