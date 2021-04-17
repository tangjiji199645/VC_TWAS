#############################################################
# Import packages needed

import pandas as pd
import subprocess


#############################################################
DPR_path = "/sample_data/simulation/DPR"
bimbam_pth = "/sample_data/simulation/train_dosage.txt"
pheno_pth = "/sample_data/simulation/expression_train.txt"
snpannot_pth = "/sample_data/simulation/snp_annot.txt"


DPR_call_args = [DPR_path, 
            '-g', bimbam_pth, 
            '-p', pheno_pth, 
            '-a', snpannot_pth, 
            '-dpr', "1", 
            '-notsnp',
            '-o', 'DPR']

dpr_file_dir = "/home/stang/tangji/VC_TWAS/simulation/dpr_test/data"

subprocess.check_call(
            DPR_call_args,
            cwd=dpr_file_dir,
            stdout=subprocess.DEVNULL)

##output of DPR is under dpr_path/out 
##.param.txt is the estimated weight file