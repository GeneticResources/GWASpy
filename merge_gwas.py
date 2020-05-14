# -*- coding: utf-8 -*-
"""

merge gwas summary statistics

Author: Xikun Han
Last update: 
    23 Sep 2019: first version
    26 Sep 2019: update pattern parameter

TODO:
    1. currently only support bolt-lmm results, consider plink (different versions), saige and other format in the future.
    


# parameters
dir_path='/working/lab_stuartma/xikunH/UKB/AMD/out/FVC_test'
trait='FVC'
time='20190101'
only_one_allele=True # only keep SNP with one Allele, not insersion or delesion
pattern='revised_bolt_imputed_ukb_imp_chr*_v3_*.bgen.assoc'

# test example

python ~/script/pygwas/merge_gwas.py --dir_path /working/lab_stuartma/xikunH/UKB/AMD/out/FVC_test --trait FVC  --keep_one_allele True


"""



import argparse
import os
import glob
import subprocess
from utils import str2bool # for bool format


def merge_gwas(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.trait is None:
        raise ValueError('The --trait flag is required.')
    if args.keep_one_allele is None:
        raise ValueError('The --keep_one_allele flag is required.')
        
    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as :'+ '20190101')
    
    
    dir_path = args.dir_path
    trait = args.trait
    time = args.time
    keep_one_allele = args.keep_one_allele
    pattern = args.pattern
    
    # get the file names
    files = [os.path.basename(f) for f in glob.glob(dir_path + "/" + pattern , recursive=False)]
    
    files.sort()
    
    if len(files) != 24:
        print('Warning: only have {} files.'.format(len(files)))
    else:
        print('Have {} files.'.format(len(files)))
    
    
    
    for i, file in enumerate(files):
        print(i+1, " : ",file)
    
    print("\n")
    # create gwas summary file, merge each filtered file
    gwas_file = trait + '_GWAS_result_all_' + time
    
    
    if software == 'bolt_lmm':
        # bashcmd = "awk '{if (FNR ==1 ) {$1=$1; print $0 }}' " + dir_path + "/" + files[0] +  " > " + dir_path + "/" + gwas_file
        # double braces as one in format
        bashcmd = "awk '{{if (FNR ==1 ) {{$1=$1; print $0 }} }}' {p1}/{p2}  >  {p3}/{p4}".format(p1= dir_path, p2=files[0], p3=dir_path, p4=gwas_file)
        
        print(bashcmd)
        subprocess.run(bashcmd, shell =True)
        
        for file in files:
            if keep_one_allele:
                bashcmd = "awk '{if ($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.3 && length($5) ==1 && length($6) ==1 ) {$12=sqrt(($11)^2)/sqrt($15); print $0 }}' " + dir_path + "/" + file +  " >> " + dir_path + "/" + gwas_file
            else:
                bashcmd = "awk '{if ($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.3) {$12=sqrt(($11)^2)/sqrt($15); print $0 }}' " + dir_path + "/" + file +  " >> " + dir_path + "/" + gwas_file
            print(bashcmd)
            subprocess.run(bashcmd, shell =True)
    elif software != 'bolt_lmm':
        print('Current only support bolt-lmm. need to do software.')



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Merge GWAS summary statistics.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--trait', type=str, required=True, help='trait name')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    parser.add_argument('--keep_one_allele', type=str2bool, required=True, help='if Ture, only keep SNP with one Allele, not insersion or delesion.')
    parser.add_argument('--pattern', type=str, default="revised_bolt_imputed_ukb_imp_chr*_v3_*.bgen.assoc", help='association pattern')
    parser.add_argument('--software', type=str, default="bolt_lmm", help='merge summary statistics from which software? bolt_lmm, plink, saige.')
    
    args = parser.parse_args()
    
    merge_gwas(args)
    
    
   














#####
