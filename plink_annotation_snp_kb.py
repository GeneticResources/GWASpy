# -*- coding: utf-8 -*-
"""

Annoatate SNPs using plink.

Last update: 23 Sep 2019
Author; Xikun Han



# parameters

dir_path='/working/lab_stuartma/xikunH/UKB/AMD/out/S3'
trait='S3'
snp_file='indep_GWhits_S3_1000kb.r2_0.01_p1_5e-8_p2_5e-8_20190101'

snp_list = dir_path + '/' + snp_file

dir_plink_annotation='/working/lab_stuartma/xikunH/UKB/OCT/result/old_imputation/post_gwas/plink_annotation/ALLsnp_new'
snp_index = '1'


# test example

python plink_annotation_snp_kb.py --dir_path /working/lab_stuartma/xikunH/UKB/AMD/out/FVC --snp_file indep_GWhits_FVC.filtered.1000kb.r20.01.boltlmm.set_20181028 --trait FVC



"""



import argparse
import os
import subprocess
import pandas as pd
from functools import reduce


def plink_annotation_snp_kb(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.snp_file is None:
        raise ValueError('The --snp_file flag is required.')
    if args.trait is None:
        raise ValueError('The --trait flag is required.')

    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as :'+ '20190101')
    
    # snp_list = args.dir_path + '/' + args.snp_file
    
    snp_list = args.dir_path + '/' + args.snp_file if os.path.basename(args.snp_file) == args.snp_file else args.snp_file
    
    snp_index = args.snp_index
    dir_path = args.dir_path
    
    dir_plink_annotation='/working/lab_stuartma/xikunH/UKB/OCT/result/old_imputation/post_gwas/plink_annotation/ALLsnp_new'
    
    # annotation 0
    annotation=dir_plink_annotation + '/ALLsnp_new_0.annot'
    bashcmd = "awk 'NR==FNR {a[$" + snp_index + "]; next;} {if(($" + snp_index + " in a) || (FNR==1)) print $0}' " + snp_list+ " "+  annotation + "   > " + dir_path + "/tmp_kb_0 "
    subprocess.run(bashcmd, shell =True)


    # ALLsnp_200.annot
    annotation=dir_plink_annotation + '/ALLsnp_new_200.annot'
    bashcmd = "awk 'NR==FNR {a[$" + snp_index + "]; next;} {if(($" + snp_index + " in a) || (FNR==1)) print $1,$6}' " + snp_list+ " "+  annotation + "   > " + dir_path + "/tmp_kb_200 "
    subprocess.run(bashcmd, shell =True)
    
    
    # ALLsnp_400.annot
    annotation=dir_plink_annotation + '/ALLsnp_new_400.annot'
    bashcmd = "awk 'NR==FNR {a[$" + snp_index + "]; next;} {if(($" + snp_index + " in a) || (FNR==1)) print $1,$6}' " + snp_list+ " "+  annotation + "   > " + dir_path + "/tmp_kb_400 "
    subprocess.run(bashcmd, shell =True)
    
    
    # ALLsnp_1000.annot
    annotation=dir_plink_annotation + '/ALLsnp_new_1000.annot'
    bashcmd = "awk 'NR==FNR {a[$" + snp_index + "]; next;} {if(($" + snp_index + " in a) || (FNR==1)) print $1,$6}' " + snp_list+ " "+  annotation + "   > " + dir_path + "/tmp_kb_1000 "
    subprocess.run(bashcmd, shell =True)
    
    
    
    df_P = pd.read_csv(snp_list, delim_whitespace = True)
    
    if not df_P.columns[int(snp_index) - 1] == 'SNP':
        print('Warning: the ' + snp_index + ' column is not SNP. Will be changed to SNP.')
        df_P.rename(columns={df_P.columns[int(snp_index) - 1]: "SNP"}, inplace = True)
    
    df_0 = pd.read_csv(dir_path + '/tmp_kb_0',delim_whitespace = True)
    df_0.rename(columns={'ANNOT': "withingen"}, inplace = True)
    
    df_200 = pd.read_csv(dir_path + '/tmp_kb_200',delim_whitespace = True)
    df_200.rename(columns={'ANNOT': "within200kb"}, inplace = True)
    
    df_400 = pd.read_csv(dir_path + '/tmp_kb_400',delim_whitespace = True)
    df_400.rename(columns={'ANNOT': "withi400kb"}, inplace = True)
    
    df_1000 = pd.read_csv(dir_path + '/tmp_kb_1000',delim_whitespace = True)
    df_1000.rename(columns={'ANNOT': "within1000kb"}, inplace = True)
    
    # remove tmp files
    os.remove(dir_path + '/tmp_kb_0')
    os.remove(dir_path + '/tmp_kb_200')
    os.remove(dir_path + '/tmp_kb_400')
    os.remove(dir_path + '/tmp_kb_1000')
    
    df_merge = reduce(lambda  left,right: pd.merge(left,right,on=['SNP'], how='outer'), [df_P,df_0,df_200,df_400,df_1000])

    df_merge.rename(columns={'CHR_x':'CHR', 'BP_x':'BP'}, inplace=True)
    df_merge.sort_values(by = ['CHR', 'BP'], ascending=[1, 1], inplace = True)
    df_merge.to_csv(dir_path + '/'+ args.trait + '_plink_annotation_out_'+ args.time + '.txt', sep = '\t', index = False)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Annotate SNPs with genes using plink annotation.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--snp_file', type=str, required=True, help='file name of SNP list having `SNP` column.')
    parser.add_argument('--snp_index', type=str, default='1', help='snp column name in snp file, default is the first column')
    parser.add_argument('--trait', type=str, required=True, help='trait name')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    
    args = parser.parse_args()
    
    plink_annotation_snp_kb(args)
    
    
   


























#####