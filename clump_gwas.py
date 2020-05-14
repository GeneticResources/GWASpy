# -*- coding: utf-8 -*-
"""

Last update: 25 Sep 2019
Author: Xikun Han



A pipeline to clump gwas.

Support all kinds of gwas summary statistics with SNP and P value columns.





# parameters

dir_path='/working/lab_stuartma/xikunH/UKB/AMD/out/S3'
gwas_summary_file='S3_GWAS_result_all_20181028'
trait='S3'
time='20190101'

LDref_bfile='/reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/LD_ref_201803_noLongSNP'

#PLINK Clumping thresholds
kbd='1000'
r2level='0.01'
p1thres='5e-8'
p2thres='5e-8'
clump_snp='SNP'
clump_field='P_BOLT_LMM'
# plink options
threads='1'
memory='4000mb'


# test example
python clump_gwas.py --dir_path /working/lab_stuartma/xikunH/UKB/AMD/out/S3 --gwas_summary_file S3_GWAS_result_all_20181028 --trait S3 --clump_snp SNP --clump_field P_BOLT_LMM


"""



import argparse
import os
import subprocess


def clump_gwas(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.gwas_summary_file is None:
        raise ValueError('The --gwas_summary_file flag is required.')
    if args.trait is None:
        raise ValueError('The --trait flag is required.')
    
    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as :'+ '20190101')
    
    
    dir_path = args.dir_path
    gwas_summary_file = args.gwas_summary_file
    
    
    gwas_summary_file_full_path = dir_path + '/'+ gwas_summary_file if os.path.basename(gwas_summary_file) == gwas_summary_file else gwas_summary_file
    out_file = args.dir_path +'/'+ args.trait + '_GWAS_result_all_'+ args.time
    
    # plink clump 
    bashcmd = '/software/plink/plink-1.90b4.1/plink --bfile ' + args.LDref_bfile + ' --clump ' + gwas_summary_file_full_path + ' --clump-snp-field ' + args.clump_snp +  ' --clump-field ' + args.clump_field + ' --clump-p1 ' + args.p1thres + ' --clump-p2 ' + args.p2thres + ' --clump-kb ' + args.kbd + ' --clump-r2 ' + args.r2level + ' --out  ' + out_file + ' --threads ' + args.threads + ' --memory ' + args.memory
    print(bashcmd + '\n')
    subprocess.run(bashcmd, shell =True)
    
    clump_OK = os.path.exists(out_file+'.clumped')
    clumped_file=out_file+'.clumped'
    
    # which column is SNP
    with open(gwas_summary_file_full_path, 'r') as f:
        f_header = f.readline().split()
    
    snp_index = f_header.index(args.clump_snp) + 1
    print('Column ' + str(snp_index) + ' is ' + args.clump_snp)
    
    
    if not clump_OK:
        print('Error in clump or no significant SNP.')
    else:
        print('Clump successed.')
        independent_file=args.dir_path +'/indep_GWhits_' + args.trait + '_' + args.kbd + 'kb.r2_' + args.r2level + '_p1_'+ args.p1thres + '_p2_' + args.p1thres + '_' + args.time
        bashcmd = "awk  'NR==FNR{a[$3]; next;}{if(FNR == 1 || $" + str(snp_index) + " in a) print $0}' "+ clumped_file + " "+ gwas_summary_file_full_path + " > " + independent_file
        print('\nStart select independent SNPs using command:\n')
        print(bashcmd)
        subprocess.run(bashcmd, shell =True)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Clump GWAS summary statistics.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--gwas_summary_file', type=str, required=True, help='file name of GWAS summary statistics.')
    parser.add_argument('--trait', type=str, required=True, help='trait name')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    parser.add_argument('--clump_snp', type=str, default='SNP', help='snp column name in GWAS summary statistics')
    parser.add_argument('--clump_field', type=str, default='P', help='P value column name in GWAS summary statistics')
    parser.add_argument('--LDref_bfile', type=str, default='/reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/LD_ref_201803_noLongSNP', help='LD referebce file')
    parser.add_argument('--kbd', type=str, default='1000', help='clump parameter, kb')
    parser.add_argument('--r2level', type=str, default='0.01', help='clump parameter, r2')
    parser.add_argument('--p1thres', type=str, default='5e-8', help='clump parameter, p1 thres')
    parser.add_argument('--p2thres', type=str, default='5e-8', help='clump parameter, p2 thres')
    parser.add_argument('--threads', type=str, default='1', help='threads in plink (v1.9)')
    parser.add_argument('--memory', type=str, default='4000mb', help='memory in plkink (v1.9)')
    
    args = parser.parse_args()
    
    clump_gwas(args)
    
    
   
























#####