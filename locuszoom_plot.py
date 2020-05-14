# -*- coding: utf-8 -*-
"""

Annoatate SNPs using plink.

Last update: 23 Sep 2019
Author; Xikun Han



# parameters


dir_path='/working/lab_stuartma/xikunH/UKB/AMD/out/FVC_test'
trait='FVC'
time='20190101'

snp_file='FVC_plink_annotation_out_20190101_locuszoom_test.txt'
snp = 'SNP'
chromo = 'CHR'
pos = 'BP'

metal_file = 'FVC_GWAS_result_all_20190101'
delim = 'space'
pvalcol = 'P_BOLT_LMM'
markercol = 'SNP'

qsub = False



# test example


python ~/script/pygwas/locuszoom_plot.py --dir_path '/working/lab_stuartma/xikunH/UKB/AMD/out/FVC_test' --trait 'FVC' --time '20190101' --snp_file='FVC_plink_annotation_out_20190101_locuszoom_test.txt' --snp 'SNP' --chromo 'CHR' --pos 'BP' --metal_file 'FVC_GWAS_result_all_20190101' --delim 'space' --pvalcol 'P_BOLT_LMM' --markercol 'SNP' --qsub False


"""



import argparse
import os
import pandas as pd
import numpy as np
import subprocess
import itertools
import glob
from utils import str2bool
import time



def locuszoom_plot(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.trait is None:
        raise ValueError('The --trait flag is required.')
    if args.snp_file is None:
        raise ValueError('The --snp_file flag is required.')
    if args.metal_file is None:
        raise ValueError('The --metal_file flag is required.')
    if args.qsub is None:
        raise ValueError('The --qsub flag is required.')

    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as :'+ '20190101')
    
    
    dir_path = args.dir_path
    trait = args.trait
    time = args.time
    snp_file = args.snp_file
    snp = args.snp
    chromo = args.chromo
    pos = args.pos
    
    metal_file = args.metal_file
    delim = args.delim
    pvalcol = args.pvalcol
    markercol = args.markercol
    qsub =  args.qsub
    

    snp_file = dir_path + '/' + snp_file if os.path.basename(snp_file) == snp_file else snp_file
    metal_file = dir_path + '/' + metal_file if os.path.basename(metal_file) == metal_file else metal_file
    
    locuszoom_folder = dir_path + '/locuszoom_' + trait + '_' + time 
    
    os.makedirs(locuszoom_folder,exist_ok=True)
    
    # df = pd.read_csv(snp_file, delim_whitespace = True)
    try:
        df = pd.read_csv(snp_file, delim_whitespace = True)
    except Exception as e:
        print(e)
        df = pd.read_csv(snp_file, sep = '\t')
    print(df.head)

    df.sort_values(by = [chromo, pos], ascending=[1, 1], inplace = True)
    snp_name = df[snp]
    print('Locuszoom for: \n ' + str(snp_name) + '\n')
    
    print('qsub: ', qsub)
    
    if qsub:
        print('\n*********************** qsub jobs ***********************\n')
        # another good ref: https://github.com/GeneticResources/h2-GRE/blob/master/gwas_sim/qsub.py
        
        locuszoom_software='/working/lab_stuartma/xikunH/software/locuszoom/locuszoom/bin/locuszoom'
        for SNP in snp_name:
            print(SNP)
            
            bashcmd = '#PBS -l ncpus=1 \n'
            bashcmd += '#PBS -l mem=100mb \n'
            bashcmd += '#PBS -l walltime=4:00:00 \n'
            bashcmd += '#PBS -e ' +  locuszoom_folder +'/'+ SNP.replace(':', '_') + '.e\n'
            bashcmd += '#PBS -o ' +  locuszoom_folder +'/'+ SNP.replace(':', '_') + '.o\n'
            
            bashcmd += 'module load locuszoom/1.4 \n'
            
            bashcmd += locuszoom_software + ' --plotonly --plotonly --flank 400kb --prefix EUR  --source 1000G_March2012 --build hg19 --pop EUR ' 
            bashcmd += ' --metal ' + metal_file +  ' --delim ' + delim + ' --pvalcol ' + pvalcol + ' --markercol ' + markercol + ' --refsnp ' + SNP
            bashcmd += " --snpset NULL geneFontSize=0.8 smallDot=0.3 largeDot=0.9 format=pdf  legend=auto metalRug="+ SNP + " showRecomb=TRUE ldCol=rsquare"
            bashcmd += ' --rundir ' + locuszoom_folder
            
            bashScript = locuszoom_folder  + '/single_locuszoom_' + SNP.replace(':', '_') + '.sh'
            with open(bashScript, "w") as fh:
                fh.write(bashcmd)
            subprocess.run('qsub ' + bashScript, shell =True)
            
    else:
        print('\n*********************** merge pdfs ***********************\n')
        df['locuszoom'] = 0
        for SNP in snp_name:
            locuszoom_pdf = glob.glob(locuszoom_folder + "/EUR_*" + SNP.replace(':', '_') + '.pdf')
            if len(locuszoom_pdf) != 1:
                print('Warning: find ' + str(len(locuszoom_pdf)) + ' pdf for SNP ' + SNP)
                continue
            df['locuszoom'] = np.where(df[snp] == SNP, 1, df['locuszoom'])
            locuszoom_pdf = locuszoom_pdf[0]
            bashcmd = 'pdfseparate -f 1 -l 1 ' + locuszoom_pdf + ' ' + locuszoom_pdf.replace('.pdf', '_1.pdf')
            bashcmd
            subprocess.run(bashcmd, shell =True)
            
        # create a locuszoom index file
        df_locuszoom = df[df.locuszoom == 1].copy()
        # merge pdf as one page per pdf
        locuszoom_merge_per_pdf = locuszoom_folder + '/locuszoom_merge_' + trait + '_' + time + '.pdf'
        locuszoom_merge_per_pdf_tmp = locuszoom_folder + '/locuszoom_merge_' + trait + '_' + time + '_tmp.pdf'
        
        if os.path.exists(locuszoom_merge_per_pdf):
            os.remove(locuszoom_merge_per_pdf)
        
        
        for SNP in snp_name:
            locuszoom_pdf = glob.glob(locuszoom_folder + "/EUR_*" + SNP.replace(':', '_') + '.pdf')
            if len(locuszoom_pdf) != 1:
                print('Warning: find ' + str(len(locuszoom_pdf)) + ' pdf for SNP ' + SNP + '\n')
                continue
            locuszoom_pdf = locuszoom_pdf[0]
            locuszoom_pdf_1 = locuszoom_pdf.replace('.pdf', '_1.pdf')
            bashcmd = 'pdfseparate -f 1 -l 1 ' + locuszoom_pdf + ' ' + locuszoom_pdf_1
            bashcmd
            subprocess.run(bashcmd, shell =True)
            
            if not os.path.exists(locuszoom_merge_per_pdf):
                bashcmd = 'cp ' + locuszoom_pdf_1 + ' ' + locuszoom_merge_per_pdf
                subprocess.run(bashcmd, shell =True)
            else:
                bashcmd = 'pdfunite ' + locuszoom_merge_per_pdf + ' ' + locuszoom_pdf_1 + ' ' + locuszoom_merge_per_pdf_tmp
                subprocess.run(bashcmd, shell =True)
                
                bashcmd = 'mv ' + locuszoom_merge_per_pdf_tmp   + ' ' +  locuszoom_merge_per_pdf
                subprocess.run(bashcmd, shell =True)
        
        # merge pdf as one page per pdf
        locuszoom_merge_9_pdf = locuszoom_folder + '/locuszoom_merge_3by3_' + trait + '_' + time + '.pdf'
        
        if os.path.exists(locuszoom_merge_9_pdf):
            os.remove(locuszoom_merge_9_pdf)
        
        
        pdfjam_software = '/working/lab_stuartma/xikunH/software/pdfjam/bin/pdfjam'
        
        def get_file_name_for_SNP(x):
            # x is the list of SNPs
            file_list = []
            for i in x:
                i_file = glob.glob(locuszoom_folder + "/EUR_*" + i.replace(':', '_') + '_1.pdf')
                i_file = i_file[0]
                file_list.append(i_file)
            return file_list
            
        file_list = get_file_name_for_SNP(df_locuszoom[snp])
        bashcmd = pdfjam_software + ' ' + ' '.join(file_list) + ' --nup 3x3 --landscape --outfile ' + locuszoom_merge_9_pdf
        subprocess.run(bashcmd, shell =True)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='LocusZoom plot.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--trait', type=str, required=True, help='trait name')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    
    parser.add_argument('--snp_file', type=str, required=True, help='file name of SNPs with snp name, chr, pos (for sorting).')
    parser.add_argument('--snp', type=str, default='SNP', help='snp name')
    parser.add_argument('--chromo', type=str, default='CHR', help='chromosome')
    parser.add_argument('--pos', type=str, default='BP', help='position')
    
    parser.add_argument('--metal_file', type=str, required=True, help='GWAS summary statistics as backgroud reference.')
    parser.add_argument('--delim', type=str, default='space', help='delim in metal_file, space, tab...')
    parser.add_argument('--pvalcol', type=str, default='P', help='p value cololumn in metal_file.')
    parser.add_argument('--markercol', type=str, default='SNP', help='SNP cololumn in metal_file.')
    
    parser.add_argument('--qsub', type=str2bool,required=True,default=False, help='True if qsub jobs, False to merge pdfs.')
        
    args = parser.parse_args()
    
    # if qsub is True, run merge pdfs after ~ 1 hours. 
    if args.qsub:
        locuszoom_plot(args)
        time.sleep(3600)
        args.qsub = False
        locuszoom_plot(args)
    else:
        locuszoom_plot(args)
        


    
    






#####
