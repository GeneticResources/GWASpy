# -*- coding: utf-8 -*-
"""

Last update: 27 Sep 2019
Author: Xikun Han



find all h2 and ldsc log files, and extract results.



# parameters

dir_path='/working/lab_stuartma/xikunH/UKB/AMD/result/post_gwas/LDSCORE'
suffix='merge_LDSC_summary'
time='20190101'




# test example
python ~/script/pygwas/format_LDSC.py --dir_path /working/lab_stuartma/xikunH/UKB/AMD/result/post_gwas/LDSCORE --suffix merge_AMD_LDSC_summary --time 20190927 


"""





import argparse
import glob
import os 
import pandas as pd
import re


def format_LDSC(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
    if args.suffix is None:
        raise ValueError('The --suffix flag is required.')
    
    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as :'+ '20190101')


    dir_path = args.dir_path
    suffix = args.suffix
    time = args.time 
    
    if args.dir_out is not None:
        dir_out = args.dir_out
    else:
        dir_out = dir_path
    
    # read heritability logs
    file_h2 = glob.glob(dir_path + '/*_h2.log', recursive = False)
    
    file_name = []
    Total_scale_h2 = []
    Lambda_GC = []
    Mean_Chisquare = []
    Intercept = []
    Ratio = []
    
    print('\nStart h2 result\n')
    v_i = 1
    for file in file_h2:
        print(v_i, ": ", os.path.basename(file))
        v_i +=1
        with open (file, 'r') as fh:
            fhlist = fh.readlines()
            file_name.extend([os.path.basename(file).strip()])
            Total_scale_h2.extend([word.split(':')[1].strip() for word in fhlist if "scale h2" in word])
            Lambda_GC.extend([word.split(':')[1].strip() for word in fhlist if "Lambda GC" in word])
            Mean_Chisquare.extend([word.split(':')[1].strip() for word in fhlist if "Mean Chi" in word])
            Intercept.extend([word.split(':')[1].strip() for word in fhlist if "Intercept" in word])
            Ratio.extend([ word.split(':')[1].strip() if ":" in word else word.strip() for word in fhlist if "Ratio" in word ]) # not : some are < 1
    
    # print(Ratio)
    df_h2 = pd.DataFrame({'file_name':file_name ,'Total_scale_h2':Total_scale_h2, 'Lambda_GC':Lambda_GC, 'Mean_Chisquare':Mean_Chisquare, 'Intercept':Intercept,'Ratio':Ratio })
    df_h2.to_csv(dir_out + '/'+ suffix + '_h2_'+ time + '.txt', sep = '\t', index = False)
    
    
    print('\nStart LDSC result\n')
    # genetic correlation logs
    file_rg = glob.glob(dir_path + '/*_ldsc.log', recursive = False)

    df_rg = pd.Series(['p1', 'p2','rg', 'se', 'z', 'p', 'h2_obs', 'h2_obs_se', 'h2_int', 'h2_int_se', 'gcov_int', 'gcov_int_se'])
    v_i = 1
    for file in file_rg:
        print(v_i, ": ", os.path.basename(file))
        v_i +=1
        with open (file, 'r') as fh:
            fhlist = fh.readlines()
            result = [word for word in fhlist if re.search('sumstats.gz.+sumstats.gz', word)]
            result = [os.path.basename(word.strip()) for word in result[len(result)-1].split(' ') if word != '']
            # print(result)
            df_rg = pd.concat([df_rg, pd.Series(result)], axis = 1)
    
    df_rg = df_rg.transpose()
    df_rg.to_csv(dir_out + '/'+ suffix + '_rg_'+ time + '.txt', sep = '\t', index = False, header = False)
        
   

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Clump GWAS summary statistics.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for GWAS summary statistics.')
    parser.add_argument('--suffix', type=str, default='merge_LDSC_summary', help='file name suffix')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    parser.add_argument('--dir_out', type=str, default=None, help='the folder for out files, the default is --dir_path')
    
    args = parser.parse_args()
    
    format_LDSC(args)






















#####