# -*- coding: utf-8 -*-
"""


sum PRS scores from plink

Author: Xikun Han
Last update: 
    03 Oct 2019: first version
    


# parameters

dir_path = '/working/lab_stuartma/xikunH/UKB/AMD/PRS/prediction/mtag_AMD2015dbGap_kaiser_2012GWAS_3trait_UKBB_BMES/PRS_plink'
suffix = 'test'




# test example

# BMES < 1 minutes
python ~/script/pygwas/PRS_sum_score.py --dir_path '/working/lab_stuartma/xikunH/UKB/AMD/PRS/prediction/mtag_AMD2015dbGap_kaiser_2012GWAS_3trait_UKBB_BMES/PRS_plink'  --suffix test

# UKB 791 11 ~ 1.5 h
python ~/script/pygwas/PRS_sum_score.py --dir_path '/working/lab_stuartma/xikunH/UKB/AMD/PRS/prediction/mtag_AMD2015GWAS_kaiser_2012GWAS_3trait_UKBB_UKB/PRS_plink'  --suffix test_UKB



"""



import argparse
import os
import glob
import pandas as pd
import numpy as np
    




def PRS_sum_score(args):
    if args.dir_path is None:
        raise ValueError('The --dir_path flag is required.')
        
    if args.suffix is None:
        raise ValueError('The --suffix flag is required.')
        
    if args.dir_out is None:
        args.dir_out = os.path.dirname(args.dir_path)
        print('The --dir_out flag is set as: '+ args.dir_out)
        
    if args.time is None:
        args.time='20190101'
        print('The --time flag is set as : '+ '20190101')
    
    dir_path = args.dir_path
    dir_out = args.dir_out
    suffix = args.suffix
    time = args.time
    
    pattern='.profile'
    files = [os.path.basename(f) for f in glob.glob(dir_path + "/*" + pattern , recursive=False)]
    print('\n{} files in {}'.format(len(files), dir_path))
    
    def f_sum_files(files):
        for i, file in enumerate(files):
            # print(i, file)
            if(i == 0):
                df = pd.read_csv(dir_path + '/' + file, delim_whitespace = True, usecols = ['FID', 'SCORESUM'])
            else:
                df_one = pd.read_csv(dir_path + '/' + file, delim_whitespace = True, usecols = ['FID', 'SCORESUM'])
                df['SCORESUM'] = df['SCORESUM'] + df_one['SCORESUM']
        return df
    
    
    for i in range(11):
        k = i +1
        print('\n')
        print(k)
        files = [os.path.basename(f) for f in glob.glob(dir_path + '/*.S'+ str(k)  + pattern , recursive=False)]
        print('{} files for S{}.'.format(len(files), k))
        if k ==1:
            df = f_sum_files(files)
            df.rename(columns = {'SCORESUM':'SCORESUM_S'+str(k)}, inplace = True)
        else:
            df_one = f_sum_files(files)
            df_one.rename(columns = {'SCORESUM':'SCORESUM_S'+str(k)}, inplace = True)
            df = df.merge(df_one, on = 'FID', copy = False)
    
    df.to_csv(dir_out + '/'+ suffix + '_'+ time + '.txt', sep = '\t', index = False)



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='sum PRS scores from plink.')
    parser.add_argument('--dir_path', type=str, required=True, help='path to the folder for plink sum score files.')
    parser.add_argument('--suffix', type=str, required=True, help='output file name (eg. out_df_PRS_sum_score)')
    parser.add_argument('--time', type=str, default='20190101', help='time as suffix')
    parser.add_argument('--dir_out', type=str, help='output path.')

    args = parser.parse_args()
    
    PRS_sum_score(args)
    













#####
