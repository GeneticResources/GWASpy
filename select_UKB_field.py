# -*- coding: utf-8 -*-
"""

select UKBB fields 

Last update: 27 Sep 2019
Author; Xikun Han



# parameters


file_out = '~/test.txt'
fid = '4689,6119,131188'

# UKB phenotype data
# fh = '/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas/BOLT_LMM/bioMarker/data/20190418_phenotypes_test.txt'
fh = '/reference/data/UKBB_500k/20190805_update/ukb34945.tab'



# test example


python ~/script/pygwas/select_UKB_field.py --file_out '~/todelete/test' --input_field '4689,6119,131188' 


"""



import argparse
import subprocess


def select_UKB_field(args):
    '''
    fo: output path for the selected IDs
    fid: a list of id to select
    fh: all UKBB phenotype file
    '''
    
    if args.file_out is None:
        raise ValueError('The --file_out flag is required.')
    if args.input_field is None:
        raise ValueError('The --input_field flag is required.')
    if args.UKB_phenotype_file is None:
        args.UKB_phenotype_file='/reference/data/UKBB_500k/20190805_update/ukb34945.tab'
        print('The --UKB_phenotype_file flag is set as: '+ args.UKB_phenotype_file)
    
    fo = args.file_out
    fid = args.input_field
    fh = args.UKB_phenotype_file
    
    # read the first line from phenotype file
    with open(fh,'r') as f:
        first_line = f.readline()
    col = [x for x in first_line.split("\t")]
    
    # create the field ID list
    fid = [index.strip() for index in fid.split(',')]
    fid = ['f.'+ x +'.' for x in fid]
    fid.append('f.eid')
    print('\nThe fields to select: ')
    print(fid)
    
    # function to test contain
    def f_contain(x, list_):
        for i in list_:
            if(i in x):
                return True
        return False
    
    
    # create col index to select
    col_index=[]
    index = 0
    for f in col:
        index = index + 1
        if(f_contain(f, fid)):
            col_index.append(index)
    
    col_index_str = ",".join( str(item) for item in col_index) 
    
    bashcmd = 'cut -f '+ col_index_str + ' ' + fh + ' >  ' + fo
    print('\nRun CMD: \n' + bashcmd)
    subprocess.run(bashcmd, shell =True)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Select UKB fields from big phenotype file.')
    parser.add_argument('--file_out', type=str, required=True, help='path to the file for selected fileds.')
    parser.add_argument('--input_field', type=str, required=True, help='fields ids to select (split by comma), eg. 4689,6119,131188')
    parser.add_argument('--UKB_phenotype_file', type=str, help='path of UKBB phenotype file, default file is the latest one.')
    
    args = parser.parse_args()
    
    select_UKB_field(args)


    








#####