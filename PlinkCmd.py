import os
import pandas as pd
import argparse

from subprocess import Popen, PIPE
from os import listdir, getcwd

def path_split(s):
    ss = s.split('/')
    file = ss[-1]
    if len(ss)>1:
        base = '/'.join(ss[:-1])
    else:
        base = getcwd()
    return base,file
        

def parse_genotype(g):
    #list file
    gdir,gfile = path_split(g)
    gfiles = [f for f in listdir(gdir) if gfile in f]
    print(gfiles)
    #match file type (bim fam G, vcf, ped map)
    if len(gfiles)==1:
        if '.vcf' in gfile:
            plink_path=vcf2plink(g)
        elif '.bim' in gfile or '.fam' in gfile:
            plink_path=g[:-4]
        elif '.ped' in gfile or '.map' in gfile:
            plink_path=ped2plink(g)
        else:
            print('error')
    else:
        if '.'.join([gfile,'bim']) in gfiles:
            plink_path=g
        elif '.'.join([gfile,'ped']) in gflies:
            plink_path=ped2plink(g)
        else:
            print('input error')
    #gt = read_plink(plink_path)
    return plink_path

    
def vcf2plink(vcf_file):
    out_put=vcf_file.split('.vcf')[0]
    cmd = ['plink','--vcf',vcf_file,'--biallelic-only','strict','--make-bed','--impute-sex','ycount','--out',out_put]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stderr)
    return out_put

def ped2plink(ped_file):
    out_put = ped_file
    cmd = ['plink','--file',ped_file,'--biallelic-only','strict','--make-bed','--impute-sex','ycount','--out',out_put]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stderr)
    return out_put

def guess_gender(gpath):
    out_put = gpath
    cmd = ['plink','--bfile',gpath,'--impute-sex','ycount','--make-bed','--out',out_put]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stderr)
    return out_put
   
def read_gender(g):
    gprefix = g[:-9]
    if not os.path.isfile(g):
        print('.sexcheck file do not exist. running plink to guess gender')
        out_put = guess_gender(gprefix)
    tt_gender = pd.read_csv(g,delim_whitespace=True)
    tt_gender['gender']=(tt_gender.F<0.5)+1
    tt_gender['nationality']=0
    tt_gender.index = list(tt_gender.IID)
    return tt_gender

def init_args(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype', help='Genotype input file')
    args = parser.parse_args(arguments)
    return args

def main(args=None):
    args = init_args(args)
    parse_genotype(args.genotype)
    
        
if __name__ == '__main__':
    main()