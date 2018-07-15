import os.path
import pandas as pd
import argparse

from subprocess import Popen, PIPE
from pandas_plink import read_plink


def liftover_plink(gpath,chain):
    """
    input
        - gpath: plink file(bim,fam,G) to make a bed file. score is the index of snp in bim
        - chain: hg38ToHg19 or hg19ToHg38
    out_put: a liftover file by chain, such as 'hg38ToHg19.over.chain.gz'
    """
    tt_plink = read_plink(gpath)
    tt_bim = tt_plink[0]
    bed = bim2bed(tt_bim)
    out_put = gpath
    out_suffix = chain.split('To')
    bed.to_csv('.'.join([out_put,out_suffix[0]]),sep='\t',index=False,header=False)
    
    chain_path='lib/'+chain+'.over.chain.gz'
    if not os.path.isfile(chain_path):
        print('chain file do not exist.')
    cmd = ['lib/liftOver','.'.join([out_put,out_suffix[0]]),chain_path, '.'.join([out_put,out_suffix[1]]), out_put+'.unmaped']
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stderr)
    return '.'.join([out_put,out_suffix[1]])

def bim2bed(bim):
    """
    sex chromosome from number(23,24,25) to X and Y
    chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    The 9 additional optional BED fields are:

    name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
    shade
    strand - Defines the strand. Either "." (=no strand) or "+" or "-".

    """
    chrom = bim.chrom 
    chrom = 'chr' + chrom.astype('str')
    chrom = chrom.replace('chr23', 'chrX')
    bed = pd.DataFrame({'chr':chrom,'start':bim.pos-1,'end':bim.pos,'name':bim.snp,'score':bim.i,'strand':'.'},columns=['chr','start','end','name','score','strand'])
    #bed = pd.merge(bed,bim.iloc[:,[0,1,3,4,5,6]], left_index=True, right_index=True)
    return bed

def init_args(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype', help='Genotype input file, the prefix of plink files (.bim,.fam,.bed)')
    parser.add_argument('-c', '--chain', help='Liftover Chain file',choices=['hg38ToHg19','hg19ToHg38'],default='hg38ToHg19')
    args = parser.parse_args(arguments)
    return args

def main(args=None):
    args = init_args(args)
    liftover_plink(args.genotype,args.chain)

        
if __name__ == '__main__':
    main()