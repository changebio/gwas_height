import argparse

from PlinkCmd import parse_genotype
from LiftoverCmd import liftover_plink
from snp2h import height_prediction


def height_pipeline(g,meta=None,chain='hg38ToHg19',ms=['linear_regression'],snp='top200_all',path='full_model'):
    plink_path = parse_genotype(g)
    print('plink_path',plink_path)
    if chain is None:
        lift_file =None
    else:
        lift_file = liftover_plink(plink_path,chain)
        print('liftover_plink',lift_file)
    print('height prediction')
    height_prediction(plink_path,lift_file,ms=ms,snp=snp,path=path)

def init_args(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype', help='Genotype input file, the prefix of plink files (.bim,.fam,.bed)')
    parser.add_argument('-s', '--specimen', help='the sample information')
    #parser.add_argument('-l', 'liftover', help='the output file with suffix .Hg19 or .Hg38 from LiftoverCmd.py',default=None)
    parser.add_argument('-c', '--chain', help='Liftover Chain file',choices=['hg38ToHg19','hg19ToHg38',None],default=None)
    parser.add_argument('-m', '--model', help='a list of prediction models',choices=['linear_regression','logistic_avg','logistic_25perc','logistic_15perc','rfc_avg','rfc_25perc','rfc_15perc'],default=['linear_regression'],nargs='+')
    parser.add_argument('--snp-list', help='Choose type of snp list', choices=['top200','top200_all','top300','top300_all','top500','top500_all'], default='top200_all')
    parser.add_argument('--model-path',help='Choose the way of model generation',choices=['full_model','rep10cv10_model'], default= 'full_model')

    args = parser.parse_args(arguments)
    return args

def main(args=None):
    args = init_args(args)
    height_pipeline(args.genotype,args.specimen,args.chain,ms=args.model,snp=args.snp_list,path=args.model_path)
        
if __name__ == '__main__':
    main()