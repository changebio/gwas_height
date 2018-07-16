import numpy as np
import pandas as pd
import crowdai
import argparse

from sklearn import linear_model
from collections import OrderedDict as odict
from pandas_plink import read_plink

from LiftoverCmd import liftover_plink
from PlinkCmd import read_gender
from models import load_model
import warnings
warnings.filterwarnings("ignore")

def height_prediction(plink,lift,meta=None,ms=['linear_regression'],snp='top200_all',path='full_model'):
    """
    for example:
        height_prediction('ddat/crowdAI/temp','ddat/crowdAI/temp.Hg19')
    """
    #input files genotype,liftover,refbim, sample info
    refbim_path = '/'.join([path,snp,'ref.bim'])
    refbim = read_bim(refbim_path)
    if lift is not None:
        tt_lift = pd.read_csv(lift,delim_whitespace=True,header=None,names=['chr','start','end','name','score','strand'])
    else:
        tt_lift = None
        
    tt_plink = read_plink(plink)
    
    df=height_feature(refbim,tt_plink,tt_lift)
    tt_gender = read_gender(plink+'.sexcheck')
    if meta is not None:
        tt_meta = pd.read_csv(meta,index_col=0)
        if all(list(tt_gender.gender)==list(tt_meta.gender)):
            print('warning input genders do not match with Guess Biological Sex')
            print('Using input genders')
            tt_gender.gender = list(tt_meta.gender)
            if 'height' in tt_meta.columns:
                tt_gender['height']=list(tt_meta.height)

    #plt.hist(tt_gender.F)

    
    ph=height_predict(tt_gender,df,ms,snp,path)
    print(ph)
    ph.to_csv(plink+'.csv')
    return ph

def height_predict(meta,gt,ms=['linear_regression'],snp='top200_all',path='full_model'):
    """
    input
        meta: information of sample,such as gender
        gt: genotype data after height_feature, rows are samples, columns are snps.
        ml: list of prediction model (linear, logistic...)
        snp: the ways of snp selection (top200, top200_all, top300,top300_all,top500,...)
        path: model directory, default=full_model (full_model,10cv_model)
        -full_model/
            -gender_residual.model
            top200/
                -pca_95.model
                -linear_regression.model
                -logistic_avg
                -logistic_25perc
            top200_all/
    output
        dataframe
    """
    #load models
    mr = load_model('gender_residual',path)
    pca = load_model('pca_95','/'.join([path,snp]))
    #prediction
    pbase = model_pred(mr,[meta[['gender','nationality']].values])[0]
    Xt = pca.transform(gt)
    pred = odict()
    for m in ms:
        ml = load_model(m,'/'.join([path,snp]))
        if 'linear' in m:
            pr = model_pred(ml,[Xt])[0]
            ph = pbase + pr
            pred[m] = ph
            if 'height' in meta.columns:
                print('R2',metrics.r2_score(meta.height,ph),'MSE',metrics.mean_squared_error(meta.height,ph))
                #plot scatter 
        else:
            pbi = model_pred(ml,[Xt])[0]
            pred[m] = pbi
            
    
    return pd.DataFrame(pred,index=meta.index)

def height_feature(refb,p,lift=None):
    """
    prepare input data for prediction
    input: reference bim, prediction plink, p_sub_bim over file from hg38 to hg19
    output: the feature of genotype for prediction model
    """
    #23 and 25 both are chrX
    p_bim = p[0]
    if '25' not in set(p_bim.chrom):
        refb['chrom'].replace('25','23',inplace=True)
    if lift is None:
        bim_filter = p_bim.apply(select_bim, axis=1)
        p=plink_slice(p,boo2idx(bim_filter))
        p_sub_bim=p[0]
        p_sub_bim['i']=list(range(len(p_sub_bim)))
    else:
        p_sub_bim = p_bim.iloc[lift.score,:]
        p_sub_bim['pos'] = list(lift.end)
        p_sub_bim['i'] = list(lift.score)

    print('snp match ratio :',len(set(refb.pos).intersection(p_sub_bim.pos))/len(refb)) #higher is better
    #print('higher is better. if the ration is less than 0.8, the result is not reliable')

    result = pd.merge(refb, p_sub_bim, how='left', on=['chrom','pos'])
    tt_sub = plink_slice(p,list(result.i_y.dropna().astype('int')))

    idx_notna = list(result.index[result.i_y.notna()])
    idx_na = list(result.index[result.i_y.isna()])
    refb_in=refb.iloc[idx_notna,:]
    g_in=alignp2b_a0_a1(refb_in,tt_sub)
    g_in.fillna(-1,inplace=True)
    g_in.index = idx_notna
    #filled with mix genotype 1.
    g_out = pd.DataFrame(1.0, index=idx_na, columns=np.arange(g_in.shape[1]))
    df = g_in.append(g_out)
    df = df.sort_index()  # sorting by index
    return df.T

def select_bim(arr):
    if arr[0]==0 or arr[0]==26:
        return 0
    elif arr[4]=='A' and arr[5]=='T':
        return 0
    elif arr[4]=='T' and arr[5]=='A':
        return 0
    elif arr[4]=='G' and arr[5]=='C':
        return 0
    elif arr[4]=='C' and arr[5]=='G':
        return 0
    else:
        return 1

def alignp2b_a0_a1(refb,p):
    refbim=refb.reset_index(drop=True)
    pbim=p[0].reset_index(drop=True)
    diff_idx=[x!=y for x,y in zip(refbim.a0,pbim.a0)]
    print('the ratio of genotype shift',sum(diff_idx)/len(refbim))
    if all([x==y for x,y in zip(refbim.a0[diff_idx],pbim.a1[diff_idx])]):
        print("Success genotype shift")
    else:
        print("warning: some genotype (a0 and a1) not match")
    pbed = pd.DataFrame(p[2].compute())
    tmp = pbed.copy()
    #tmp1 = pbed.copy()
    tmp = tmp.iloc[diff_idx,:]
    #tmp1 = tmp1.iloc[diff_idx,:]
    tmp.replace([0.0,2.0],[2.0,0.0],inplace=True)
    pbed.loc[diff_idx]=tmp
    return pbed


def rsqure(yh,y):
    yhat = yh                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    r2 = ssreg / sstot
    return r2


def clamp(n, smallest, largest): return max(smallest, min(n, largest))
def triclip(pr,pb,pt):
    new_pr=[]
    for i,r in enumerate(pr):
        if pb[i]:
            if pt[i]==1:
                new_pr.append(max(5,r))
            elif pt[i]==0:
                new_pr.append(max(0,r))
            else:
                new_pr.append(0)
        else:
            if pt[i]==1:
                new_pr.append(0)
            elif pt[i]==0:
                new_pr.append(min(0,r))
            else:
                new_pr.append(min(-5,r))
    return np.array(new_pr)
def biclip(pr,pb):
    new_pr=[]
    for i,r in enumerate(pr):
        if pb[i]:
            new_pr.append(max(0,r))
        else:
            new_pr.append(min(0,r))
    return np.array(new_pr)

def model_pred(ml,dt,p=False):
    pr= []
    if p:
        for d in dt:
            pr.append(ml.predict_proba(d)[:,1])
    else:
        for d in dt:
            pr.append(ml.predict(d))
    return pr

def boo2idx(t):
    return [i for i, x in enumerate(t) if x]

def gn_residual(meta,test=None,t=True):
    #Load training data
    x_train = meta[['gender','nationality']]
    y_train = meta.height

    # Instantiate a linear model
    regr = linear_model.LinearRegression()
    regr.fit(x_train, y_train)
    if t:
        x_test = test[['gender','nationality']]
        # Predict the heights for the test set
        heights = regr.predict(x_test)
        submit(heights)
    return regr

def plink_slice(p,pb=None,pf=None):
    """
    p: list of bim, fam, bed
    pb: index of bim
    pf: index of fam
    """
    (bim,fam,bed)=p
    if pb:
        bim = bim.iloc[pb]
        bed = bed[pb,:]
    if pf:
        fam = fam.iloc[pf]
        bed = bed[:,pf]
    return(bim,fam,bed)

def select_snp(arr):
    if arr[0]==0 or arr[0]==26:
        return 0
    elif arr[4]<0.05 or arr[4]>0.95:
        return 0
    elif arr[2]=='A' and arr[3]=='T':
        return 0
    elif arr[2]=='T' and arr[3]=='A':
        return 0
    elif arr[2]=='G' and arr[3]=='C':
        return 0
    elif arr[2]=='C' and arr[3]=='G':
        return 0
    else:
        return 1


def _read_csv(fn, header):
    return pd.read_csv(
        fn,
        delim_whitespace=True,
        header=None,
        names=header.keys(),
        dtype=header,
        compression=None,
        engine="c",
    )


def read_bim(fn):
    header = odict(
        [
            ("chrom", bytes),
            ("snp", bytes),
            ("cm", float),
            ("pos", int),
            ("a0", bytes),
            ("a1", bytes),
        ]
    )
    df = _read_csv(fn, header)

    df["chrom"] = df["chrom"].astype("category")
    df["a0"] = df["a0"].astype("category")
    df["a1"] = df["a1"].astype("category")
    df["i"] = range(df.shape[0])
    return df


def read_fam(fn):
    header = odict(
        [
            ("fid", str),
            ("iid", str),
            ("father", str),
            ("mother", str),
            ("gender", bytes),
            ("trait", str),
        ]
    )

    df = _read_csv(fn, header)

    df["gender"] = df["gender"].astype("category")
    df["i"] = range(df.shape[0])
    return df
    
def write_bim(bim,path):
    bim.iloc[:,0:6].to_csv(path,sep='\t',index=False,header=False)

def write_fam(fam,fn):
    pass    
    
def submit(x):
    challenge = crowdai.Challenge("OpenSNPChallenge2017", "02f5a9824254b85a068fbdfbb1298b8e")
    challenge.submit(x.tolist())
    challenge.disconnect()
    
    
def init_args(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype', help='Genotype input file, the prefix of plink files (.bim,.fam,.bed)')
    parser.add_argument('-s', '--specimen', help='the sample information')
    parser.add_argument('-l', '--liftover', help='the output file with suffix .Hg19 or .Hg38 from LiftoverCmd.py',default=None)
    #parser.add_argument('-c', '--chain', help='Liftover Chain file',choices=['hg38ToHg19','hg19ToHg38'],default='hg38ToHg19')
    parser.add_argument('-m', '--model', help='a list of prediction models',choices=['linear_regression','logistic_avg','logistic_25perc','rfc_avg','rfc_25perc','rfc_15perc'],default=['linear_regression'],nargs='+')
    parser.add_argument('--snp-list', help='Choose type of snp list', choices=['top200','top200_all','top300','top300_all','top500','top500_all'], default='top200_all')
    parser.add_argument('--model-path',help='Choose the way of model generation',choices=['full_model','rep10cv10_model'], default= 'full_model')
    args = parser.parse_args(arguments)
    return args

def main(args=None):
    args = init_args(args)
    height_prediction(args.genotype,args.liftover,args.specimen,ms=args.model,snp=args.snp_list,path=args.model_path)
        
if __name__ == '__main__':
    main()