import numpy as np
import pandas as pd
from sklearn import linear_model


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
    
def submit(x):
    challenge = crowdai.Challenge("OpenSNPChallenge2017", "02f5a9824254b85a068fbdfbb1298b8e")
    challenge.submit(x.tolist())
    challenge.disconnect()