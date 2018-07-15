
# coding: utf-8

# In[251]:


from multiprocessing import Pool
import subprocess
import sys
import os
import gzip
import pandas as pd
import time
import re


# In[252]:


file2 = sys.argv[1] ## dir to sample genome vcf
file3 = 'GWAS_snps.txt' ### files contained snps you need
file1 = re.sub('.*/',  '', file2)
file1 = re.sub('.vcf.gz', '_GWAS.txt', file1)


# In[254]:


if os.path.isfile(file1)==False:
    pass
else:
    os.remove(file1)


# In[255]:


fp1 = gzip.open(file2, 'r')
for line in fp1:
    line=line.decode()
    if "##" in line:
        pass
    elif "CHROM" in line:
        header = line
    else:
        break
fp1.close()
fp1 = open(file1,'a+')
fp1.write(header)
fp1.close()


# In[256]:


fp = pd.read_csv(file3,sep='\t',header=None, index_col=None)
find_idx = fp[1].astype(str).str.cat(fp[2].astype(str), sep=':')
find_idx = find_idx.str.cat(fp[2].astype(str), sep='-')


# In[257]:


def write_in(x):
    with open(file1, 'a+') as f:
        f.write(x)


def tabix_find(index):
    new_line = subprocess.check_output(['tabix',file2,index])
    new_line = new_line.decode()
    return new_line

if __name__ == '__main__':
    e1 = time.time()
    pool = Pool(8)

    for i in find_idx:
        pool.apply_async(tabix_find, (i,), callback=write_in)

    pool.close()
    pool.join()
    e2 = time.time()
    print(float(e2-e1))


# In[258]:


import pandas as pd
import numpy as np
"""
read in sample genome vcf
df_snp: the information of each snps
df_sample: samples' genotypes for each snps
0: ref homozygous; 1: heterozygous; 2: alt homozygous; -1: missing
"""

fname = file1

df = pd.read_csv(fname, sep='\t')


# In[260]:


df_snp = df.iloc[:,[2,0,1,3,4]]
df_snp = df_snp.copy()
df_snp.columns = ['ID','CHROM', 'POS', 'REF', 'ALT']
df_snp.loc[:,'identifier_s'] = df_snp['CHROM'].astype(str).str.cat(df_snp['POS'].values.astype(np.str),sep=':')
df_sample = df.iloc[:,9:]
df_sample = df_sample.copy()
df_sample.index = df_snp['identifier_s'].values
def split_str(ele):
    ele = ele.split(':')[0]
    if ele=='0|0':
        a = 0
    elif ele=='0|1' or ele=='1|0':
        a = 1
    elif ele=='1|1':
        a=2
    else:
        a=-1
    return a

df_sample = df_sample.applymap(split_str)
df_sample = df_sample.T
rowname = df_sample.index
colname = df_sample.columns


# In[262]:


"""
Impute missing allele (most frequent)
if all missing, then 1 will be return
"""
for i in df_sample.columns:
    count = df_sample[i].value_counts()
    best = count.argmax()
    if best == -1:
        best = 1
    df_sample.loc[df_sample[i]==-1, i] = best


# In[264]:


"""
compare the differences between sample snps information and markers' info
1. whether the a1, a2 matched
2. whether the genotype needs to change
"""

marker = pd.read_csv('GWAS_snps.txt',sep='\t',
                     low_memory=False, index_col=None ,header=None)


# In[265]:


marker.columns = ['ID','chr','pos','Allele1', 'Allele2','Freq']


# In[266]:


marker.loc[:,'identifier_m'] = marker['chr'].astype(str).str.cat(marker['pos'].astype(str),sep=':')

marker = marker.merge(df_snp, how='left', left_on = 'identifier_m', right_on = 'identifier_s')


# In[226]:


missing_snp = (marker.loc[np.isnan(marker['CHROM']),'identifier_m'].values).tolist()


# In[268]:


if len(missing_snp)!=0:
    for i in missing_snp:
        df_sample[i] = 1


# In[269]:


marker = marker[['identifier_m','ID_y','chr','POS','Allele1','Allele2','REF','ALT']]
marker = marker.copy()
marker.loc[:,'matched'] = 0
marker.loc[:,'need_to_convert'] = 0


# In[271]:


def determine_alleles(arr):
    cal1_AG = int(arr['Allele1'] in ['A','G']) + int(arr['Allele2'] in ['A','G'])
    cal1_AC = int(arr['Allele1'] in ['A','C']) + int(arr['Allele2'] in ['A','C'])
    cal1_GT = int(arr['Allele1'] in ['T','G']) + int(arr['Allele2'] in ['T','G'])
    cal1_CT = int(arr['Allele1'] in ['C','T']) + int(arr['Allele2'] in ['C','T'])
    cal2_AG = int(arr['REF'] in ['A','G']) + int(arr['ALT'] in ['A','G'])
    cal2_AC = int(arr['REF'] in ['A','C']) + int(arr['ALT'] in ['A','C'])
    cal2_GT = int(arr['REF'] in ['T','G']) + int(arr['ALT'] in ['T','G'])
    cal2_CT = int(arr['REF'] in ['C','T']) + int(arr['ALT'] in ['C','T'])
    if cal1_AG==2:
        if (cal2_AG==2) or (cal2_CT==2):
            return 'matched'
        else:
            return 'not_matched'
    elif cal1_AC==2:
        if (cal2_AC==2) or (cal2_GT==2):
            return 'matched'
        else:
            return 'not_matched'
    elif cal1_GT==2:
        if (cal2_GT==2) or (cal2_AC==2):
            return 'matched'
        else:
            return 'not_matched'
    elif cal1_CT==2:
        if (cal2_CT==2) or (cal2_AG==2):
            return 'matched'
        else:
            return 'not_matched'
    else:
        return 'not_matched'


# In[272]:


marker.loc[:,'matched'] = marker.apply(determine_alleles, axis=1)


# In[273]:


def need_to_convert(arr):
    if arr['matched']=='matched':
        if arr['Allele1']==arr['ALT']:
            return 'no_need'
        elif arr['Allele1']==arr['REF']:
            return 'need'
        elif arr['Allele1']=='A' and arr['ALT']=='T':
            return 'no_need'
        elif arr['Allele1']=='T' and arr['ALT']=='A':
            return 'no_need'
        elif arr['Allele1']=='C' and arr['ALT']=='G':
            return 'no_need'
        elif arr['Allele1']=='G' and arr['ALT']=='C':
            return 'no_need'
        else:
            return 'need'
    elif arr['matched']=='not_matched':
        return 'not_valid'


# In[274]:


marker.loc[:,'need_to_convert'] = marker.apply(need_to_convert, axis=1)

not_matched_snps = marker.loc[marker.need_to_convert == 'not_valid', 'identifier_m'].values.tolist()
if len(not_matched_snps)!=0:
    for i in not_matched_snps:
        df_sample[i] = 1


# In[275]:


df_sample = df_sample[marker.identifier_m]
df_sample = df_sample.copy()
# In[279]:


for i in range(df_sample.shape[1]):
    if marker.iloc[i,9] == 'need':
        df_sample.is_copy = None
        df_sample.iloc[:,i] = 2 - df_sample.iloc[:,i]


# In[283]:

""" Model Part """
import pickle
from sklearn.model_selection import cross_val_predict, cross_val_score,StratifiedKFold,train_test_split
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression, LinearRegression


# In[284]:
predicted_proba = pd.DataFrame(columns=['top1'])
predicted = pd.DataFrame(columns=['top1'])

for i in range(1,11):
    clfCV_m = pickle.load(open('logistic_m_rmXchrom_top'+str(i)+'.pickle', 'rb'))
    s = 'top' + str(i)
    predicted_proba[s] = clfCV_m.predict_proba(df_sample)[:,1]
    predicted[s] = clfCV_m.predict(df_sample)

# In[290]:

predicted_proba['final'] = (predicted_proba['top1']*0.7+predicted_proba['top2']*0.67+\
predicted_proba['top3']*0.65+predicted_proba['top4']*0.65+predicted_proba['top5']*0.64+\
predicted_proba['top6']*0.63+predicted_proba['top7']*0.62+predicted_proba['top8']*0.61+\
predicted_proba['top9']*0.61+predicted_proba['top10']*0.60)/(0.7+0.67+0.65+0.65+0.64+0.63+0.62+\
                                                             0.61+0.61+0.6)

predicted_proba['final_class'] = (predicted_proba['final'] > 0.5).astype(int)
predicted_proba.index = df_sample.index
predicted.index = df_sample.index


# In[293]:

if not os.path.exists('result'):
    os.mkdir('result')

if not os.path.exists('result_class'):
    os.mkdir('result_class')

result_file = 'result/' + file1
predicted_proba.to_csv(result_file,header=True,sep='\t' ,index=True)
result_file = 'result_class/' + file1
predicted.to_csv(result_file,header=True,sep='\t' ,index=True)
