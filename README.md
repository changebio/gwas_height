Predict Height
==============

Predict height by gender and genotypes.
The models are trained on hg19 data.
It can also work well on hg38 data
These are many kinds of models by different snp cutoff(top200,top200_all,...) and training methods(full_model). 
predict_height/
├── full_model
│   ├── gender_residual.model
│   ├── top200
│   │   ├── linear_regression.model
│   │   ├── logistic_25perc.model
│   │   ├── logistic_avg.model
│   │   ├── pca_95.model
│   │   └── ref.bim
│   ├── top200_all
│   ├         …
├── lib
│   ├── hg19ToHg38.over.chain.gz
│   ├── hg38ToHg19.over.chain.gz
│   ├── liftOver
│   └── maf.csv
├── LiftoverCmd.py
├── PlinkCmd.py
├── predhpipeline.py
├── snp2h.py
├── README.md

# Usage

### Prepare data ( for example, using temp as file prefix )

#### Genotype 
1. temp.vcf or temp.vcf.gz
2. temp.ped and temp.map 
3. temp.bim, temp.fam, and temp.bed

#### LiftOver
temp.Hg19

#### Specimen
sample information is not necessary. if there is this information, it file should have at least one column (`gender`, 1 as male, 2 as female) and separate by Tab. The order of the rows must be the same as the samples in Genotype file. If height column is given, model will output R2 and MSE for regression or AUC and F1 score for classification.
for example, meta.csv 
Specimen file as below: 1 column (gender)
```
	gender
TL2257	1
TL0067	2
TL0139	2
TL0279	2
TL1090	2
```
or Specimen file as below: 2 columns (gender, height)
```
	gender	height
TL2257	1	177
TL0067	2	155
TL0139	2	143.5
TL0279	2	161.9
TL1090	2	155.2
```

### Prediction
There are two ways to use the scripts: pipeline and step by step.

#### Pipeline
input is genotype file (for example: temp.vcf.gz)
output is a csv file with the same prefix as the genotype file (temp.csv)
```
python predhpipeline.py -g <genotype_file> -c <liftover_chain>
```
for example 
1.pipeline without liftover
```
python predhpipeline.py -g ddat/crowdAI/temp.vcf.gz 
#or
python predhpipeline.py -g ddat/crowdAI/temp.vcf.gz -m linear_regression logistic_avg logistic_25perc --snp-list top200_all --model-path full_model
#or for plink file(.ped and .map) or binary plink file(.bim, .fam and .bed). temp is the prefix of plink files.
python predhpipeline.py -g ddat/crowdAI/temp
```

2.pipeline with liftover
```
python predhpipeline.py -g ddat/crowdAI/temp.vcf.gz -c hg38ToHg19
#or
python predhpipeline.py -g ddat/crowdAI/temp.vcf.gz -c hg38ToHg19 -m linear_regression --snp-list top200_all --model-path full_model
#or for plink file(.ped and .map) or binary plink file(.bim, .fam and .bed). temp is the prefix of plink files.
python predhpipeline.py -g ddat/crowdAI/temp -c hg38ToHg19 -m linear_regression --snp-list top200_all --model-path full_model
```

#### Step by step
1. Genotype convert from .vcf or .ped and .pad to binary plink (.bim,.fam,.bed)
```
#for vcf
python PlinkCmd.py -g ddat/crowdAI/temp.vcf.gz
#for plink file(.ped and .pad) input is prefix
python PlinkCmd.py -g ddat/crowdAI/temp
```
2. Liftover hg38 to hg19. If It is hg19 data, skip this step and go to 3.1 Height prediction without liftover
```
python LiftoverCmd.py -g ddat/crowdAI/temp -c hg38ToHg19
```
3. Height prediction
```
python snp2h.py -g ddat/crowdAI/temp -l ddat/crowdAI/temp.Hg19
#or
python snp2h.py -g ddat/crowdAI/temp -l ddat/crowdAI/temp.Hg19 -m linear_regression --snp-list top200_all --model-path full_model
```
3.1 Height prediction without liftover
```
python snp2h.py -g ddat/crowdAI/temp
#or
python snp2h.py -g ddat/crowdAI/temp -m linear_regression logistic_avg logistic_25perc --snp-list top200_all --model-path full_model
```

# Requirements
Python 3.6.5
sklearn=0.18.2
numpy=1.14.3=py36h28100ab_2
pandas=0.22.0=py36_1
pandas-plink=1.2.25=py36_0
plink=1.90b4=1
scipy=1.1.0=py36hfc37229_0

install plink
```
conda install -c bioconda plink
```
