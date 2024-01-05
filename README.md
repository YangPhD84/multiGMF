# multiGMF
To enable better integration of various drug and disease similarity data for mining potential drug-disease associations, we propose a new approach based on multi-similarity geometric matrix factorization (multiGMF). It does not simply concatenate and fuse multiple similarities but instead incorporates each similarity into a mathematical model for optimization. First, we initiate the process by utilizing the weighted k-nearest neighbors (WKNN) technique to increase the overall density of the drug-disease association matrix, consequently aiding in predictive analysis. Subsequently, we apply the matrix factorization technique to break down the updated association matrix into feature matrices for drugs and diseases. Leveraging these feature matrices, we incorporate soft regularization and graph regularization to linearly fuse the multiple similarity matrices. It yields a more reliable feature matrix for drugs and diseases. Finally, the unknown drug- disease associations can be inferred by multiplying drug features and disease features.

# Requirements
* Matlab >= 2020

# Installation
multiGMF can be downloaded by
```
git clone https://github.com/YangPhD84/multiGMF
```
Installation has been tested in a Windows platform.

# Datasets Description
* Gold_standard_dataset.mat: this file contains information about the gold standard dataset;
* CTDdataset2023.mat: this file contains information about our newly organized CTDdataset2023 dataset;
* drug_ChemS: chemical structure similarity matrix;
* drug_AtcS: drug's ATC code similarity matrix;
* drug_SideS: side-effect similarity matrix;
* drug_DDIS: drug-drug interaction similarity matrix;
* drug_TargetS: drug's target profile similarity matrix;
* disease_PhS: disease phenotype similarity matrix;
* disease_DoS: disease ontology similarity matrix;
* didr: disease-drug association matrix.

# Instructions
We provide detailed step-by-step instructions for running multiGMF model.

**Step 1**: add datasets paths
```
addpath('Datasets');
```

**Step 2**: load datasets with association matirx and similarity matrices
```
load  Gold_standard_dataset
Wrr1 = drug_ChemS;
Wrr2 = drug_AtcS;
Wrr3 = drug_SideS;
Wrr4 = drug_DDIS;
Wrr5 = drug_TargetS;
Wrr=(Wrr1+Wrr2+Wrr3+Wrr4+Wrr5)/5;
R = {Wrr1, Wrr2, Wrr3, Wrr4, Wrr5};
Wdd1 = disease_PhS;
Wdd2 = disease_DoS;
Wdd=(Wdd1+Wdd2)/2;
D = {Wdd1, Wdd2};
Wdr = didr;
Wrd = Wdr';
[dn, dr] = size(Wdr);
min_mn = min(dn, dr);
```

**Step 3**: parameter Settings

The hyper-parameters are fixed.
```
lambda1 = 0.0001;
lambda2 = lambda1;
lambda3 = 1;
r = 0.9;
k = 10;
tau = 0.7;
MaxIter = 300;
tol1 = 2*1e-3;
tol2 = 1*1e-4;
```

**Step 4**: run the multi-similarity geometric matrix factorization (multiGMF)
```
[W, H, iter] = fmultiGMF(Wrd,R,D,tau,Wrr,Wdd,r,k,MaxIter,lambda1,lambda2,lambda3,tol1,tol2);
M_recovery = W * H';
```

# A Quickstart Guide
Users can immediately start playing with multiGMF running ``` Demo_multiGMF.m ``` in matlab.
* ```Demo_multiGMF.m```: it demonstrates a process of predicting drug-disease associations on the gold standard dataset (Gold_standard_dataset) by multiGMF algorithm.

# Contact
If you have any questions or suggestions with the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
