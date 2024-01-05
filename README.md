# multiGMF
To enable better integration of various drug and disease similarity data for mining potential drug-disease associations, we propose a new approach based on multi-similarity geometric matrix factorization (multiGMF). It does not simply concatenate and fuse multiple similarities but instead incorporates each similarity into a mathematical model for optimization. First, we initiate the process by utilizing the weighted k-nearest neighbors (WKNN) technique to increase the overall density of the drug-disease association matrix, consequently aiding in predictive analysis. Subsequently, we apply the matrix factorization technique to break down the updated association matrix into feature matrices for drugs and diseases. Leveraging these feature matrices, we incorporate soft regularization and graph regularization to linearly fuse the multiple similarity matrices. It yields a more reliable feature matrix for drugs and diseases. Finally, the unknown drug- disease associations can be inferred by multiplying drug features and disease features.

# Requirements
* Matlab >= 2020

# Installation
ITRPCA can be downloaded by
```
git clone https://github.com/YangPhD84/multiGMF
```
Installation has been tested in a Windows platform.

# Datasets Description
* Gold_standard_dataset.mat: this file contains information about the gold standard dataset;
* CTDdataset2023.mat: this file contains information about the gold standard dataset;
* drug_ChemS: chemical structure similarity matrix;
* drug_AtcS: drug's ATC code similarity matrix;
* drug_SideS: side-effect similarity matrix;
* drug_DDIS: drug-drug interaction similarity matrix;
* drug_TargetS: drug's target profile similarity matrix;
* disease_PhS: disease phenotype similarity matrix;
* disease_DoS: disease ontology similarity matrix;
* didr: disease-drug association matrix.
