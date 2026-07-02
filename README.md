# multiGMF
To enable better integration of various drug and disease similarity data for mining potential drug-disease associations, we propose a new approach based on multi-similarity geometric matrix factorization (multiGMF). It does not simply concatenate and fuse multiple similarities but instead incorporates each similarity into a mathematical model for optimization. First, we use the weighted k-nearest neighbors (WKNN) technique to increase the density of the drug-disease association matrix. Subsequently, we apply matrix factorization to decompose the updated association matrix into latent feature matrices for drugs and diseases. Soft regularization and graph regularization are then used to integrate multiple similarity matrices. Finally, unknown drug-disease associations can be inferred by multiplying the learned drug and disease feature matrices.

# Requirements
* Matlab >= 2014

# Installation
multiGMF can be downloaded by
```
git clone https://github.com/YangPhD84/multiGMF
```
Installation has been tested on a Windows platform.

# Repository Structure
* `Datasets/`: benchmark datasets used in the manuscript, including `Fdataset.mat`, `Cdataset.mat`, and `CTDdataset2023.mat`.
* `Demo_multiGMF.m`: demo script for final prediction and case-study candidate ranking after model evaluation and parameter selection.
* `multiGMF_10CV.m`: script for 10-fold cross-validation evaluation.
* `multiGMF_Denovo.m`: script for drug-side cold-start evaluation.
* `fmultiGMF.m`: main implementation of the full multiGMF model.
* `no_soft_fmultiGMF.m`: implementation of the no-soft ablation variant.
* `WKNN.m`: WKNN preprocessing function.
* `Fun_Auc3.m`: function for calculating AUC, AUPR, and Precision in the evaluation scripts.
* `Baseline/`: contains 5 baseline methods implemented in MATLAB. Each baseline script can be run for both 10-fold cross-validation and cold-start (denovo) experiments to generate comparative metrics against multiGMF.

# Evaluation Scripts vs Final-Prediction Demo

The evaluation scripts and the final-prediction demo serve different purposes.

* `multiGMF_10CV.m` is used to reproduce the repeated 10-fold cross-validation experiments reported in the manuscript. In this script, test associations are masked before WKNN preprocessing, model training, hyperparameter selection, and final evaluation.
* `multiGMF_Denovo.m` is used to reproduce the drug-side cold-start experiments reported in the manuscript. The script name is retained for compatibility, but the experiment corresponds to the drug-side cold-start setting.
* `Demo_multiGMF.m` is used only for final full-data matrix completion and case-study candidate ranking after model evaluation and parameter selection. It should not be used to reproduce the repeated 10-fold cross-validation or drug-side cold-start results reported in the manuscript.

# Datasets Description
* `Fdataset.mat`: the Fdataset used in the manuscript.
* `Cdataset.mat`: the Cdataset used in the manuscript.
* `CTDdataset2023.mat`: the CTDdataset2023 dataset used in the manuscript.
* `drug_ChemS`: chemical structure similarity matrix.
* `drug_AtcS`: drug ATC code similarity matrix.
* `drug_SideS`: side-effect similarity matrix.
* `drug_DDIS`: drug-drug interaction similarity matrix.
* `drug_TargetS`: drug target profile similarity matrix.
* `disease_PhS`: disease phenotype similarity matrix.
* `disease_DoS`: disease ontology similarity matrix.
* `didr`: disease-by-drug association matrix.

# Instructions
We provide step-by-step instructions for running the multiGMF model.

**Step 1**: add the dataset path
```
addpath('Datasets');
```

**Step 2**: load one benchmark dataset and similarity matrices, and prepare association matrices
```
load Fdataset
% load Cdataset
% load CTDdataset2023

Wrr1 = drug_ChemS;
Wrr2 = drug_AtcS;
Wrr3 = drug_SideS;
Wrr4 = drug_DDIS;
Wrr5 = drug_TargetS;
Wrr = (Wrr1 + Wrr2 + Wrr3 + Wrr4 + Wrr5) / 5;
R = {Wrr1, Wrr2, Wrr3, Wrr4, Wrr5};

Wdd1 = disease_PhS;
Wdd2 = disease_DoS;
Wdd = (Wdd1 + Wdd2) / 2;
D = {Wdd1, Wdd2};

Wdr = didr;     % disease-by-drug association matrix
Wrd = Wdr';     % drug-by-disease association matrix used by multiGMF
[dn, dr] = size(Wdr);
min_mn = min(dn, dr);
```

**Step 3**: parameter settings
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

**Step 4**: run the multi-similarity geometric matrix factorization model
```
rankk = floor(min_mn * tau);

Input.A = R;
Input.B = D;
Input.X = Wdr;
Input.WInit = rand(dr, rankk);
Input.HInit = rand(dn, rankk);
Input.kk = k;
Input.Wrr = Wrr;
Input.Wdd = Wdd;
Input.P_TMat = Wrd;

Options.MaxIter = MaxIter;
Options.lambda_soft = 1;
Options.lambda1 = lambda1;
Options.lambda3 = lambda3;
Options.mu1 = 1;
Options.mu2 = 1;
Options.tol1 = tol1;
Options.tol2 = tol2;

Output = fmultiGMF(Input, Options);

W = Output.W;
H = Output.H;
iter = Output.t;

M_recovery = W * H';
```

# Reproducing Experiments
* Run `multiGMF_10CV.m` to reproduce the 10-fold cross-validation experiments.
* Run `multiGMF_Denovo.m` to reproduce the drug-side cold-start experiments. The script name is retained for compatibility, but the experiment corresponds to the cold-start setting in the manuscript.
* AUC, AUPR, and Precision are calculated in the evaluation scripts using `Fun_Auc3.m`.
* The no-soft ablation can be reproduced by calling `no_soft_fmultiGMF.m` in the 10-fold cross-validation or cold-start script.
* The no-GrLap ablation can be reproduced by setting the graph-regularization weights to zero, i.e., lambda1 = lambda2 = 0, while retaining WKNN preprocessing, multi-source similarity integration, soft coupling, and all other settings.
* Single-similarity ablation experiments can be reproduced by following the commented settings in `multiGMF_10CV.m` or `multiGMF_Denovo.m`.
* Case-study candidate rankings are generated using `Demo_multiGMF.m` after parameter selection and model evaluation. Known associations in the benchmark matrix should be excluded from the ranked candidate list.
* Baseline methods in `Baseline/` can be run similarly for 10CV and cold-start experiments.
* The full Holm-adjusted p-value tables are provided in Supplementary Tables S2-S7.

# Mapping between Manuscript Tables and Scripts

| Manuscript item | Dataset | Evaluation setting | Script or setting | Key parameters |
|---|---|---|---|---|
| Table 2 | Fdataset | Internal validation for lambda1 and tau | Parameter search within `multiGMF_10CV.m` | k = 10; tau in {0.1, 0.3, 0.5, 0.7, 0.9}; lambda1 = lambda2 in {0.1, 0.01, 0.001, 0.0001, 0.00001} |
| Table 3 | Fdataset | Internal validation for k | Parameter search within `multiGMF_10CV.m` | tau = 0.7; lambda1 = lambda2 = 0.0001; k in {1, 5, 10, 15, 20} |
| Table 4 | Fdataset | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m` | tau = 0.7; lambda1 = lambda2 = 0.0001; k = 10 |
| Table 5 | Fdataset | Ablation and single-source sensitivity | `multiGMF_10CV.m`; `multiGMF_Denovo.m`; `no_soft_fmultiGMF.m`; single-source settings | Same fold-wise masking protocol |
| Table 6 | Fdataset | Case-study candidate ranking | `Demo_multiGMF.m` | Final prediction after parameter selection; known associations excluded from candidate ranking |
| Table 7 | Cdataset | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m` | tau = 0.7; lambda1 = lambda2 = 0.0001; k = 10 |
| Table 8 | CTDdataset2023 | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m` | tau = 0.7; lambda1 = lambda2 = 0.01; k = 10 |
| Supplementary Table S1 | All datasets | no-GrLap ablation | Set graph-regularization weights to zero while retaining WKNN, multi-source integration, and soft coupling | lambda1 = lambda2 = 0 |
| Supplementary Tables S2-S7 | All datasets | Holm-corrected paired t-tests | Supplementary Information | Paired two-sided t-tests with Holm correction |

# A Quickstart Guide
Users can immediately start using multiGMF by running `Demo_multiGMF.m` in MATLAB.
* `Demo_multiGMF.m`: demonstrates final prediction and candidate ranking on the Fdataset using the multiGMF algorithm.

# Contact
If you have any questions or suggestions about the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
