# multiGMF
To enable better integration of various drug and disease similarity data for mining potential drug-disease associations, we propose a new approach based on multi-similarity geometric matrix factorization (multiGMF). It does not simply concatenate and fuse multiple similarities but instead incorporates each similarity into a mathematical model for optimization. First, we use the weighted k-nearest neighbors (WKNN) technique to increase the density of the drug-disease association matrix. Subsequently, we apply matrix factorization to decompose the updated association matrix into latent feature matrices for drugs and diseases. Soft regularization and graph regularization are then used to integrate multiple similarity matrices. Finally, unknown drug-disease associations can be inferred by multiplying the learned drug and disease feature matrices.

# Requirements

- MATLAB R2017b or later.
- Statistics and Machine Learning Toolbox, required by `Ttest_10CV.m` and `Ttest_Denovo.m` for the paired `ttest` function.

# Installation
multiGMF can be downloaded by
```bash
git clone https://github.com/YangPhD84/multiGMF
```
Installation has been tested on a Windows platform.

# Repository Structure
* `Datasets/`: benchmark datasets used in the manuscript, including `Fdataset.mat`, `Cdataset.mat`, and `CTDdataset2023.mat`.
* `Demo_multiGMF.m`: demo script for final full-data matrix completion after model evaluation and parameter selection. The resulting score matrix can be used for case-study candidate ranking after known associations are excluded.
* `multiGMF_10CV.m`: script for repeated 10-fold cross-validation evaluation.
* `multiGMF_Denovo.m`: script for drug-side cold-start evaluation.
* `fmultiGMF.m`: main implementation of the full multiGMF model.
* `no_soft_fmultiGMF.m`: implementation of the no-soft ablation variant.
* `WKNN.m`: WKNN preprocessing function.
* `Fun_Auc3.m`: function for calculating AUC, AUPR, and Precision in the evaluation scripts.
* `Baseline/`: contains 5 baseline methods implemented in MATLAB. Each baseline script can be run for both 10-fold cross-validation and cold-start (denovo) experiments to generate comparative metrics against multiGMF.
* `t-test/Ttest_10CV.m`: script for generating the Holm-corrected paired t-test results for repeated 10-fold cross-validation.
* `t-test/Denovo/Ttest_Denovo.m`: script for generating the Holm-corrected paired t-test results for drug-side cold-start evaluation.
* `t-test/` and `t-test/Denovo/`: contain the MAT result files and the detailed and summary CSV files used to prepare Supplementary Tables S2-S7.

# Evaluation Scripts vs Final-Prediction Demo

The evaluation scripts and the final-prediction demo serve different purposes.

* `multiGMF_10CV.m` is used to reproduce the repeated 10-fold cross-validation experiments reported in the manuscript. In this script, test associations are masked before WKNN preprocessing, model training, and final evaluation.
* `multiGMF_Denovo.m` is used to reproduce the drug-side cold-start experiments reported in the manuscript. The script name is retained for compatibility, but the experiment corresponds to the drug-side cold-start setting. The full `fmultiGMF.m` model is used as the default configuration.
* `Demo_multiGMF.m` is used only for final full-data matrix completion after model evaluation and parameter selection. It should not be used to reproduce the repeated 10-fold cross-validation or drug-side cold-start results reported in the manuscript.

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
```matlab
addpath('Datasets');
```

**Step 2**: load one benchmark dataset and similarity matrices, and prepare association matrices
```matlab
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
```matlab
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

Here, `lambda1` and `lambda2` control drug-side and disease-side graph regularization, respectively; `lambda3` controls Tikhonov regularization; `r` is the WKNN rank-decay factor; `k` is the number of nearest neighbors; and `tau` controls the latent dimension through `floor(tau * min(dn,dr))`.

**Step 4**: run the multi-similarity geometric matrix factorization model
```matlab
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
Options.lambda2 = lambda2;
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
* Run `multiGMF_10CV.m` to reproduce the repeated 10-fold cross-validation experiments.
* Run `multiGMF_Denovo.m` to reproduce the drug-side cold-start experiments. The script name is retained for compatibility, but the experiment corresponds to the cold-start setting in the manuscript.
* AUC, AUPR, and Precision are calculated in the evaluation scripts using `Fun_Auc3.m`.
* The no-soft ablation can be reproduced by calling `no_soft_fmultiGMF.m` in the 10-fold cross-validation or cold-start script.
* The no-GrLap ablation can be reproduced by setting the graph-regularization weights to zero, i.e., `lambda1 = lambda2 = 0`, while retaining WKNN preprocessing, multi-source similarity integration, soft coupling, and all other settings.
* Single-similarity ablation experiments can be reproduced by following the commented settings in `multiGMF_10CV.m` or `multiGMF_Denovo.m`.
* Case-study candidate rankings are derived from the `M_recovery` score matrix generated by `Demo_multiGMF.m` after parameter selection and model evaluation. Known associations in the benchmark matrix must be excluded before unknown drug-disease pairs are ranked.
* Baseline methods in `Baseline/` can be run similarly for 10CV and drug-side cold-start experiments.

The principal result-file naming conventions are:

* multiGMF 10CV: `<dataset>_multiGMF_10CV_fold_results.mat`.
* Baseline 10CV: `<dataset>_STresult_<method>_10CV_fold_results.mat`.
* multiGMF and baseline cold-start: `<dataset>_STresult_<method>_Denovoone.mat`.

## Regenerating MAT and CSV Outputs

All top-level evaluation scripts use Fdataset as the default active dataset. Run `Demo_multiGMF.m`, `multiGMF_10CV.m`, and `multiGMF_Denovo.m` from the repository root. Run each baseline evaluation script from its corresponding method directory under `Baseline/`, because the baseline scripts use method-specific relative paths.

The evaluation scripts save their MAT result files in the current working directory. To regenerate the statistical-analysis outputs, copy the repeated 10-fold cross-validation MAT files to `t-test/` and the drug-side cold-start MAT files to `t-test/Denovo/`, retaining the documented filenames. Then run `Ttest_10CV.m` and `Ttest_Denovo.m`, respectively.

# Parameter Tuning

The parameter analyses reported in Tables 2 and 3 were performed on Fdataset using one fixed, complete 10-fold cross-validation run. Specifically, the evaluation code performs 10 repetitions of 10-fold cross-validation (`Count_CV = 10` and `nCV = 10`), and the first complete 10-fold cross-validation run (`num = 1`) was used for parameter selection.

For every candidate parameter setting, the same 10-fold partition generated in this fixed run was retained. Each candidate setting was evaluated by completing all 10 folds, with one fold held out in turn and the remaining nine folds used for model fitting. In every fold, the held-out associations were masked before WKNN preprocessing and model training. The run-level `AUC + AUPR` obtained after completing all 10 folds, corresponding to `A_AUC_values(1) + m_A_AUPR_values(1)` in `multiGMF_10CV.m`, was used to compare candidate settings together with consideration of latent-dimensional complexity.

The fold assignments, random seed, association masking procedure, WKNN preprocessing procedure, and evaluation order were kept identical across all candidate parameter settings. Therefore, the parameter analyses were based on one complete 10-fold cross-validation run rather than on the result of a single fold or the average of all 10 repeated runs. The selected parameter settings were subsequently fixed for the reported repeated 10-fold cross-validation and drug-side cold-start evaluations.

For Table 2:

- The first complete 10-fold cross-validation run (`num = 1`) on Fdataset was used for parameter selection.
- `k` was fixed at `10`.
- `lambda1 = lambda2` was selected from `{0.1, 0.01, 0.001, 0.0001, 0.00001}`.
- `tau` was selected from `{0.1, 0.3, 0.5, 0.7, 0.9}`.
- Every candidate combination was evaluated by completing all 10 folds of the same fixed 10-fold partition.
- The same fold assignments, random seed, test-association masking procedure, and WKNN preprocessing procedure were retained across all candidate combinations.
- In each fold, the held-out associations were masked before WKNN preprocessing and model fitting.
- The run-level `AUC + AUPR`, calculated as `A_AUC_values(1) + m_A_AUPR_values(1)` after all 10 folds had been completed, was used to compare candidate settings together with consideration of latent-dimensional complexity.
- The selected values were `lambda1 = lambda2 = 0.0001` and `tau = 0.7`.

For Table 3:

- The same fixed, complete 10-fold cross-validation run (`num = 1`) used for Table 2 was retained.
- `lambda1 = lambda2 = 0.0001` and `tau = 0.7` were fixed.
- `k` was selected from `{1, 5, 10, 15, 20}`.
- Every candidate value of `k` was evaluated by completing all 10 folds of the same fixed 10-fold partition.
- The fold assignments, random seed, masked associations, WKNN preprocessing procedure, and evaluation order were identical across all candidate values.
- The parameter-selection criterion was the run-level `AUC + AUPR`, calculated as `A_AUC_values(1) + m_A_AUPR_values(1)` after all 10 folds had been completed.
- The selected value was `k = 10`, which achieved the highest run-level `AUC + AUPR`.

The final settings used in the reported evaluations are:

| Dataset | lambda_soft | lambda1 | lambda2 | lambda3 | tau | k | WKNN r | MaxIter | tol1 | tol2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fdataset | 1 | 0.0001 | 0.0001 | 1 | 0.7 | 10 | 0.9 | 300 | 0.002 | 0.0001 |
| Cdataset | 1 | 0.0001 | 0.0001 | 1 | 0.7 | 10 | 0.9 | 300 | 0.002 | 0.0001 |
| CTDdataset2023 | 1 | 0.01 | 0.01 | 1 | 0.7 | 10 | 0.9 | 300 | 0.002 | 0.0001 |

# Statistical Analysis for Supplementary Tables S2-S7

Supplementary Tables S2-S7 were prepared from the CSV files generated by `Ttest_10CV.m` and `Ttest_Denovo.m`. The entries in the detailed CSV files provide the method pair, Holm-adjusted p-value, and significance result used in the corresponding Supplementary Table. The star-summary CSV files provide the table-ready mean, standard deviation, and star annotation.

For each dataset, evaluation setting, and metric, the method with the highest reported mean was compared with each of the other five methods using a two-sided paired t-test. The raw p-values from these five comparisons were adjusted using the Holm procedure with `alpha = 0.05`.

* For repeated 10-fold cross-validation, pairing is based on the same repeated data split across methods. `Ttest_10CV.m` reads `A_AUC_values`, `m_A_AUPR_values`, and `Precision_values` from each MAT result file.
* For drug-side cold-start evaluation, pairing is based on the same held-out target drug across methods. `Ttest_Denovo.m` uses `RowAucValue`, `n_RowAuPRValue` (or the documented `RowAuPRValue` fallback), and `RowPrecisionValue` for the paired tests.
* In the cold-start analysis, the reported AUC, AUPR, and Precision point estimates are read from `Result_AUC_value`, `R_m_A_AUPR_value`, and `Result_Precision_value`, respectively. The point estimates and the target-drug-level paired vectors should not be interchanged.
* A comparison is marked significant only when the Holm-adjusted p-value is below 0.05 and the paired mean difference favors the method listed as the best method.
* A star is added to the best-performing method only when it is significantly better than all five compared baseline methods after Holm correction.
* Before running either script, all paired vectors must have equal lengths and preserve the same split or target-drug order across methods.

Run the 10CV analysis from `t-test/`:

```matlab
run('Ttest_10CV.m');
```

The generated files are:

```text
<dataset>_10CV_best_method_paired_ttest_holm_detail.csv
<dataset>_10CV_best_method_star_summary_holm.csv
```

Run the drug-side cold-start analysis from `t-test/Denovo/`:

```matlab
run('Ttest_Denovo.m');
```

The generated files are:

```text
<dataset>_Denovo_best_method_paired_ttest_holm_detail.csv
<dataset>_Denovo_best_method_star_summary_holm.csv
```

Activate only one dataset block at a time in each T-test script.

# Mapping between Manuscript Tables and Scripts

| Manuscript item | Dataset | Evaluation setting | Script or setting | Key parameters | Primary output |
|---|---|---|---|---|---|
| Table 2 | Fdataset | Parameter analysis using one fixed, complete 10-fold cross-validation run | First complete 10-fold CV run (`num = 1`) in `multiGMF_10CV.m` | `k = 10`; `tau` in `{0.1, 0.3, 0.5, 0.7, 0.9}`; `lambda1 = lambda2` in `{0.1, 0.01, 0.001, 0.0001, 0.00001}` | Run-level `AUC + AUPR` after completing all 10 folds |
| Table 3 | Fdataset | WKNN neighborhood-size analysis using the same fixed, complete 10-fold cross-validation run as Table 2 | First complete 10-fold CV run (`num = 1`) in `multiGMF_10CV.m` | `tau = 0.7`; `lambda1 = lambda2 = 0.0001`; `k` in `{1, 5, 10, 15, 20}` | Run-level `AUC + AUPR` after completing all 10 folds |
| Table 4 | Fdataset | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m`; corresponding baseline scripts | tau = 0.7; lambda1 = lambda2 = 0.0001; k = 10 | Fdataset multiGMF and baseline 10CV/cold-start MAT files |
| Table 5 | Fdataset | Ablation and single-source sensitivity | `multiGMF_10CV.m`; `multiGMF_Denovo.m`; `no_soft_fmultiGMF.m`; single-source and no-WKNN settings | Same fold-wise masking protocol and selected Fdataset parameters | AUC, AUPR, and Precision for each model variant |
| Table 6 | Fdataset | Case-study candidate ranking | `Demo_multiGMF.m` | Final prediction after parameter selection; known associations excluded from candidate ranking | `M_recovery` and the ranked unknown drug-disease pairs |
| Table 7 | Cdataset | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m`; corresponding baseline scripts | tau = 0.7; lambda1 = lambda2 = 0.0001; k = 10 | Cdataset multiGMF and baseline 10CV/cold-start MAT files |
| Table 8 | CTDdataset2023 | Repeated 10-fold CV and drug-side cold-start | `multiGMF_10CV.m`; `multiGMF_Denovo.m`; corresponding baseline scripts | tau = 0.7; lambda1 = lambda2 = 0.01; k = 10 | CTDdataset2023 multiGMF and baseline 10CV/cold-start MAT files |
| Supplementary Table S1 | All datasets | no-GrLap ablation | Set graph-regularization weights to zero while retaining WKNN, multi-source integration, and soft coupling | lambda1 = lambda2 = 0 | AUC, AUPR, and Precision for multiGMF and multiGMF (w/o GrLap) |
| Supplementary Table S2 | Fdataset | Repeated 10-fold CV paired tests | `t-test/Ttest_10CV.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `Fdataset_10CV_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |
| Supplementary Table S3 | Fdataset | Drug-side cold-start paired tests | `t-test/Denovo/Ttest_Denovo.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `Fdataset_Denovo_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |
| Supplementary Table S4 | Cdataset | Repeated 10-fold CV paired tests | `t-test/Ttest_10CV.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `Cdataset_10CV_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |
| Supplementary Table S5 | Cdataset | Drug-side cold-start paired tests | `t-test/Denovo/Ttest_Denovo.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `Cdataset_Denovo_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |
| Supplementary Table S6 | CTDdataset2023 | Repeated 10-fold CV paired tests | `t-test/Ttest_10CV.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `CTDdataset2023_10CV_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |
| Supplementary Table S7 | CTDdataset2023 | Drug-side cold-start paired tests | `t-test/Denovo/Ttest_Denovo.m` | Two-sided paired t-test; Holm correction; alpha = 0.05 | `CTDdataset2023_Denovo_best_method_paired_ttest_holm_detail.csv` and star-summary CSV |

# Reproducibility Notes

- The evaluation code performs 10 repetitions of 10-fold cross-validation, with `Count_CV = 10` and `nCV = 10`.
- The parameter analyses reported in Tables 2 and 3 use the first complete 10-fold cross-validation run (`num = 1`) on Fdataset.
- Parameter selection is based on the complete result obtained after all 10 folds of this fixed run have been evaluated, rather than on the result of a single fold.
- Keep the 10-fold partition, fold assignments, random seed, association masking procedure, WKNN preprocessing procedure, and evaluation order identical across all candidate parameter settings.
- For every fold, mask the held-out associations before WKNN preprocessing and model fitting.
- Compare candidate settings using the run-level `AUC + AUPR`, corresponding to `A_AUC_values(1) + m_A_AUPR_values(1)`, together with consideration of latent-dimensional complexity.
- After parameter selection, keep the selected settings fixed for the reported repeated 10-fold cross-validation and drug-side cold-start evaluations.
- Use identical repeated splits or an identical held-out target-drug order across all compared methods.
- Use distinct MAT filenames for different datasets, methods, evaluation settings, and ablation variants.
- Regenerate the relevant CSV files whenever a source MAT result file is changed.
- Generate the CSV display values directly from the verified MAT result files rather than from manually entered manuscript values.

# A Quickstart Guide
Users can immediately start using multiGMF by running `Demo_multiGMF.m` in MATLAB.
* `Demo_multiGMF.m`: demonstrates final full-data matrix completion on the Fdataset using the multiGMF algorithm. Candidate rankings are subsequently obtained from `M_recovery` after known associations are removed.

# Contact
If you have any questions or suggestions about the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
