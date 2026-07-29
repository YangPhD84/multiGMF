clear all
addpath('Datasets');
rand('state', 2023); % fix random seed
%% 1. Load Datesets
load  Fdataset
% load Cdataset
% load CTDdataset2023

drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

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
Wdr = didr; % disease-by-drug association matrix in the dataset
Wrd = Wdr';% drug-by-disease association matrix used by fmultiGMF
[dn, dr] = size(Wdr);
min_mn = min(dn, dr);

%% 2. multiGMF algorithm
lambda1 = 0.0001;
lambda2 = lambda1;   %
lambda3 = 1;
r = 0.9;             %
k = 10;
tau = 0.7;
MaxIter = 300;
tol1 = 2*1e-3;
tol2 = 1*1e-4;

%% ==== fmultiGMF =====

rankk = floor(min_mn * tau);

% Input.A
Input.A = R;

% Input.B
Input.B = D;

% Input.X
Input.X = Wdr;

% 
Input.WInit = rand(dr, rankk);
Input.HInit = rand(dn, rankk);

% WKNN 
Input.kk = k;
Input.Wrr = Wrr;
Input.Wdd = Wdd;
Input.P_TMat = Wrd;

% Options 
Options.MaxIter = MaxIter;
Options.lambda_soft = 1;
Options.lambda1 = lambda1;
Options.lambda2 = lambda2;
Options.lambda3 = lambda3;
Options.mu1 = 1;
Options.mu2 = 1;
Options.tol1 = tol1;
Options.tol2 = tol2;

%% ===== fmultiGMF =====
Output = fmultiGMF(Input, Options);

W = Output.W;
H = Output.H;
iter = Output.t;

M_recovery = W * H';
