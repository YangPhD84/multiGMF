clear all
addpath('Datasets');
rand('state', 2023); % fix random seed
%% 1. Load Datesets
load  Gold_standard_dataset
% load CTDdataset2023
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

%% 2. multiGMF algorithm
lambda1 = 0.0001;
lambda2 = lambda1;
lambda3 = 1;
r = 0.9;
k = 10;
tau = 0.7;
MaxIter = 300;
tol1 = 2*1e-3;
tol2 = 1*1e-4;
[W, H, iter] = multiGMF(Wrd,R,D,tau,Wrr,Wdd,r,k,MaxIter,lambda1,lambda2,lambda3,tol1,tol2);
M_recovery = W * H';