clear all
addpath('Datasets');
rand('state', 2023); % fix random seed
%% 1. Load Datesets
load  Fdataset
% load Cdataset
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
Wdr = didr; % disease-by-drug association matrix in the dataset
Wrd = Wdr';% drug-by-disease association matrix used by fmultiGMF
[dn, dr] = size(Wdr);
min_mn = min(dn, dr);

%% 2. multiGMF algorithm
lambda1 = 0.0001;
lambda2 = lambda1;   % 保留变量，方便记录；新版 fmultiGMF 中主要使用 lambda_soft、lambda1、lambda3
lambda3 = 1;
r = 0.9;             % WKNN 衰减系数已在新版 fmultiGMF 内部固定为 0.9，这里保留用于参数记录
k = 10;
tau = 0.7;
MaxIter = 300;
tol1 = 2*1e-3;
tol2 = 1*1e-4;

%% ==== fmultiGMF 所需的 Input 和 Options =====

rankk = floor(min_mn * tau);

% Input.A 是 drug-side 多个相似性矩阵
Input.A = R;

% Input.B 是 disease-side 多个相似性矩阵
Input.B = D;

% Input.X 使用 disease-by-drug 关联矩阵
Input.X = Wdr;

% 初始化低秩矩阵
% WInit 对应 drug-by-rank
% HInit 对应 disease-by-rank
Input.WInit = rand(dr, rankk);
Input.HInit = rand(dn, rankk);

% WKNN 相关输入
Input.kk = k;
Input.Wrr = Wrr;
Input.Wdd = Wdd;
Input.P_TMat = Wrd;

% Options 参数
Options.MaxIter = MaxIter;
Options.lambda_soft = 1;
Options.lambda1 = lambda1;
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