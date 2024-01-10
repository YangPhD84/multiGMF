function [W, H, iter] = multiGMF(M, R, D, tau, Wrr, Wdd, r, k, MaxIter, lambda1, lambda2, lambda3, tol1, tol2)
% MSGrNMF: Drug repositioning based on multi-similarity graph regularised non-negative matrix decomposition
% Usage: [U, V, iter] = MSGrNMF(M, R, D, W, H, Wrr, Wdd, r, k, MaxIter, lambda1, lambda2, lambda3, tol1, tol2)
%
% Inputs:
%        M                  - the target matrix with only known entries and the unobserved entries are 0.
%        R                  - drug similarity matrix
%        D                  - disease similarity matrix
%        tau                - the latent dimension of matrix factorization.
%        Wrr                - mean similarity matrix of drugs.
%        Wdd                - mean similarity matrix of diseases.
%        r                  - decay term for weighted k-nearest neighbours.
%        k                  - neighborhood size of weighted k-nearest neighbours
%        MaxIter            - maximum number of iterations.
%        lambda1,2,3        - parameters needed to give.
%        tol1, tol2         - tolerance of termination conditions.
%
% Outputs:
%        W, H               - two latent low-rank matrices of the completed matrix.
%        iter               - the number of iterations.
%
% Written by: Bin Yang
% Email: yangbin_9966@163.com
% Created: August 12, 2023

X = WKNN( M, Wrr, Wdd, k, r );

WH = X;
[dr, dn] = size(X);
min_mn = min(dn, dr);
rankk = floor(min_mn*tau);
H = rand(dn,rankk);
W = rand(dr,rankk);
stop1 = 1;
stop2 = 1;

ER = eye(size(W, 1));
ED = eye(size(H, 1));

GraphNum1 = size(R, 2);
GraphNum2 = size(D, 2);
Rweight = ones(1, GraphNum1) / GraphNum1;
Dweight = ones(1, GraphNum2) / GraphNum2;

for k_r = 1 : GraphNum1
    Lr{k_r} = diag(sum(R{k_r},2)) - R{k_r};
end

for k_d = 1 : GraphNum2
    Ld{k_d} = diag(sum(D{k_d},2)) - D{k_d};
end


for iter = 1 : MaxIter
    
    for k_r = 1 : GraphNum1
        R{k_r} = inv(Rweight(k_r) * ER + lambda1 * Lr{k_r}) * (Rweight(k_r) * W);
    end
    
    for k_d = 1 : GraphNum2
        D{k_d} = inv(Dweight(k_d) * ED + lambda2 * Ld{k_d}) * (Dweight(k_d) * H);
    end
    
    Sr = zeros(size(R{1}));
    for k_r = 1 : GraphNum1
        Sr = Sr + Rweight(k_r) * R{k_r};
    end
    
    Sw = zeros(size(W));
    for k_r = 1 : GraphNum1
        Sw = Sw + Rweight(k_r) * W;
    end
    
    Sd=zeros(size(D{1}));
    for k_d = 1 : GraphNum2
        Sd = Sd + Dweight(k_d) * D{k_d};
    end
    
    Sh = zeros(size(H));
    for k_d = 1 : GraphNum2
        Sh = Sh + Dweight(k_d) * H;
    end
    
    W = (X * H + Sr) ./ (W * H' * H + Sw + lambda3 * W) .*  W + 1e-8;
    H = (X'* W + Sd) ./ (H * W' * W + Sh + lambda3 * H) .*  H + 1e-8;
    
    stop1_0 = stop1;
    stop1 = norm(W * H' - WH, 'fro')/norm(WH, 'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    WH = W * H';
    
    if stop1 < tol1 && stop2 < tol2
        break
    end
end
end
