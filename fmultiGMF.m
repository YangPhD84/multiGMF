function Output = fmultiGMF(Input,Options)
% fmultiGMF: full multiGMF model with WKNN preprocessing, graph regularization, and soft coupling.
%Input.A,B,X,W,H,,  drug5 disease2 association matrix 
%Options.lambda_soft lambda1 lambda3 MaxIter

%% 对drug-disease关联矩阵整体进行KNN处理
r = 0.9;%WKNN 衰减系数
KK = Input.kk;

% WKNN is applied to the masked drug-by-disease training matrix.
% Input.Wrr is the drug-by-drug mean similarity matrix.
% Input.Wdd is the disease-by-disease mean similarity matrix.
X=WKNN( Input.P_TMat, Input.Wrr, Input.Wdd, KK, r );

% %% 采用随机初始化低秩矩阵H、W
H=Input.HInit; %初始化的低秩矩阵H
W=Input.WInit; %初始化的低秩矩阵W

M  = Input.X';
omega=double(M~=0);
omega=ones(size(omega))-omega;

WH = X;
stop1 = 1;
stop2 = 1;

ER = eye(size(W,1));
ED = eye(size(H,1));

% Keep the original similarity matrices unchanged.
% Rsim and Dsim are used only to construct fixed graph Laplacians.
Rsim = Input.A;
Dsim = Input.B;

GraphNum1 = size(Rsim, 2);
GraphNum2 = size(Dsim, 2);
Rweight = ones(1,GraphNum1)/GraphNum1;%ones(1,GraphNum1)/GraphNum1
Dweight = ones(1,GraphNum2)/GraphNum2;%ones(1,GraphNum2)/GraphNum2

for k=1:GraphNum1
    Lr{k} = diag(sum(Rsim{k}, 2)) - Rsim{k};
end

for k=1:GraphNum2
    Ld{k} = diag(sum(Dsim{k}, 2)) - Dsim{k};
end

% Auxiliary representations updated during optimization.
% They are separated from Rsim and Dsim to avoid overwriting input similarities.
R_aux = cell(1, GraphNum1);
D_aux = cell(1, GraphNum2);

%%迭代过程
for t=1:Options.MaxIter
    %Update Input.A{k}
    for k = 1:GraphNum1
        % Update drug-side auxiliary representations.
        R_aux{k} = inv(Options.mu1*Options.lambda_soft*Rweight(k)*ER + Options.lambda1*Lr{k})*(Options.mu1*Rweight(k)*W);
    end
    
    %Update Input.B{k}
    for k = 1:GraphNum2
        % Update disease-side auxiliary representations.
        D_aux{k} = inv(Options.mu2*Options.lambda_soft*Dweight(k)*ED + Options.lambda1*Ld{k})*(Options.mu2*Dweight(k)*H);
    
    end
    
    %Update W,H
    Sr=zeros(size(Input.A{1}));
    for k=1:GraphNum1
        Sr = Sr + Rweight(k) * R_aux{k};
    end
    
    % Since sum(Rweight) = 1, Sw is mathematically equal to W.
    Sw=zeros(size(W));
    for k=1:GraphNum1
        Sw=Sw+Rweight(k)*W;
    end
    
    Sd=zeros(size(Input.B{1}));
    for k=1:GraphNum2
        Sd = Sd + Dweight(k) * D_aux{k};
    end

    % Since sum(Dweight) = 1, Sh is mathematically equal to H.
    Sh=zeros(size(H));
    for k=1:GraphNum2
        Sh=Sh+Dweight(k)*H;
    end

    W = (X*H+Options.mu1*Options.lambda_soft*Sr)./(W*H'*H+Options.mu1*Options.lambda_soft*Sw+Options.lambda3*W).*W+1e-8;
    H = (X'*W+Options.mu2*Options.lambda_soft*Sd)./(H*W'*W+Options.mu2*Options.lambda_soft*Sh+Options.lambda3*H).*H+1e-8;
    
    if stop1<Options.tol1 & stop2<Options.tol2
        break
    end
    
end

Output.H=H;
Output.W=W;
Output.t=t;
Output.Rweight=Rweight;
Output.Dweight=Dweight;
end

