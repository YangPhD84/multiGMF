function Output = no_soft_fmultiGMF(Input,Options)
% no_soft_fmultiGMF: no-soft ablation variant of multiGMF. This function corresponds to the multiGMF without soft coupling experiment reported.
%Input.A,B,X,W,H,,  drug5 disease2 association matrix 
%Options.lambda_soft lambda1 lambda3 MaxIter

%% 对drug-disease关联矩阵整体进行KNN处理
r = 0.9;
KK = Input.kk;

% WKNN preprocessing is retained in the no-soft ablation.
% The input matrix is the masked drug-by-disease training matrix.
X=WKNN( Input.P_TMat, Input.Wrr, Input.Wdd, KK, r );%KNN处理后得到关联矩阵

%% 采用随机初始化低秩矩阵H、W
H=Input.HInit; %初始化的低秩矩阵H
W=Input.WInit; %初始化的低秩矩阵W

% [dn,dr] = size(X);
% H=rand(dr,dr);
% W=rand(dn,dn);

M  = Input.X';
omega=double(M~=0);
omega=ones(size(omega))-omega;

WH = X;
stop1 = 1;
stop2 = 1;

ER = eye(size(W,1));
ED = eye(size(H,1));

GraphNum1=size(Input.A,2);
GraphNum2=size(Input.B,2);
Rweight = ones(1,GraphNum1)/GraphNum1;%ones(1,GraphNum1)/GraphNum1
Dweight = ones(1,GraphNum2)/GraphNum2;%ones(1,GraphNum2)/GraphNum2

for k=1:GraphNum1
    Lr{k}=diag(sum(Input.A{k},2))-Input.A{k};
end

for k=1:GraphNum2
    Ld{k}=diag(sum(Input.B{k},2))-Input.B{k};
end

%%迭代过程
for t=1:Options.MaxIter

    %Update W,H
    Sr=zeros(size(Input.A{1}));
    for k=1:GraphNum1
        Sr=Sr+Input.A{k};
    end
    
    Su=zeros(size(Input.A{1}));
    for k=1:GraphNum1
        Su=Su+diag(sum(Input.A{k},2));
    end
    
    Sd=zeros(size(Input.B{1}));
    for k=1:GraphNum2
        Sd=Sd+Input.B{k};
    end
    
    So=zeros(size(Input.B{1}));
    for k=1:GraphNum2
        So=So+diag(sum(Input.B{k},2));
    end
    
    % No-soft ablation update.
    % The soft coupling terms used in full multiGMF are removed.
    % Graph Laplacian regularization and Tikhonov regularization are retained.
    W = (X*H+Options.lambda1*Sr*W)./(W*H'*H+Options.lambda1*Su*W+Options.lambda3*W).*W+1e-8;
    H = (X'*W+Options.lambda1*Sd*H)./(H*W'*W+Options.lambda1*So*H+Options.lambda3*H).*H+1e-8;
    
    
    stop1_0 = stop1;
    stop1=norm(W*H'-WH,'fro')/norm(WH,'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    WH=W*H';
    
    
%      if stop2<Options.tol
%         break
%     end
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

