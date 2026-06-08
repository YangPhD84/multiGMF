function [X, Y, i,stop1,stop2]=fMSNMF2_R_mu(M,D,R,lambda1,lambda2,lambda3,lambda33,lambda4,lambda44,k,tol1,tol2,maxiter)
rand('state',2019); % fix random seed
omega=double(M~=0);
omega_=ones(size(omega))-omega;
U=rand(size(M,1),k);
V=rand(size(M,2),k);
P=rand(size(D,2),k);
Q=rand(size(R,2),k);
X=U;
Y=V;
Z=M;
W1=zeros(size(U));
W2=zeros(size(V));
XY=M;

rho = 1.05;
mu = 1e-4;
max_mu = 1e20;

stop1 = 1;
stop2 = 1;
for i=1:maxiter
    
    U=(Z*V+lambda1*D*P-W1+mu*X)*inv(V'*V+lambda1*P'*P+(lambda3+mu)*eye(k));
    
    V=(Z'*U+lambda2*R*Q-W2+mu*Y)*inv(U'*U+lambda2*Q'*Q+(lambda33+mu)*eye(k));
    
    P=(lambda1*D'*U)*inv(lambda1*U'*U+lambda4*eye(k));
    
    Q=(lambda2*R'*V)*inv(lambda2*V'*V+lambda44*eye(k));
    
    X=U+(1/mu)*W1;
    X(X < 0) = 0;
    
    Y=V+(1/mu)*W2;
    Y(Y < 0) = 0;
    
    Z=M.*omega+(U*V').*omega_;
    
    W1=W1+mu*(U-X);
    
    W2=W2+mu*(V-Y);
    
    stop1_0 = stop1;
    stop1=norm(X*Y'-XY,'fro')/norm(XY,'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    XY=X*Y';
    if stop1<tol1 && stop2<tol2
        break
    else
        mu = min(mu * rho, max_mu);
    end
    
end
end