function [X,W,iter,stop1,stop2] = fBNNR_econ(alpha,beta,T,trIndex,tol1,tol2,maxiter,a,b)
%This algorithm is called BNNR.
%    Description
%        alpha,beta  - parameters needed to give.
%        T           - the target matrix with only known entries when the unobserved entries are 0.
%        trIndex     - a matrix recording the observed positions in the target matrix.
%        tol1,tol2   - tolerance of termination conditions.
%        maxiter     - maximum number of iterations.
%
% Written by: Mengyun Yang
% Email: mengyunyang@csu.edu.cn
% Created: January 16, 2019

% % a=-10^10;b=10^10;  %No bounded
% % a=0;b=1;  %Bounded
[n1,n2]=size(T);
X=T;
W=X;
Y=X;

i=1;
stop1=1;
stop2=1;
while(stop1>tol1 || stop2>tol2 )
    %the process of getting W
    tran=(1/beta)*(Y+alpha*(T.*trIndex))+X;
    W=tran-(alpha/(alpha+beta))*(tran.*trIndex);
    W(W<a)=a;
    W(W>b)=b;
    
    %the process of getting X
    X_1=mysvt_econ(W-1/beta*Y,1/beta);
    
    %the process of getting Y
    Y=Y+1*beta*(X_1-W);
    
    stop1_0=stop1;
    stop1=norm(X_1-X,'fro')/norm(X,'fro');
    stop2=abs(stop1-stop1_0)/max(1,abs(stop1_0));
    
    X=X_1;
    i=i+1;
    
    if i<maxiter
        iter=i-1;
    else
        iter=maxiter;
        warning('reach maximum iteration~~do not converge!!!');
        break
    end
    
end


