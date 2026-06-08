function Rt=fMBiRW(alpha,l,r,RWrr,RWdd,RWrd,Wrname,Wdname)

A = RWrd;

% % alpha = 0.3;
% % l = 2;
% % r = 2;
d = log(9999);

dn = size(RWdd,1);
dr = size(RWrr,1);

newWrr = sharematrix(A);
newWdd =  sharematrix(A');
LWrr = RWrr;
LWdd = RWdd;
cr = setparFun(RWrd,LWrr);
cd = setparFun(RWrd',LWdd);

LWrr = 1./(1+exp(cr*LWrr+d));
LWdd = 1./(1+exp(cd*LWdd+d));

[RWrr,RWdd] = nManiCluester(LWrr,LWdd,newWrr,newWdd,Wrname,Wdname);

normWrr = normFun(RWrr);
normWdd = normFun(RWdd);

R0 = A/sum(A(:));%？？？？
Rt = R0;

%%
for t=1:max(l,r)
    ftl = 0;
    ftr = 0;
    
    if(t<=l)
        nRtleft = alpha * normWrr*Rt + (1-alpha)*R0;
        ftl = 1;
    end
    if(t<=r)
        nRtright = alpha * Rt * normWdd + (1-alpha)*R0;
        ftr = 1;
    end
    Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
end


