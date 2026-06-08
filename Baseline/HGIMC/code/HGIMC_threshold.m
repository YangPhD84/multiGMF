function [T_recovery, iter] = HGIMC_threshold(alpha, beta, HGI_alpha, threshold, Wdd, Wrr, T,  tol1, tol2, maxiter)

trIndex = double(T ~= 0);
[BMC_M,iter] = fBMC(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
TT=BMC_M.*double(BMC_M > threshold);

TTT = TT + T;
TTT(TTT>1) = 1;

T_recovery= fHGI(HGI_alpha,Wdd,Wrr,TTT);

end



