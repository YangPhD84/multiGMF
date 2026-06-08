function [ P_ResultMat ] = Fun_Methods(alpha,fWrr,fWdd,fWdr)
%     alpha = 0.4;
    normWrr = myFunNorm(fWrr);
    normWdd = myFunNorm(fWdd);

    Wrd= fWdr';

    Wrd0 = Wrd;
    Wrdn = Wrd0;
    SWrdn = alpha*normWrr*Wrdn*normWdd + (1-alpha)*Wrd0;

    while(max(max(abs(SWrdn - Wrdn)))>10^-10)
        Wrdn = SWrdn;
        SWrdn = alpha*normWrr*Wrdn*normWdd + (1-alpha)*Wrd0;
    end

    Wdr_t = SWrdn';
    P_ResultMat = Wdr_t;
end
