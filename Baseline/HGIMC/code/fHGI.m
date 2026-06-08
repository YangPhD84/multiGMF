function [Result_Wdr] = fHGI(alpha,Wdd,Wrr,Wdr)

normWdd = fNorm(Wdd);
normWrr = fNorm(Wrr);
Wdr0=Wdr;
Wdr_i=Wdr0;
Wdr_I = alpha*normWdd*Wdr_i*normWrr + (1-alpha)*Wdr0;

while(max(max(abs(Wdr_I - Wdr_i)))>10^-10)
    Wdr_i=Wdr_I;
    Wdr_I=alpha*normWdd*Wdr_i*normWrr + (1-alpha)*Wdr0;
end

Result_Wdr = Wdr_I;

end
