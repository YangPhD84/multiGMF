
Wrr = importdata('DrugSimMat');
Wdd = importdata('DiseaseSimMat');

Wtt = importdata('TargetSimMat.mat');
Wdr = importdata('DiDrAMat');
Wrt = importdata('DrTaAMat.mat');

alpha = 0.1;
thresh_value = 0.5;

for i=1:dn
    for j=1:dn
        if(Wdd(i,j)<thresh_value)
            Wdd(i,j) = 0;
        end
    end
end
for i=1:dr
    for j=1:dr
        if(Wrr(i,j)<thresh_value)
            Wrr(i,j) = 0;
        end
    end
end
for i=1:dt
    for j=1:dt
        if(Wtt(i,j)<thresh_value)
            Wtt(i,j) = 0;
        end
    end
end


Wdr0 = Wdr;
Wrt0 = Wrt;
	
	
Wdrn = Wdr0;
Wrtn = Wrt0;
    
fWdr = Wrr*Wrt0*Wtt*Wrt0';
norm_fWdr = rc_normFun(fWdr);
SWdrn = alpha*Wdrn*norm_fWdr + (1-alpha)*Wdr0;
    
fWrt = Wdr0'*Wdd*Wdr0*Wrr;
norm_fWrt = rc_normFun(fWrt);
SWrtn = alpha*norm_fWrt*Wrtn + (1-alpha)*Wrt0;
    
while(max(max(abs(SWdrn - Wdrn)))>10^-10)
    Wdrn = SWdrn;
    Wrtn = SWrtn;
        
    fWdr = Wrr*Wrtn*Wtt*Wrtn';
    norm_fWdr = rc_normFun(fWdr);
    SWdrn = alpha*Wdrn*norm_fWdr + (1-alpha)*Wdr0;
        
    fWrt = Wdrn'*Wdd*Wdrn*Wrr;
    norm_fWrt = rc_normFun(fWrt);
    SWrtn = alpha*norm_fWrt*Wrtn + (1-alpha)*Wrt0;
end