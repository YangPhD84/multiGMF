%%======================== De novo ========================%%

clear

addpath('code');
%% 1.Load Datesets       
load Fdataset;
% load Cdataset;
% load CTDdataset2023

sigma=0.5;
Wrr=(drug_ChemS+drug_AtcS+drug_TargetS)/3;
Wdd=(disease_PhS+disease_DoS)/2;
Wrr=fRBFkernel([Wrr],sigma);
Wdd=fRBFkernel([Wdd],sigma);
 weight=1;%----------- weight ----------------
 Wdr=didr*weight;
 Wrd = Wdr';

maxiter = 300; tol1 = 2*1e-3;   tol2 = 1*1e-5;
alpha=1; beta=10; HGBI_alpha=0.1; threshold=0.1;

[dn,dr] = size(Wdr);
NumAS = sum(Wdr);
p_drugPos = find(NumAS==weight);%-------------- for one -------------
% % p_drugPos = find(NumAS~=0);%-------------- for all -------------
p_drugLen = size(p_drugPos,2);
DrugResult = zeros(dn,1);
RowAucValue = zeros(1,p_drugLen);
RowAuPRValue = zeros(1,p_drugLen);
RowPrecisionValue = zeros(1,p_drugLen);   % 
Row_R_m_A_AUPR_value = zeros(1,p_drugLen);

n_RowAucValue = zeros(1,p_drugLen);
n_RowAuPRValue = zeros(1,p_drugLen);
n_RowPrecisionValue = zeros(1,p_drugLen); %

A_DresultMat_TPR = zeros(p_drugLen,dn);
A_DresultMat_FPR = zeros(p_drugLen,dn);
A_DresultMat_Pre = zeros(p_drugLen,dn);

TIME_ePos = zeros(1,p_drugLen);           %
t_all = tic;                              % 

for num = 1:p_drugLen
    num
    t_epos = tic;   %
    test_r_index = p_drugPos(num);
    ePos = find(Wdr(:,test_r_index)==weight);%-------------- for one -------------
% %     ePos = find(Wdr(:,test_r_index)~=0);%-------------- for all -------------
    Tfnum = length(ePos);

    Wdr(ePos,test_r_index)= 0;%
		
    P_TMat = Wdr;
    
M_ResultMat=HGIMC_threshold(alpha, beta, HGBI_alpha, threshold, Wdd, Wrr, P_TMat,  tol1, tol2, maxiter);       

    DrugResult = M_ResultMat(:,test_r_index);
    Qvalue = M_ResultMat(ePos,test_r_index);
    
    SdrugResult = sort(DrugResult,'descend');
    rds_len = size(ePos,1);
    
    RresultSort = zeros(dn,1);
    
    for k=1:rds_len
        eQvalue = Qvalue(k);
        TfindposMat =  find(SdrugResult==eQvalue);
        TfindposMatlen = size(TfindposMat,1);
        Tfindpos = TfindposMat(TfindposMatlen);
        t_Tfindpos = Tfindpos;
        while(RresultSort(t_Tfindpos)==1)
            t_Tfindpos = t_Tfindpos-1;
        end
        if(t_Tfindpos==0)
            t_Tfindpos = Tfindpos + 1;
        end
        Tfindpos = t_Tfindpos;
        
        RresultSort(Tfindpos) = 1;

    end
    
    TPRArray = zeros(1,dn);
    FPRArray = zeros(1,dn);
    PrecisionArray = zeros(1,dn);

    CountP =  rds_len;
    CountN =  dn - rds_len;

    Tpnum = 0;
    Fpnum = 0;
    
    for m =1:dn
        if(RresultSort(m)==1)
            Tpnum = Tpnum + 1;
        else
            Fpnum = Fpnum + 1;
        end
        TPRArray(m) = Tpnum/CountP;
        FPRArray(m) = Fpnum/CountN;
        PrecisionArray(m) = Tpnum/(Tpnum+Fpnum);
    end
    
    n_TPRArray = [0,TPRArray];
    n_FPRArray = [0,FPRArray];
    n_PrecisionArray = [1,PrecisionArray];
    
    RowAucValue(num) = trapz(FPRArray,TPRArray);
    RowAuPRValue(num) = trapz(TPRArray,PrecisionArray);
    RowPrecisionValue(num) = PrecisionArray(1);
    
    n_RowAucValue(num) = trapz(n_FPRArray,n_TPRArray);
    n_RowAuPRValue(num) = trapz(n_TPRArray,n_PrecisionArray);
    n_RowPrecisionValue(num) = n_PrecisionArray(2);  
    Row_R_m_TPRArray = [0, TPRArray];
    Row_R_m_PreArray = [1, PrecisionArray];
    Row_R_m_A_AUPR_value(num) = trapz(Row_R_m_TPRArray, Row_R_m_PreArray);
    
    A_DresultMat_TPR(num,:) = TPRArray;
    A_DresultMat_FPR(num,:) = FPRArray;
    A_DresultMat_Pre(num,:) = PrecisionArray;
    
    TIME_ePos(num) = toc(t_epos);   

    Wdr(ePos,test_r_index) = weight;
    
end
    runningtime = toc(t_all);
    Result_TPRArray = mean(A_DresultMat_TPR);%VV
    Result_FPRArray = mean(A_DresultMat_FPR);%VV
    Result_PreArray = mean(A_DresultMat_Pre);%VV
    
    Result_AUC_value = trapz(Result_FPRArray,Result_TPRArray);%VV
    Result_AUPR_value = trapz(Result_TPRArray,Result_PreArray);%V
    Result_Precision_value = Result_PreArray(1);%VV
    
    R_m_TPRArray = [0,Result_TPRArray];%V
    R_m_PreArray = [1,Result_PreArray];%V
    R_m_A_AUPR_value = trapz(R_m_TPRArray,R_m_PreArray);%V
    
    Denovo_AUC_mean = mean(RowAucValue);
    Denovo_AUC_SD   = std(RowAucValue);
    
    % 
    Denovo_AUPR_mean = mean(n_RowAuPRValue);
    Denovo_AUPR_SD   = std(n_RowAuPRValue);     
    
    Denovo_Precision_mean = mean(RowPrecisionValue);
    Denovo_Precision_SD   = std(RowPrecisionValue);
    
    %
    Denovo_raw_AUPR_mean = mean(RowAuPRValue);
    Denovo_raw_AUPR_SD   = std(RowAuPRValue);
    R_m_A_AUPR_SD = Denovo_AUPR_SD;
    R_m_A_AUPR_mean_from_rows = Denovo_AUPR_mean;
    
    R_m_A_AUPR_mean_SD_text = sprintf('%.4f ± %.4f', ...
        R_m_A_AUPR_value, R_m_A_AUPR_SD);
	Denovo_AUPR_mean = R_m_A_AUPR_value;
    
    DenovoStats.AUC_mean = Denovo_AUC_mean;
    DenovoStats.AUC_SD = Denovo_AUC_SD;
    
    DenovoStats.AUPR_mean = Denovo_AUPR_mean;
    DenovoStats.AUPR_SD = Denovo_AUPR_SD;
    DenovoStats.Row_R_m_A_AUPR_value = Row_R_m_A_AUPR_value;
    
    DenovoStats.Precision_mean = Denovo_Precision_mean;
    DenovoStats.Precision_SD = Denovo_Precision_SD;
    
    DenovoStats.raw_AUPR_mean = Denovo_raw_AUPR_mean;
    DenovoStats.raw_AUPR_SD = Denovo_raw_AUPR_SD;
    
    %
    DenovoStats.AUC_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUC_mean, Denovo_AUC_SD);
    DenovoStats.AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUPR_mean, Denovo_AUPR_SD);
    DenovoStats.Precision_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_Precision_mean, Denovo_Precision_SD);
    DenovoStats.raw_AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_raw_AUPR_mean, Denovo_raw_AUPR_SD);

    DenovoStats.R_m_A_AUPR_value = R_m_A_AUPR_value;
    DenovoStats.R_m_A_AUPR_SD = R_m_A_AUPR_SD;
    DenovoStats.R_m_A_AUPR_mean_from_rows = R_m_A_AUPR_mean_from_rows;
    DenovoStats.R_m_A_AUPR_mean_SD_text = R_m_A_AUPR_mean_SD_text;
    
    runningtime = toc(t_all);
        
    parameters.HGBI_alpha = HGBI_alpha;
    parameters.BNNR_alpha = alpha;
    parameters.BNNR_beta = beta;
    parameters.threshold = threshold;
    parameters.sigma = sigma;
    parameters.weight = weight;
    parameters.maxiter = maxiter;
    parameters.tol1 = tol1;
    parameters.tol2 = tol2;
    
    
save Fdataset_STresult_HGIMC_Denovoone ...
    p_drugLen p_drugPos parameters ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value ...
    R_m_TPRArray R_m_PreArray R_m_A_AUPR_value ...
    RowAucValue RowAuPRValue RowPrecisionValue ...
    n_RowAucValue n_RowAuPRValue n_RowPrecisionValue ...
    Row_R_m_A_AUPR_value ...
    R_m_A_AUPR_SD R_m_A_AUPR_mean_from_rows R_m_A_AUPR_mean_SD_text ...
    Denovo_AUC_mean Denovo_AUC_SD ...
    Denovo_AUPR_mean Denovo_AUPR_SD ...
    Denovo_raw_AUPR_mean Denovo_raw_AUPR_SD ...
    Denovo_Precision_mean Denovo_Precision_SD ...
    DenovoStats TIME_ePos runningtime ...
    alpha beta HGBI_alpha threshold sigma weight tol1 tol2 maxiter

 
