%%======================== De novo ========================%%
% %注意：1.weight的取值
%        2.denovo for one和de novo for all
%        3.save 名称（方法_方案_数据）
clear
%% 1.Load Datesets
load Fdataset
% load Cdataset
% load CTDdataset2023

Wrr1=drug_ChemS;
Wrr2=drug_AtcS;
Wrr5=drug_TargetS;
Wrr=[Wrr1,Wrr2,Wrr5];

Wdd1=disease_PhS;
Wdd2=disease_DoS;
Wdd=[Wdd1,Wdd2];
Wdr=didr;
Wrd = Wdr';

%% 2.参数赋值                                           2. 参数赋值
lambda1=0.01;
lambda2=lambda1;
lambda3=0.1;
lambda33=lambda3;
lambda4=lambda1;
lambda44=lambda4;
numk=0.7;
tol1=2*1e-3;
tol2=1e-4;
maxiter=300;
[dn,dr] = size(Wdr);
rankk=floor(min(dn,dr)*numk);

NumAS = sum(Wdr);
p_drugPos = find(NumAS==1);%-------------- for one -------------
% % p_drugPos = find(NumAS~=0);%-------------- for all -------------
p_drugLen = size(p_drugPos,2);
DrugResult = zeros(dn,1);
%% ===== 每个 ePos/test drug 的指标记录 =====
RowAucValue = zeros(1,p_drugLen);
RowAuPRValue = zeros(1,p_drugLen);
RowPrecisionValue = zeros(1,p_drugLen);
Row_R_m_A_AUPR_value = zeros(1,p_drugLen);

n_RowAucValue = zeros(1,p_drugLen);
n_RowAuPRValue = zeros(1,p_drugLen);
n_RowPrecisionValue = zeros(1,p_drugLen);

A_DresultMat_TPR = zeros(p_drugLen,dn);
A_DresultMat_FPR = zeros(p_drugLen,dn);
A_DresultMat_Pre = zeros(p_drugLen,dn);

iter_ePos = zeros(1,p_drugLen);           % 每个ePos/test drug的迭代次数
TIME_ePos = zeros(1,p_drugLen);           % 每个ePos/test drug运行时间
t_all = tic;                              % Denovo总运行时间
for num = 1:p_drugLen
    num
    t_epos = tic;
    test_r_index = p_drugPos(num);
    ePos = find(Wdr(:,test_r_index)==1);%-------------- for one -------------
    % %     ePos = find(Wdr(:,test_r_index)~=0);%-------------- for all -------------
    Tfnum = length(ePos);
    
    Wdr(ePos,test_r_index)= 0;%行列坐标
    
    P_TMat = Wdr;
    
    %% MSNMF
    [U,V,iter,stop1,stop2]=fMSNMF2_R_mu(P_TMat,Wdd,Wrr,lambda1,lambda2,lambda3,lambda33,lambda4,lambda44,rankk,tol1,tol2,maxiter);
    M_ResultMat=U*V';
    iter_ePos(num) = iter;
    
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
    
    %% ===== 当前 ePos/test drug 的 AUC、AUPR、Precision =====
    RowAucValue(num) = trapz(FPRArray,TPRArray);
    RowAuPRValue(num) = trapz(TPRArray,PrecisionArray);
    RowPrecisionValue(num) = PrecisionArray(1);
    
    % 带起点修正版本
    n_RowAucValue(num) = trapz(n_FPRArray,n_TPRArray);
    n_RowAuPRValue(num) = trapz(n_TPRArray,n_PrecisionArray);
    n_RowPrecisionValue(num) = n_PrecisionArray(2);  % 等价于 PrecisionArray(1)

     % 当前药物单独计算 R_m_A_AUPR_value
    Row_R_m_TPRArray = [0, TPRArray];
    Row_R_m_PreArray = [1, PrecisionArray];
    Row_R_m_A_AUPR_value(num) = trapz(Row_R_m_TPRArray, Row_R_m_PreArray);
    
    A_DresultMat_TPR(num,:) = TPRArray;
    A_DresultMat_FPR(num,:) = FPRArray;
    A_DresultMat_Pre(num,:) = PrecisionArray;
    
    TIME_ePos(num) = toc(t_epos);
    
    Wdr(ePos,test_r_index)= 1;
    
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

%% ===== 新增：Denovo ePos-level mean ± SD 统计 =====
Denovo_AUC_mean = mean(RowAucValue);
Denovo_AUC_SD   = std(RowAucValue);

% AUPR建议使用带起点修正的 n_RowAuPRValue
Denovo_AUPR_mean = mean(Row_R_m_A_AUPR_value);
Denovo_AUPR_SD   = std(Row_R_m_A_AUPR_value);

Denovo_Precision_mean = mean(RowPrecisionValue);
Denovo_Precision_SD   = std(RowPrecisionValue);

% 同时保留未修正AUPR
Denovo_raw_AUPR_mean = mean(RowAuPRValue);
Denovo_raw_AUPR_SD   = std(RowAuPRValue);

Denovo_iter_mean = mean(iter_ePos);
Denovo_iter_SD   = std(iter_ePos);

DenovoStats.AUC_mean = Denovo_AUC_mean;
DenovoStats.AUC_SD = Denovo_AUC_SD;
DenovoStats.AUPR_mean = Denovo_AUPR_mean;
DenovoStats.AUPR_SD = Denovo_AUPR_SD;
DenovoStats.Row_R_m_A_AUPR_value = Row_R_m_A_AUPR_value; % 新增
DenovoStats.Precision_mean = Denovo_Precision_mean;
DenovoStats.Precision_SD = Denovo_Precision_SD;
DenovoStats.raw_AUPR_mean = Denovo_raw_AUPR_mean;
DenovoStats.raw_AUPR_SD = Denovo_raw_AUPR_SD;
DenovoStats.iter_mean = Denovo_iter_mean;
DenovoStats.iter_SD = Denovo_iter_SD;

% 小数点后4位
DenovoStats.AUC_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUC_mean, Denovo_AUC_SD);
DenovoStats.AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUPR_mean, Denovo_AUPR_SD);
DenovoStats.Precision_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_Precision_mean, Denovo_Precision_SD);
DenovoStats.raw_AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_raw_AUPR_mean, Denovo_raw_AUPR_SD);
DenovoStats.iter_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_iter_mean, Denovo_iter_SD);

runningtime = toc(t_all);

parameters.lambda1 = lambda1;
parameters.lambda2 = lambda2;
parameters.lambda3 = lambda3;
parameters.lambda33 = lambda33;
parameters.lambda4 = lambda4;
parameters.lambda44 = lambda44;
parameters.numk = numk;
parameters.rankk = rankk;
parameters.maxiter = maxiter;
parameters.iter_last = iter;
parameters.iter_mean = Denovo_iter_mean;
parameters.iter_SD = Denovo_iter_SD;
parameters.tol1 = tol1;
parameters.tol2 = tol2;
fprintf('Curve-mean result: AUC=%4.4f, AUPR=%4.4f, Precision=%4.4f.\n', ...
    Result_AUC_value, R_m_A_AUPR_value, Result_Precision_value);

fprintf('ePos-level mean ± SD: AUC=%s, AUPR=%s, Precision=%s.\n', ...
    DenovoStats.AUC_mean_SD_text, ...
    DenovoStats.AUPR_mean_SD_text, ...
    DenovoStats.Precision_mean_SD_text);

%% 4 测试结果保存和记录

%根据不同数据集替换 Fdataset  Cdataset CTDdataset2023
save Fdataset_STresult_MSBMF_Denovoone ...
    RowAucValue RowAuPRValue RowPrecisionValue ...
    n_RowAucValue n_RowAuPRValue n_RowPrecisionValue ...
    Row_R_m_A_AUPR_value ...
    Denovo_AUC_mean Denovo_AUC_SD ...
    Denovo_AUPR_mean Denovo_AUPR_SD ...
    Denovo_raw_AUPR_mean Denovo_raw_AUPR_SD ...
    Denovo_Precision_mean Denovo_Precision_SD ...
    Denovo_iter_mean Denovo_iter_SD ...
    DenovoStats TIME_ePos iter_ePos runningtime ...
    lambda1 lambda2 lambda3 lambda33 lambda4 lambda44 numk rankk tol1 tol2 maxiter