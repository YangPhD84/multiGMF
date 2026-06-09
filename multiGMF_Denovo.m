%%======================== De novo ========================%%
% %注意：1.weight的取值
%        2.denovo for one和de novo for all
%        3.save 名称（方法_方案_数据）
clear
addpath('code');
rng('default');
%% 1.载入数据                                            
% load Fdataset
load Cdataset 
% load CTDdataset2023

drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

Wdr=didr;
Wrd = Wdr';
tic
[dn,dr] = size(Wdr);

 Wrr=(drug_ChemS+drug_AtcS+drug_TargetS)/3;
 Wdd=(disease_PhS+disease_DoS)/2;

%% 2.参数赋值                                           2. 参数赋值 


%% Fdataset  knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.0001 lambda3= 1
% Options.MaxIter=300; %迭代次数
% lambda_soft= 1;
% lambda1= 0.0001;%0.0001
% lambda3= 1;
% tau = 0.7;
% KNN_K = 10;%10

%% Cdataset  knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.0001 lambda3= 1
Options.MaxIter=300; 
lambda_soft= 1;
lambda1= 0;%0.0001
lambda3= 1;
tau = 0.7;
KNN_K = 10;

%% CTDdataset2023 knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.01 lambda3= 1
% Options.MaxIter=300; 
% lambda_soft= 1;
% lambda1= 0.01;%0.01
% lambda3= 1;
% tau = 0.7;
% KNN_K = 10;

out_t = [];%记录每次迭代终止的次数

Options.lambda_soft=lambda_soft;
Options.lambda1=lambda1;
Options.lambda3=lambda3;
Options.mu1 = 1;
Options.mu2 = 1;

tol1 = 2*1e-3;
tol2 = 1*1e-4;
Options.tol1 = tol1;%5*1e-3  2*1e-3
Options.tol2 = tol2;%2*1e-3  1*1e-4

[dn,dr] = size(Wdr);

NumAS = sum(Wdr);
p_drugPos = find(NumAS==1);%-------------- for one -------------
% % p_drugPos = find(NumAS~=0);%-------------- for all -------------
p_drugLen = size(p_drugPos,2);
DrugResult = zeros(dn,1);
%% ===== 每个 ePos/test drug 的指标记录 =====
RowAucValue = zeros(1,p_drugLen);
RowAuPRValue = zeros(1,p_drugLen);
RowPrecisionValue = zeros(1,p_drugLen);   % 新增：每个 ePos/test drug 的 Precision

n_RowAucValue = zeros(1,p_drugLen);
n_RowAuPRValue = zeros(1,p_drugLen);
n_RowPrecisionValue = zeros(1,p_drugLen); % 新增：带起点修正版本的 Precision

A_DresultMat_TPR = zeros(p_drugLen,dn);
A_DresultMat_FPR = zeros(p_drugLen,dn);
A_DresultMat_Pre = zeros(p_drugLen,dn);

TIME_ePos = zeros(1,p_drugLen);           % 新增：每个 ePos/test drug 的运行时间
t_all = tic;                              % 新增：Denovo 总运行时间

for num = 1:p_drugLen
    num
    t_epos = tic;   % 新增：记录当前 ePos/test drug 的运行时间
    test_r_index = p_drugPos(num);
    ePos = find(Wdr(:,test_r_index)==1);%-------------- for one -------------
    % %     ePos = find(Wdr(:,test_r_index)~=0);%-------------- for all -------------
    Tfnum = length(ePos);
    
    Wdr(ePos,test_r_index)= 0;%行列坐标
    
    P_TMat = Wdr;
    

%% fmultiGrLapNMF
Input.A={drug_ChemS,drug_AtcS,drug_TargetS};
Input.B={disease_PhS,disease_DoS};

Input.Wdr = Wdr;

%低秩矩阵初始化

min_mn = min(dn,dr);
rankk = floor(min_mn*tau);
Input.HInit = rand(dn,rankk);
Input.WInit = rand(dr,rankk);

%WKNN处理

Input.kk = KNN_K;

Input.X = Wdr;
Input.Wrr = Wrr;
Input.Wdd = Wdd;
Input.P_TMat = P_TMat';

% % 单相似性对比 R1\R2\R5
% Input.A = {drug_TargetS};%{drug_ChemS,drug_AtcS,drug_SideS,drug_DDIS,drug_TargetS}
% Input.B = {disease_PhS};%{disease_PhS,disease_DoS}
% % % 
% Input.Wrr = drug_TargetS;%{drug_ChemS,drug_AtcS,drug_SideS,drug_DDIS,drug_TargetS}
% Input.Wdd = disease_PhS;%{disease_PhS,disease_DoS}

%执行迭代函数
Output = fmultiGMF(Input, Options);
% Output = no_soft_fmultiGrLapNMF(Input, Options);

%矩阵填充复原
WW=Output.W*Output.H';
M_ResultMat=WW';   

    
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

    A_DresultMat_TPR(num,:) = TPRArray;
    A_DresultMat_FPR(num,:) = FPRArray;
    A_DresultMat_Pre(num,:) = PrecisionArray;

    TIME_ePos(num) = toc(t_epos);   % 新增：当前 ePos/test drug 的运行时间

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
% 说明：
% RowAucValue / RowAuPRValue / RowPrecisionValue 是每一个 ePos/test drug 的独立评价结果。
% 因此这里统计的是所有 ePos/test drug 结果的 mean ± SD。

Denovo_AUC_mean = mean(RowAucValue);
Denovo_AUC_SD   = std(RowAucValue);

% AUPR 建议使用 n_RowAuPRValue，因为它加入了起点 [0,1] 修正
Denovo_AUPR_mean = R_m_A_AUPR_value;
Denovo_AUPR_SD   = std(Row_R_m_A_AUPR_value);

Denovo_Precision_mean = mean(RowPrecisionValue);
Denovo_Precision_SD   = std(RowPrecisionValue);

% 同时保留未修正 AUPR 的 mean ± SD，方便你后续对比
Denovo_raw_AUPR_mean = mean(RowAuPRValue);
Denovo_raw_AUPR_SD   = std(RowAuPRValue);

% 保存为结构体，后续 load .mat 后读取更方便
DenovoStats.AUC_mean = Denovo_AUC_mean;
DenovoStats.AUC_SD = Denovo_AUC_SD;
DenovoStats.AUPR_mean = Denovo_AUPR_mean;
DenovoStats.AUPR_SD = Denovo_AUPR_SD;
DenovoStats.Precision_mean = Denovo_Precision_mean;
DenovoStats.Precision_SD = Denovo_Precision_SD;

DenovoStats.raw_AUPR_mean = Denovo_raw_AUPR_mean;
DenovoStats.raw_AUPR_SD = Denovo_raw_AUPR_SD;

% 字符串形式：默认小数点后4位，便于直接复制到论文表格
DenovoStats.AUC_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUC_mean, Denovo_AUC_SD);
DenovoStats.AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUPR_mean, Denovo_AUPR_SD);
DenovoStats.Precision_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_Precision_mean, Denovo_Precision_SD);
DenovoStats.raw_AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_raw_AUPR_mean, Denovo_raw_AUPR_SD);

parameters.lambda_soft = lambda_soft;
parameters.lambda1 = lambda1;
parameters.lambda3 = lambda3;
parameters.tau = tau;
parameters.KNN_K = KNN_K;
parameters.tol1 = tol1;
parameters.tol2 = tol2;
parameters.MaxIter = Options.MaxIter;

fprintf('Curve-mean result: AUC=%4.4f, Precision=%4.4f, AUPR=%4.4f.\n', ...
    Result_AUC_value, Result_Precision_value, R_m_A_AUPR_value);

fprintf('ePos-level mean ± SD: AUC=%s, AUPR=%s, Precision=%s.\n', ...
    DenovoStats.AUC_mean_SD_text, ...
    DenovoStats.AUPR_mean_SD_text, ...
    DenovoStats.Precision_mean_SD_text);

% 不建议继续使用 save matlab，避免额外生成 matlab.mat 并保存大量无关变量
% save matlab
%% 4 测试结果保存和记录

save Cdataset_no_GrLap_STresult_multiGrNMF_Denovoone ...
    p_drugLen p_drugPos parameters ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value ...
    R_m_TPRArray R_m_PreArray R_m_A_AUPR_value ...
    RowAucValue RowAuPRValue RowPrecisionValue ...
    n_RowAucValue n_RowAuPRValue n_RowPrecisionValue ...
    Denovo_AUC_mean Denovo_AUC_SD ...
    Denovo_AUPR_mean Denovo_AUPR_SD ...
    Denovo_raw_AUPR_mean Denovo_raw_AUPR_SD ...
    Denovo_Precision_mean Denovo_Precision_SD ...
    DenovoStats TIME_ePos runningtime ...
    lambda_soft lambda1 lambda3 tau KNN_K tol1 tol2
