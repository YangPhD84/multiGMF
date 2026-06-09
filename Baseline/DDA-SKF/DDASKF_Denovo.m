%%======================== De novo ========================%%
% %注意：1.weight的取值
%        2.denovo for one和de novo for all
%        3.save 名称（方法_方案_数据）
clear
addpath('code');
rng('default');
%% 1.载入数据                                            
% load Fdataset;

% load CTDdataset2023
load Cdataset;


drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

Input.drugChemS = drug_ChemS;
% Input.drugDDIS = drug_DDIS;
Input.drugAtcS = drug_AtcS;
% Input.drugSideS = drug_SideS;
Input.drugTargetS = drug_TargetS;
Input.diseaseDoS = disease_DoS;
Input.diseasePhS = disease_PhS;
Wdr = didr;%didr


Input.beta = 0.4;
Input.lamuda = 2^(-16);
% Input.predictSimMatdcChemical = predictSimMatdcChemical;
% Input.predictSimMatdg = predictSimMatdg;
% Input.predictSimMatdcGo = predictSimMatdcGo;
% Input.predictSimMatdcDomain = predictSimMatdcDomain;
% Wdr = predictAdMatdgc;%didr

%C数据集测试
% Input.predictSimMatdcChemical = C_Simmatdc_chemical;
% Input.predictSimMatdg = C_Simmatdc_dg;
% Input.predictSimMatdcGo = C_Simmatdc_Go;
% Input.predictSimMatdcDomain = C_Simmatdc_domain;
% Wdr = C_mat_dgc;%didr

[dn,dr] = size(Wdr);
Wrd = Wdr';

NumAS = sum(Wdr);
p_drugPos = find(NumAS==1);%-------------- for one -------------
% % p_drugPos = find(NumAS~=0);%-------------- for all -------------
p_drugLen = size(p_drugPos,2);
DrugResult = zeros(dn,1);
%% complete ten-fold cross validation ten times;
%% ===== 每个 ePos/test drug 的指标记录 =====
RowAucValue = zeros(1,p_drugLen);
RowAuPRValue = zeros(1,p_drugLen);
RowPrecisionValue = zeros(1,p_drugLen);   % 新增：每个 ePos/test drug 的 Precision
Row_R_m_A_AUPR_value = zeros(1,p_drugLen);

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
    

%% DDA-SKF
score_matrix = DDASKF(Input, P_TMat);

M_ResultMat = score_matrix;

    
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
Row_R_m_A_AUPR_value(num) = R_m_A_AUPR_value; % 记录每个目标药物的 R_m_A_AUPR_value

%% ===== 新增：Denovo ePos-level mean ± SD 统计 =====
% RowAucValue / RowAuPRValue / RowPrecisionValue：
% 每一个元素对应一个 ePos/test drug 的独立评价结果。

Denovo_AUC_mean = mean(RowAucValue);
Denovo_AUC_SD   = std(RowAucValue);

% AUPR建议使用带起点修正的 n_RowAuPRValue
Denovo_AUPR_mean = R_m_A_AUPR_value;
Denovo_AUPR_SD   = std(Row_R_m_A_AUPR_value);

Denovo_Precision_mean = mean(RowPrecisionValue);
Denovo_Precision_SD   = std(RowPrecisionValue);

% 同时保留未修正AUPR，方便对比
Denovo_raw_AUPR_mean = mean(RowAuPRValue);
Denovo_raw_AUPR_SD   = std(RowAuPRValue);

DenovoStats.AUC_mean = Denovo_AUC_mean;
DenovoStats.AUPR_mean = Denovo_AUPR_mean;
DenovoStats.AUPR_SD = Denovo_AUPR_SD;
DenovoStats.Row_R_m_A_AUPR_value = Row_R_m_A_AUPR_value; % 新增，每个药物的 AUPR
DenovoStats.AUPR_SD = Denovo_AUPR_SD;
DenovoStats.Precision_mean = Denovo_Precision_mean;
DenovoStats.Precision_SD = Denovo_Precision_SD;

DenovoStats.raw_AUPR_mean = Denovo_raw_AUPR_mean;
DenovoStats.raw_AUPR_SD = Denovo_raw_AUPR_SD;

% 字符串形式，小数点后4位，便于论文表格复制
DenovoStats.AUC_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUC_mean, Denovo_AUC_SD);
DenovoStats.AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_AUPR_mean, Denovo_AUPR_SD);
DenovoStats.Precision_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_Precision_mean, Denovo_Precision_SD);
DenovoStats.raw_AUPR_mean_SD_text = sprintf('%.4f ± %.4f', Denovo_raw_AUPR_mean, Denovo_raw_AUPR_SD);

runningtime = toc(t_all);

% 不建议继续使用 save matlab，避免额外生成 matlab.mat 并保存大量无关变量
% save matlab

%% 4 测试结果保存和记录

%根据不同数据集替换 Fdataset  Cdataset CTDdataset2023
save Cdataset_STresult_DDASKF_Denovoone ...
    p_drugLen p_drugPos ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value ...
    R_m_TPRArray R_m_PreArray R_m_A_AUPR_value ...
    RowAucValue RowAuPRValue RowPrecisionValue ...
    n_RowAucValue n_RowAuPRValue n_RowPrecisionValue ...
    Denovo_AUC_mean Denovo_AUC_SD ...
    Denovo_AUPR_mean Denovo_AUPR_SD ...
    Denovo_raw_AUPR_mean Denovo_raw_AUPR_SD ...
    Denovo_Precision_mean Denovo_Precision_SD ...
    DenovoStats TIME_ePos runningtime
