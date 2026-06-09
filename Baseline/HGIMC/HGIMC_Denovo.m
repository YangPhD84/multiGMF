%%======================== De novo ========================%%
% %注意：1.weight的取值
%        2.denovo for one和de novo for all
%        3.save 名称（方法_方案_数据）
clear
% addpath('Luohm data');
addpath('code');
%% 1.Load Datesets       
% load Fdataset;
% load Cdataset;
load CTDdataset2023

sigma=0.5;
Wrr=(drug_ChemS+drug_AtcS+drug_TargetS)/3;
Wdd=(disease_PhS+disease_DoS)/2;
Wrr=fRBFkernel([Wrr],sigma);
Wdd=fRBFkernel([Wdd],sigma);
 weight=1;%----------- weight ----------------
 Wdr=didr*weight;
 Wrd = Wdr';

%% 2.参数赋值                                           2. 参数赋值 
maxiter = 300; tol1 = 2*1e-3;   tol2 = 1*1e-5;
alpha=1; beta=10; HGBI_alpha=0.1; threshold=0.1;

[dn,dr] = size(Wdr);
NumAS = sum(Wdr);
p_drugPos = find(NumAS==weight);%-------------- for one -------------
% % p_drugPos = find(NumAS~=0);%-------------- for all -------------
p_drugLen = size(p_drugPos,2);
DrugResult = zeros(dn,1);
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
    ePos = find(Wdr(:,test_r_index)==weight);%-------------- for one -------------
% %     ePos = find(Wdr(:,test_r_index)~=0);%-------------- for all -------------
    Tfnum = length(ePos);

    Wdr(ePos,test_r_index)= 0;%行列坐标
		
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
    
    %% ===== 当前 ePos/test drug 的 AUC、AUPR、Precision =====
    RowAucValue(num) = trapz(FPRArray,TPRArray);
    RowAuPRValue(num) = trapz(TPRArray,PrecisionArray);
    RowPrecisionValue(num) = PrecisionArray(1);
    
    % 带起点修正版本
    n_RowAucValue(num) = trapz(n_FPRArray,n_TPRArray);
    n_RowAuPRValue(num) = trapz(n_TPRArray,n_PrecisionArray);
    n_RowPrecisionValue(num) = n_PrecisionArray(2);  % 等价于 PrecisionArray(1)

    % 新增：当前目标药物单独按照 R_m_A_AUPR_value 的计算逻辑得到的 AUPR
    Row_R_m_TPRArray = [0, TPRArray];
    Row_R_m_PreArray = [1, PrecisionArray];
    Row_R_m_A_AUPR_value(num) = trapz(Row_R_m_TPRArray, Row_R_m_PreArray);
    
    A_DresultMat_TPR(num,:) = TPRArray;
    A_DresultMat_FPR(num,:) = FPRArray;
    A_DresultMat_Pre(num,:) = PrecisionArray;
    
    TIME_ePos(num) = toc(t_epos);   % 新增：当前 ePos/test drug 的运行时间
    
    % weight=1时等价于原写法；保留weight更通用
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

    %% ===== 新增：Denovo ePos-level mean ± SD 统计 =====
    % RowAucValue / RowAuPRValue / RowPrecisionValue：
    % 每一个元素对应一个 ePos/test drug 的独立评价结果。
    % AUPR 输出统计使用每个目标药物单独按照 R_m_A_AUPR_value 逻辑计算得到的 Row_R_m_A_AUPR_value。
    % 注意：这里不改变最终整体 R_m_A_AUPR_value 的原始计算逻辑。
    
    Denovo_AUC_mean = mean(RowAucValue);
    Denovo_AUC_SD   = std(RowAucValue);
    
    % 正式输出的 Denovo AUPR mean ± SD
    Denovo_AUPR_mean = R_m_A_AUPR_value;
	Denovo_AUPR_SD   = std(Row_R_m_A_AUPR_value);
    
    Denovo_Precision_mean = mean(RowPrecisionValue);
    Denovo_Precision_SD   = std(RowPrecisionValue);
    
    % 同时保留未修正AUPR，方便对比
    Denovo_raw_AUPR_mean = mean(RowAuPRValue);
    Denovo_raw_AUPR_SD   = std(RowAuPRValue);
    
    DenovoStats.AUC_mean = Denovo_AUC_mean;
    DenovoStats.AUC_SD = Denovo_AUC_SD;
    
    DenovoStats.AUPR_mean = Denovo_AUPR_mean;
    DenovoStats.AUPR_SD = Denovo_AUPR_SD;
    DenovoStats.Row_R_m_A_AUPR_value = Row_R_m_A_AUPR_value;
    
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
        
    parameters.HGBI_alpha = HGBI_alpha;
    parameters.BNNR_alpha = alpha;
    parameters.BNNR_beta = beta;
    parameters.threshold = threshold;
    parameters.sigma = sigma;
    parameters.weight = weight;
    parameters.maxiter = maxiter;
    parameters.tol1 = tol1;
    parameters.tol2 = tol2;
    
    
%% 4 测试结果保存和记录     

%根据不同数据集替换 Fdataset  Cdataset CTDdataset2023
save CTDdataset2023_STresult_HGIMC_Denovoone ...
    p_drugLen p_drugPos parameters ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value ...
    R_m_TPRArray R_m_PreArray R_m_A_AUPR_value ...
    RowAucValue RowAuPRValue RowPrecisionValue ...
    n_RowAucValue n_RowAuPRValue n_RowPrecisionValue ...
    Row_R_m_A_AUPR_value ...
    Denovo_AUC_mean Denovo_AUC_SD ...
    Denovo_AUPR_mean Denovo_AUPR_SD ...
    Denovo_raw_AUPR_mean Denovo_raw_AUPR_SD ...
    Denovo_Precision_mean Denovo_Precision_SD ...
    DenovoStats TIME_ePos runningtime ...
    alpha beta HGBI_alpha threshold sigma weight tol1 tol2 maxiter

 
