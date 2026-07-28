%%==================== DDASKF ====================%%
% 注意：1.rng('default')和rand('state', num); %#ok<RAND>
%       2.保存名称
clear
addpath('code');
rng('default');

%% 1.载入数据                                            

% load Fdataset;
load Cdataset;
% load CTDdataset2023

drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

Input.drugChemS = drug_ChemS;
Input.drugDDIS = drug_DDIS;
Input.drugAtcS = drug_AtcS;
Input.drugSideS = drug_SideS;
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

Preci_ruselt = [];
%% 3.进行多次十倍交叉验证                           3.进行十倍交叉验证
%%%%%%%%%%%%%% set the CV parameters %%%%%%%%%%
Count_CV =10;     %几次
nCV = 10;         %十倍交叉

PosMat = find(Wdr==1);
NumAs = length(PosMat);


T_NumAs = ceil(NumAs/nCV)*nCV;

A_Result_TPRArray = zeros(Count_CV,dn);%一个Wrr对所有Wdd（dn）进行排名
A_Result_FPRArray = zeros(Count_CV,dn);
A_Result_PreArray = zeros(Count_CV,dn);
A_AUC_values = zeros(Count_CV,1);
A_AUPR_values = zeros(Count_CV,1);
m_A_AUPR_values = zeros(Count_CV,1);

%% ===== 新增：fold-level指标与运行时间记录 =====
AUC_alone = zeros(Count_CV, nCV);
AUPR_alone = zeros(Count_CV, nCV);
Precision_alone = zeros(Count_CV, nCV);

TIME_alone = zeros(Count_CV, nCV);        % 每一折运行时间
runningtime_10CV = zeros(Count_CV, 1);    % 每一次10倍交叉验证运行时间
t_all = tic;                              % 10次10CV总运行时间

for num = 1:Count_CV
    num
    rand('state', num); %#ok<RAND>
    t_cv = tic;   % 新增：记录第num次10CV运行时间

    random_indices = randperm(NumAs);
    random_indices(NumAs+1:T_NumAs) = 0;
    
    Indices_groups = reshape(random_indices(1:floor(length(random_indices)/nCV)*nCV), nCV, floor(length(random_indices)/nCV));

    C_DresultMat_TPR = zeros(NumAs,dn);%对每个样本点求TPRR
    C_DresultMat_FPR = zeros(NumAs,dn);
    C_DresultMat_Pre = zeros(NumAs,dn);
    
    ass_num = 1;
    
    for i = 1:nCV
        t_fold = tic;   % 新增：记录当前fold运行时间
        
        G_TestIds = Indices_groups(i,:);
        G_TestIds(G_TestIds==0) = [];%%排除最后几个可能的空值
%%%%%%%%%% Tfnum indicates the number of elements in each group %%%%%%%%%%
        Tfnum = length(G_TestIds);
        TestIds = PosMat(G_TestIds);%%%%%%测试集
 
        P_TMat = Wdr;
        P_TMat(TestIds) = 0;%%训练集
        
%%
pp=1;
ppp=1;

%执行迭代函数
score_matrix = DDASKF(Input, P_TMat);

M_ResultMat=score_matrix;      

%% 记录每次迭代终止的次数
[AUC_alone(num,i),Precision_alone(num,i),AUPR_alone(num,i)]=Fun_Auc3(M_ResultMat,P_TMat,TestIds)
TIME_alone(num,i) = toc(t_fold);

%% 存储每一次各参数的值以及AUC等值
AAres(1,ppp)=AUC_alone(num,i);
AAres(2,ppp)=Precision_alone(num,i);
AAres(3,ppp)=AUPR_alone(num,i);

ppp=ppp+1;


%% 方便查看每折之后的Precision值
Preci_calues = Precision_alone(num,i);
Preci_ruselt = [Preci_ruselt,Preci_calues];

%% ----------------------------------------------------------
        
        thresh_value = min(M_ResultMat(:))-10;
        thresh_value = ceil(thresh_value);
        
        A_DresultMat_TPR = zeros(Tfnum,dn);
        A_DresultMat_FPR = zeros(Tfnum,dn);
        A_DresultMat_Pre = zeros(Tfnum,dn);
%%%%%%%%%%%%%%%% DresultMat stores the predicted lables(0/1)%%%%%%%%%%%%%%%
        DresultMat = zeros(dn,dr);
        
        Qvalue = M_ResultMat(TestIds);
        TPosMat = PosMat;
        TPosMat(G_TestIds) = [];
        M_ResultMat(TPosMat)= thresh_value;
        DresultMat(TPosMat)= thresh_value;
        S_ResultMat = sort(M_ResultMat,'descend');
        DresultMat = sort(DresultMat,'descend');
        
        for k=1:Tfnum
            rdPos = TestIds(k);
            rindex = ceil(rdPos/dn);
            eQvalue = Qvalue(k);
            TfindposMat =  find(S_ResultMat(:,rindex)==eQvalue);
            TfindposMatlen = size(TfindposMat,1);
            Tfindpos = TfindposMat(TfindposMatlen);
            
            result_Mat = DresultMat(:,rindex);
            result_Mat(result_Mat==thresh_value)=[];
            result_Mat(Tfindpos) = 1;
            result_len = length(result_Mat);
            
 %%%%% Construct TPRArray and FPRArray,and then compute the AUC value %%%%%
            TPRArray = zeros(1,dn);
            FPRArray = zeros(1,dn);
            PreArray = zeros(1,dn);
            
            CountP =  1;
            CountN =  result_len-1;
            
            Tpnum = 0;
            Fpnum = 0;
            for m =1:result_len
                if(result_Mat(m)==1)
                    Tpnum = Tpnum + 1;
                else
                    Fpnum = Fpnum + 1;
                end
                TPRArray(m) = Tpnum/CountP;
                FPRArray(m) = Fpnum/CountN;
                PreArray(m) = Tpnum/(Tpnum+Fpnum);
            end
            TPRArray(result_len+1:dn) = TPRArray(result_len);
            FPRArray(result_len+1:dn) = FPRArray(result_len);
            PreArray(result_len+1:dn) = PreArray(result_len);
            
            A_DresultMat_TPR(k,:) = TPRArray;%193*313 one fold中每个关联排序一次
            A_DresultMat_FPR(k,:) = FPRArray;
            A_DresultMat_Pre(k,:) = PreArray;
        end 
  
        C_DresultMat_TPR(ass_num:ass_num+Tfnum-1,:) = A_DresultMat_TPR;%1993*313 一次10CV中 累计ten folds
        C_DresultMat_FPR(ass_num:ass_num+Tfnum-1,:) = A_DresultMat_FPR;
        C_DresultMat_Pre(ass_num:ass_num+Tfnum-1,:) = A_DresultMat_Pre;
        
        ass_num = ass_num+Tfnum;
        if i==1
            time_1 = TIME_alone(num,i);
        end
    end
    A_Result_TPRArray(num,:) = mean(C_DresultMat_TPR);%一次10CV中 平均ten folds
    A_Result_FPRArray(num,:) = mean(C_DresultMat_FPR);
    A_Result_PreArray(num,:) = mean(C_DresultMat_Pre);
   
    A_AUC_values(num,:) = trapz(A_Result_FPRArray(num,:),A_Result_TPRArray(num,:))
    A_AUPR_values(num,:) = trapz(A_Result_TPRArray(num,:),A_Result_PreArray(num,:));
    m_TPRArray = [0,A_Result_TPRArray(num,:)];
    m_PreArray = [1,A_Result_PreArray(num,:)];
    m_A_AUPR_values(num,:) = trapz(m_TPRArray,m_PreArray)
    runningtime_10CV(num) = toc(t_cv);

    %% 杨斌加
    Precision_values(num,:) = mean(Precision_alone(num,:));
    %% 查看一折十次之后各参数和AUC等值
%     parameter(pp,1:5) = [A_AUC_values(num,:),Precision_values(num,:),m_A_AUPR_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:),A_AUC_values(num,:) + Precision_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:) + Precision_values(num,:)];
end

Result_TPRArray =mean(A_Result_TPRArray);%VV%    对于只跑一次，去掉mean
Result_FPRArray =mean(A_Result_FPRArray);%VV
Result_PreArray =mean(A_Result_PreArray);%VV

Result_AUC_value = trapz(Result_FPRArray,Result_TPRArray);%VV
Result_AUPR_value = trapz(Result_TPRArray,Result_PreArray);%V
Result_Precision_value=Result_PreArray(1);%VV

R_m_TPRArray = [0,Result_TPRArray];%V
R_m_PreArray = [1,Result_PreArray];%V
R_m_A_AUPR_value = trapz(R_m_TPRArray,R_m_PreArray);%V

%% ===== 新增：10次10倍交叉验证的 mean ± SD 统计 =====
% A_AUC_values / m_A_AUPR_values / Precision_values：
% 每一行对应一次完整10-fold CV，因此共有 Count_CV=10 个结果。

AUC_10CV_mean = mean(A_AUC_values);
AUC_10CV_SD   = std(A_AUC_values);

AUPR_10CV_mean = mean(m_A_AUPR_values);
AUPR_10CV_SD   = std(m_A_AUPR_values);

Precision_10CV_mean = mean(Precision_values);
Precision_10CV_SD   = std(Precision_values);

% fold-level总体统计：10次 × 10折 = 100个fold结果
AUC_fold_mean = mean(AUC_alone(:));
AUC_fold_SD   = std(AUC_alone(:));

AUPR_fold_mean = mean(AUPR_alone(:));
AUPR_fold_SD   = std(AUPR_alone(:));

Precision_fold_mean = mean(Precision_alone(:));
Precision_fold_SD   = std(Precision_alone(:));

% 结构体保存，便于后续 load mat 文件后读取
MetricStats.AUC_10CV_mean = AUC_10CV_mean;
MetricStats.AUC_10CV_SD = AUC_10CV_SD;
MetricStats.AUPR_10CV_mean = AUPR_10CV_mean;
MetricStats.AUPR_10CV_SD = AUPR_10CV_SD;
MetricStats.Precision_10CV_mean = Precision_10CV_mean;
MetricStats.Precision_10CV_SD = Precision_10CV_SD;

MetricStats.AUC_fold_mean = AUC_fold_mean;
MetricStats.AUC_fold_SD = AUC_fold_SD;
MetricStats.AUPR_fold_mean = AUPR_fold_mean;
MetricStats.AUPR_fold_SD = AUPR_fold_SD;
MetricStats.Precision_fold_mean = Precision_fold_mean;
MetricStats.Precision_fold_SD = Precision_fold_SD;

% 字符串形式，默认小数点后4位，方便论文表格直接复制
MetricStats.AUC_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUC_10CV_mean, AUC_10CV_SD);
MetricStats.AUPR_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_10CV_mean, AUPR_10CV_SD);
MetricStats.Precision_10CV_mean_SD_text = sprintf('%.4f ± %.4f', Precision_10CV_mean, Precision_10CV_SD);

MetricStats.AUC_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUC_fold_mean, AUC_fold_SD);
MetricStats.AUPR_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_fold_mean, AUPR_fold_SD);
MetricStats.Precision_fold_mean_SD_text = sprintf('%.4f ± %.4f', Precision_fold_mean, Precision_fold_SD);

time_all = toc(t_all);

fprintf('Curve-mean result: AUC=%4.4f, AUPR=%4.4f, Precision=%4.4f.\n', ...
    Result_AUC_value, R_m_A_AUPR_value, Result_Precision_value);

fprintf('10CV-level mean ± SD: AUC=%s, AUPR=%s, Precision=%s.\n', ...
    MetricStats.AUC_10CV_mean_SD_text, ...
    MetricStats.AUPR_10CV_mean_SD_text, ...
    MetricStats.Precision_10CV_mean_SD_text);

%% 4 测试结果保存和记录    

%根据不同数据集替换 Fdataset  Cdataset CTDdataset2023
save Cdataset_STresult_DDASKF_10CV_fold_results.mat ...
    time_1 time_all runningtime_10CV TIME_alone MetricStats ...
    AUC_alone AUPR_alone Precision_alone AAres ...
    A_AUC_values A_AUPR_values m_A_AUPR_values Precision_values ...
    AUC_10CV_mean AUC_10CV_SD ...
    AUPR_10CV_mean AUPR_10CV_SD ...
    Precision_10CV_mean Precision_10CV_SD ...
    AUC_fold_mean AUC_fold_SD ...
    AUPR_fold_mean AUPR_fold_SD ...
    Precision_fold_mean Precision_fold_SD ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value R_m_A_AUPR_value