%%==================== multiGMF_10CV ====================%%
% This script performs the standard 10-fold cross-validation evaluation.

clear
addpath('Datasets');
rng('default');
%% 1.载入数据                                            
load Fdataset
% load Cdataset 
% load CTDdataset2023

%%

drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

Wdr=didr; % disease-by-drug association matrix in the dataset
Wrd = Wdr'; % drug-by-disease association matrix used by fmultiGMF
tic
[dn,dr] = size(Wdr);

 Wrr=(drug_ChemS+drug_AtcS+drug_SideS+drug_DDIS+drug_TargetS)/5;
 Wdd=(disease_PhS+disease_DoS)/2;

%% 2.参数赋值                                           2. 参数赋值 
out_t = [];%记录每次迭代终止的次数
Preci_ruselt = [];%记录一折十次的各参数值

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
    t_cv = tic;   % 新增：记录第num次10CV的运行时间

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

% Full multiGMF setting used in the main experiments.
% Drug similarities: R1 chemical, R2 ATC, R3 side effect, R4 DDI, R5 target.
% Disease similarities: D1 phenotype, D2 Disease Ontology.
Input.A={drug_ChemS,drug_AtcS,drug_SideS,drug_DDIS,drug_TargetS};
Input.B={disease_PhS,disease_DoS};


%低秩矩阵初始化
min_mn = min(dn,dr);

Input.X = Wdr;%取关联矩阵中已知关联位置需要用到Wdr
Input.Wrr = Wrr;
Input.Wdd = Wdd;
Input.P_TMat = P_TMat';

% 单相似性对比
% Single-similarity ablation experiments reported in mutliGMF.
% To reproduce a single-source variant, uncomment one selected drug similarity and one selected disease similarity below, and set the corresponding Wrr and Wdd.
% This block should remain commented when running the full multiGMF model.

% Input.A = {drug_TargetS};%{drug_ChemS,drug_AtcS,drug_SideS,drug_DDIS,drug_TargetS}
% Input.B = {disease_PhS};%{disease_PhS,disease_DoS}
% 
% Input.Wrr = drug_ChemS;%{drug_ChemS,drug_AtcS,drug_SideS,drug_DDIS,drug_TargetS}
% Input.Wdd = disease_PhS;%{disease_PhS,disease_DoS}




%% Fdataset  knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.0001 lambda3= 1 Options.MaxIter=300 Options.tol1 = 2*1e-3; Options.tol2 = 1*1e-4;
lambda_soft = 1;
lambda1 = 0.0001;%0.0001
lambda3 = 1;
KNN_K = 10;%10
tau = 0.7;
Options.MaxIter=300; 
Options.tol1 = 2*1e-3;
Options.tol2 = 1*1e-4;

%% Cdataset   knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.0001 lambda3= 1 Options.MaxIter=300 Options.tol1 = 2*1e-3; Options.tol2 = 1*1e-4;
% lambda_soft = 1;
% lambda1 = 0.0001;
% lambda3 = 1;
% KNN_K = 10;
% tau = 0.7;
% Options.MaxIter=300; 
% Options.tol1 = 2*1e-3;
% Options.tol2 = 1*1e-4;

%% CTDdataset2023   knn = 10 tau = 0.7 lambda_soft= 1 lambda1= 0.0001 lambda3= 1 Options.MaxIter=300 Options.tol1 = 2*1e-3; Options.tol2 = 1*1e-4;
% lambda_soft = 1;
% lambda1 = 0.0001;
% lambda3 = 1;
% KNN_K = 10;
% tau = 0.7;
% Options.MaxIter=300;
% Options.tol1 = 2*1e-3;
% Options.tol2 = 1*1e-4;


Options.lambda_soft=lambda_soft;
Options.lambda1=lambda1;
Options.lambda3=lambda3;

Options.lambad4 = 1;
Options.mu1 = 1;
Options.mu2 = 1;

Input.kk = KNN_K;

rankk = floor(min_mn*tau);

Input.HInit = rand(dn,rankk);
Input.WInit = rand(dr,rankk);

%执行迭代函数
Output = fmultiGMF(Input, Options);% Run the full multiGMF model. This is the default setting used for the main 10CV results.

% No-soft ablation reported in multiGMF. Uncomment the following line and comment the full model line above
% only when reproducing the multiGMF without soft coupling experiment.
% Output = no_soft_fmultiGMF(Input, Options);

%矩阵填充复原
WW=Output.W*Output.H';
M_ResultMat=WW';      

Output.t %显示最终迭代次数
%% 记录每次迭代终止的次数
out_t = [Output.t,out_t];


[AUC_alone(num,i),Precision_alone(num,i),AUPR_alone(num,i)]=Fun_Auc3(M_ResultMat,P_TMat,TestIds)
TIME_alone(num,i) = toc(t_fold);

%% 存储每一次各参数的值以及AUC等值
AAres(1,ppp)=AUC_alone(num,i);
AAres(2,ppp)=Precision_alone(num,i);
AAres(3,ppp)=AUPR_alone(num,i);
AAres(4:5,ppp)=[tau,KNN_K];
AAres(6:8,ppp)=[lambda_soft,lambda1,lambda3];

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
            time_1=toc
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

    Precision_values(num,:) = mean(Precision_alone(num,:));
    %% 查看一折十次之后各参数和AUC等值
    parameter(pp,1:3) = [lambda_soft,lambda1,lambda3];
    parameter(pp,4:9) = [A_AUC_values(num,:),Precision_values(num,:),m_A_AUPR_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:),A_AUC_values(num,:) + Precision_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:) + Precision_values(num,:)];
    num_knn(pp,1:2) = [tau,KNN_K];
    num_knn(pp,3:8) = [A_AUC_values(num,:),Precision_values(num,:),m_A_AUPR_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:),A_AUC_values(num,:) + Precision_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:) + Precision_values(num,:)];
    result_lambda1(pp,1:3) = [lambda1,tau,KNN_K];
    result_lambda1(pp,4:9) = [A_AUC_values(num,:),Precision_values(num,:),m_A_AUPR_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:),A_AUC_values(num,:) + Precision_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:) + Precision_values(num,:)];
    result_lambda3(pp,1:3) = [tau,lambda1,lambda3];
    result_lambda3(pp,4:9) = [A_AUC_values(num,:),Precision_values(num,:),m_A_AUPR_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:),A_AUC_values(num,:) + Precision_values(num,:),A_AUC_values(num,:) + m_A_AUPR_values(num,:) + Precision_values(num,:)];
    % 
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

%% ===== 新增：10次10倍交叉验证的 mean 和 SD 统计 =====
% A_AUC_values / m_A_AUPR_values / Precision_values 是每一次10CV汇总后的结果，共 Count_CV=10 个值。
% 因此这里计算的是“10次10CV结果”的 mean ± SD。

AUC_10CV_mean = mean(A_AUC_values);
AUC_10CV_SD   = std(A_AUC_values);

AUPR_10CV_mean = mean(m_A_AUPR_values);
AUPR_10CV_SD   = std(m_A_AUPR_values);

Precision_10CV_mean = mean(Precision_values);
Precision_10CV_SD   = std(Precision_values);

% 同时保留 fold-level 的总体统计，即 10 × 10 = 100 个fold结果
AUC_fold_mean = mean(AUC_alone(:));
AUC_fold_SD   = std(AUC_alone(:));

AUPR_fold_mean = mean(AUPR_alone(:));
AUPR_fold_SD   = std(AUPR_alone(:));

Precision_fold_mean = mean(Precision_alone(:));
Precision_fold_SD   = std(Precision_alone(:));

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

MetricStats.AUC_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUC_10CV_mean, AUC_10CV_SD);
MetricStats.AUPR_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_10CV_mean, AUPR_10CV_SD);
MetricStats.Precision_10CV_mean_SD_text = sprintf('%.4f ± %.4f', Precision_10CV_mean, Precision_10CV_SD);

MetricStats.AUC_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUC_fold_mean, AUC_fold_SD);
MetricStats.AUPR_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_fold_mean, AUPR_fold_SD);
MetricStats.Precision_fold_mean_SD_text = sprintf('%.4f ± %.4f', Precision_fold_mean, Precision_fold_SD);

time_all = toc(t_all);

parameters.lambda_soft=lambda_soft;
parameters.lambda1=lambda1;
fprintf('AUC=%4.4f,AUPR=%4.4f,Precision=%4.4f.\n',Result_AUC_value,R_m_A_AUPR_value,Result_Precision_value)
%% 4 测试结果保存和记录    
tol1 = Options.tol1;
tol2 = Options.tol2;

% 不同数据集更换mat文件的名称 Fdataset Cdataset CTDdataset2023
save Fdataset_multiGMF_10CV_fold_results.mat ...
    time_1 time_all runningtime_10CV TIME_alone parameters MetricStats ...
    AUC_alone AUPR_alone Precision_alone AAres ...
    A_AUC_values A_AUPR_values m_A_AUPR_values Precision_values ...
    AUC_10CV_mean AUC_10CV_SD ...
    AUPR_10CV_mean AUPR_10CV_SD ...
    Precision_10CV_mean Precision_10CV_SD ...
    AUC_fold_mean AUC_fold_SD ...
    AUPR_fold_mean AUPR_fold_SD ...
    Precision_fold_mean Precision_fold_SD ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value R_m_A_AUPR_value ...
    out_t lambda_soft lambda1 lambda3 tau KNN_K tol1 tol2
