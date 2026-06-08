%%==================== MSBMF ======总共只有三个参数==============%%
clear all
rng('default');
%% 1.载入数据                                            1.载入数据
load Fdataset
% load Cdataset
% load CTDdataset2023

Wrr1=drug_ChemS;
Wrr2=drug_AtcS;
Wrr3=drug_SideS;
Wrr4=drug_DDIS;
Wrr5=drug_TargetS;
Wrr=[Wrr1,Wrr2,Wrr3,Wrr4,Wrr5];
Wdd1=disease_PhS;
Wdd2=disease_DoS;
Wdd=[Wdd1,Wdd2];
Wdr=didr;
Wrd = Wdr';

%%
% addpath('VDdataset');
% load drugsim;
% load virussim;
% load virusdrug;
% 
%  Wrr=virussim;
%  Wdd=drugsim;
%  Wdr=virusdrug;
%  Wrd = Wdr';

%% 2.参数赋值                                           2. 参数赋值
tol1=2*1e-3;%2*1e-2; %1e-5;%            %1e-2
tol2=1e-4;
maxiter=300;          %300
%% 3.进行多次十倍交叉验证                           3.进行十倍交叉验证
%%%%%%%%%%%%%% set the CV parameters %%%%%%%%%%
Count_CV =10;     %几次
nCV = 10;         %十倍交叉

PosMat = find(Wdr==1);%全部1
NumAs = length(PosMat);
[dn,dr] = size(Wdr);

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
iter_alone = zeros(Count_CV, nCV);

TIME_alone = zeros(Count_CV, nCV);        % 每一折运行时间
runningtime_10CV = zeros(Count_CV, 1);    % 每一次10CV运行时间
t_all = tic;                              % 总运行时间

koko=1;
min_mn=min(dn,dr);
for lambda3=0.1%[ 0.001 0.01 0.1 1]
    for lambda1=0.01%[0.001 0.01 0.1 1]
        for numk=[0.7]
            lambda2=lambda1;
            lambda33=lambda3;
            lambda4=lambda1;
            lambda44=lambda4;
            rankk=floor(min_mn*numk);
            
            for num = 1:Count_CV
                num
                t_cv = tic;   % 当前第num次10CV运行时间
                rand('state', num); %#ok<RAND>
                
                random_indices = randperm(NumAs);
                random_indices(NumAs+1:T_NumAs) = 0;
                
                Indices_groups = reshape(random_indices(1:floor(length(random_indices)/nCV)*nCV), nCV, floor(length(random_indices)/nCV));
                
                C_DresultMat_TPR = zeros(NumAs,dn);%对每个样本点求TPRR
                C_DresultMat_FPR = zeros(NumAs,dn);
                C_DresultMat_Pre = zeros(NumAs,dn);
                
                ass_num = 1;
                
                for i = 1:nCV
                    t_fold = tic;   % 当前fold运行时间
                    G_TestIds = Indices_groups(i,:);
                    G_TestIds(G_TestIds==0) = [];%%排除最后几个可能的空值
                    %%%%%%%%%% Tfnum indicates the number of elements in each group %%%%%%%%%%
                    Tfnum = length(G_TestIds);
                    TestIds = PosMat(G_TestIds);%%%%%%测试集向量
                    
                    P_TMat = Wdr;
                    P_TMat(TestIds) = 0;%%训练集
                    
                    T=P_TMat;
                    [U,V,iter,stop1,stop2]=fMSNMF2_R_mu(T,Wdd,Wrr,lambda1,lambda2,lambda3,lambda33,lambda4,lambda44,rankk,tol1,tol2,maxiter);
      
                    
                    M_ResultMat=U*V';
                    iter
                    [AUC_alone(num,i),Precision_alone(num,i),AUPR_alone(num,i)]=Fun_Auc3(M_ResultMat,P_TMat,TestIds)
                    iter_alone(num,i)=iter;
                    TIME_alone(num,i) = toc(t_fold);
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
                end
                A_Result_TPRArray(num,:) = mean(C_DresultMat_TPR);%10*313 一次10CV中 平均ten folds
                A_Result_FPRArray(num,:) = mean(C_DresultMat_FPR);
                A_Result_PreArray(num,:) = mean(C_DresultMat_Pre);
                
                A_AUC_values(num,:) = trapz(A_Result_FPRArray(num,:),A_Result_TPRArray(num,:))
                A_AUPR_values(num,:) = trapz(A_Result_TPRArray(num,:),A_Result_PreArray(num,:))
                m_TPRArray = [0,A_Result_TPRArray(num,:)];
                m_PreArray = [1,A_Result_PreArray(num,:)];
                m_A_AUPR_values(num,:) = trapz(m_TPRArray,m_PreArray);
                runningtime_10CV(num) = toc(t_cv);
            end
            
            Result_TPRArray =mean(A_Result_TPRArray);%VV%    十次10CV中 平均ten folds    对于只跑一次，去掉mean
            Result_FPRArray =mean(A_Result_FPRArray);%VV
            Result_PreArray =mean(A_Result_PreArray);%VV
            %             T_TPRArray = [0,Result_TPRArray];
            %             T_PreArray = [1,Result_PreArray];
            %             T_AUPR_value=trapz(T_TPRArray,T_PreArray);
            Result_AUC_value = trapz(Result_FPRArray,Result_TPRArray);%VV
            Result_AUPR_value = trapz(Result_TPRArray,Result_PreArray);%V
            Result_Precision_value=Result_PreArray(1);%VV
            
            R_m_TPRArray = [0,Result_TPRArray];%V
            R_m_PreArray = [1,Result_PreArray];%V
            R_m_A_AUPR_value = trapz(R_m_TPRArray,R_m_PreArray);%V
            
            %% ===== 新增：10次10倍交叉验证的 mean ± SD 统计 =====
            AUC_10CV_mean = mean(A_AUC_values);
            AUC_10CV_SD   = std(A_AUC_values);
            
            AUPR_10CV_mean = mean(m_A_AUPR_values);
            AUPR_10CV_SD   = std(m_A_AUPR_values);
            
            % 每次10CV的Precision使用该次平均Precision曲线的第一个点
            Precision_values = A_Result_PreArray(:,1);
            Precision_10CV_mean = mean(Precision_values);
            Precision_10CV_SD   = std(Precision_values);
            
            % fold-level统计：10 × 10 = 100个fold
            AUC_fold_mean = mean(AUC_alone(:));
            AUC_fold_SD   = std(AUC_alone(:));
            
            AUPR_fold_mean = mean(AUPR_alone(:));
            AUPR_fold_SD   = std(AUPR_alone(:));
            
            Precision_fold_mean = mean(Precision_alone(:));
            Precision_fold_SD   = std(Precision_alone(:));
            
            Iter_fold_mean = mean(iter_alone(:));
            Iter_fold_SD   = std(iter_alone(:));
            
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
            
            MetricStats.Iter_fold_mean = Iter_fold_mean;
            MetricStats.Iter_fold_SD = Iter_fold_SD;
            
            % 小数点后4位，方便论文表格复制
            MetricStats.AUC_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUC_10CV_mean, AUC_10CV_SD);
            MetricStats.AUPR_10CV_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_10CV_mean, AUPR_10CV_SD);
            MetricStats.Precision_10CV_mean_SD_text = sprintf('%.4f ± %.4f', Precision_10CV_mean, Precision_10CV_SD);
            
            MetricStats.AUC_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUC_fold_mean, AUC_fold_SD);
            MetricStats.AUPR_fold_mean_SD_text = sprintf('%.4f ± %.4f', AUPR_fold_mean, AUPR_fold_SD);
            MetricStats.Precision_fold_mean_SD_text = sprintf('%.4f ± %.4f', Precision_fold_mean, Precision_fold_SD);
            
            MetricStats.Iter_fold_mean_SD_text = sprintf('%.4f ± %.4f', Iter_fold_mean, Iter_fold_SD);
            
            time_all = toc(t_all);

        end
    end
end
parameters.lambda1=lambda1;
parameters.lambda3=lambda3;
parameters.numk=numk;
parameters.rankk=rankk;
parameters.maxiter=maxiter;
parameters.iter=iter;
parameters.tol1=tol1;
parameters.tol2=tol2;
fprintf('Curve-mean result: AUC=%4.4f, AUPR=%4.4f, Precision=%4.4f.\n', ...
    Result_AUC_value, R_m_A_AUPR_value, Result_Precision_value);

fprintf('10CV-level mean ± SD: AUC=%s, AUPR=%s, Precision=%s.\n', ...
    MetricStats.AUC_10CV_mean_SD_text, ...
    MetricStats.AUPR_10CV_mean_SD_text, ...
    MetricStats.Precision_10CV_mean_SD_text);

%% 4 测试结果保存和记录

%根据不同数据集替换 Fdataset  Cdataset CTDdataset2023
save Fdataset_STresult_MSBMF_10CV_fold_results.mat ...
    time_all runningtime_10CV TIME_alone MetricStats ...
    parameters AUC_alone AUPR_alone Precision_alone iter_alone ...
    A_AUC_values A_AUPR_values m_A_AUPR_values Precision_values ...
    AUC_10CV_mean AUC_10CV_SD ...
    AUPR_10CV_mean AUPR_10CV_SD ...
    Precision_10CV_mean Precision_10CV_SD ...
    AUC_fold_mean AUC_fold_SD ...
    AUPR_fold_mean AUPR_fold_SD ...
    Precision_fold_mean Precision_fold_SD ...
    Iter_fold_mean Iter_fold_SD ...
    Result_TPRArray Result_FPRArray Result_PreArray ...
    Result_AUC_value Result_AUPR_value Result_Precision_value ...
    R_m_TPRArray R_m_PreArray R_m_A_AUPR_value ...
    lambda1 lambda2 lambda3 lambda33 lambda4 lambda44 numk rankk tol1 tol2 maxiter