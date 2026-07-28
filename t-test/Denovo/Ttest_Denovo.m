%% best_method_paired_ttest_holm_for_denovo_tables.m
% =========================================================================
% Purpose:
%   Denovo/cold-start setting:
%   Identify the best-performing method for each metric and test whether
%   the best method is significantly better than baselines using paired
%   two-sided t-tests + Holm/Bonferroni correction.
%
% Important for Denovo:
%   AUC main value       -> Result_AUC_value
%   AUPR main value      -> R_m_A_AUPR_value
%   Precision main value -> Result_Precision_value
%
%   The paired t-test still needs per-test-drug/ePos vectors:
%   AUC vector       -> RowAucValue / rowAucValue / n_RowAucValue
%   AUPR vector      -> n_RowAuPRValue / n_rowAuPRValue / RowAuPRValue
%   Precision vector -> RowPrecisionValue / rowPrecisionValue / n_RowPrecisionValue
%
% Output:
%   1. Command-window summary: whether to add $^{*}$ to the best method.
%   2. CSV file containing raw p-values, adjusted p-values, and significance.
%   3. A manuscript-friendly summary for LaTeX tables.
% =========================================================================

clear; clc;

%% ===================== User settings =====================

script_dir = fileparts(mfilename('fullpath'));
data_dir = script_dir;

test_name = 'Denovo';

correction_method = 'holm';

alpha = 0.05;

compare_only_strong_baselines = false;

strong_baselines = {'MSBMF', 'HGIMC', 'BNNR'};

require_significant_against_all = true;

%% ===================== Method files =====================

files = struct();

%% ---------- Fdataset ----------
dataset_name = 'Fdataset';
files.multiGMF = fullfile(data_dir, 'Fdataset_STresult_multiGMF_Denovoone.mat');
files.DDASKF   = fullfile(data_dir, 'Fdataset_STresult_DDASKF_Denovoone.mat');
files.MSBMF    = fullfile(data_dir, 'Fdataset_STresult_MSBMF_Denovoone.mat');
files.HGIMC    = fullfile(data_dir, 'Fdataset_STresult_HGIMC_Denovoone.mat');
files.BNNR     = fullfile(data_dir, 'Fdataset_STresult_BNNR_Denovoone.mat');
files.MBiRW    = fullfile(data_dir, 'Fdataset_STresult_MBiRW_Denovoone.mat');

%% ---------- Cdataset ----------
% dataset_name = 'Cdataset';
% files.multiGMF = fullfile(data_dir, 'Cdataset_STresult_multiGMF_Denovoone.mat');
% files.DDASKF   = fullfile(data_dir, 'Cdataset_STresult_DDASKF_Denovoone.mat');
% files.MSBMF    = fullfile(data_dir, 'Cdataset_STresult_MSBMF_Denovoone.mat');
% files.HGIMC    = fullfile(data_dir, 'Cdataset_STresult_HGIMC_Denovoone.mat');
% files.BNNR     = fullfile(data_dir, 'Cdataset_STresult_BNNR_Denovoone.mat');
% files.MBiRW    = fullfile(data_dir, 'Cdataset_STresult_MBiRW_Denovoone.mat');

%% ---------- CTDdataset2023 ----------
% dataset_name = 'CTDdataset2023';
% files.multiGMF = fullfile(data_dir, 'CTDdataset2023_STresult_multiGMF_Denovoone.mat');
% files.DDASKF   = fullfile(data_dir, 'CTDdataset2023_STresult_DDASKF_Denovoone.mat');
% files.MSBMF    = fullfile(data_dir, 'CTDdataset2023_STresult_MSBMF_Denovoone.mat');
% files.HGIMC    = fullfile(data_dir, 'CTDdataset2023_STresult_HGIMC_Denovoone.mat');
% files.BNNR     = fullfile(data_dir, 'CTDdataset2023_STresult_BNNR_Denovoone.mat');
% files.MBiRW    = fullfile(data_dir, 'CTDdataset2023_STresult_MBiRW_Denovoone.mat');

metric_names = {'AUC', 'AUPR', 'Precision'};

%% ===================== Load Denovo data =====================

method_names = fieldnames(files);
Data = struct();

fprintf('\n============================================================\n');
fprintf('Dataset: %s | Test: %s | Correction: %s\n', ...
    dataset_name, test_name, correction_method);
fprintf('Data directory: %s\n', data_dir);
fprintf('AUPR main value source: R_m_A_AUPR_value\n');
fprintf('Paired t-test AUPR vector source: n_RowAuPRValue or RowAuPRValue\n');
fprintf('============================================================\n');

for i = 1:numel(method_names)

    mname = method_names{i};
    fpath = files.(mname);

    if ~isfile(fpath)
        error('File not found for method %s:\n%s\nPlease check the Denovo file name in the Method files section.', ...
            mname, fpath);
    end

    fprintf('\n[LOAD] %s\n%s\n', mname, fpath);

    S = load(fpath);

    for k = 1:numel(metric_names)

        metric = metric_names{k};

        [vec, report_mean, report_sd, source_msg] = extract_denovo_metric(S, metric, mname);
        vec = vec(:);

        if isempty(vec)
            error('Cannot extract %s vector from %s.', metric, fpath);
        end

        Data.(mname).(metric).vec = vec;
        Data.(mname).(metric).report_mean = report_mean;
        Data.(mname).(metric).report_sd = report_sd;
        Data.(mname).(metric).source = source_msg;

        fprintf('  %-10s: report_mean = %.6f, report_SD = %.6f, n = %d | %s\n', ...
            metric, report_mean, report_sd, numel(vec), source_msg);
    end
end

%% ===================== Best-method paired t-test =====================

all_results = table();
star_summary = table();

fprintf('\n\n==================== Best method significance suggestions ====================\n');

for k = 1:numel(metric_names)

    metric = metric_names{k};

    % ---------- use manuscript/report mean to identify best method ----------
    report_means = nan(numel(method_names), 1);
    report_sds   = nan(numel(method_names), 1);
    ns           = nan(numel(method_names), 1);

    for i = 1:numel(method_names)
        mname = method_names{i};
        x = Data.(mname).(metric).vec;
        report_means(i) = Data.(mname).(metric).report_mean;
        report_sds(i)   = Data.(mname).(metric).report_sd;
        ns(i)           = sum(~isnan(x));
    end

    [best_mean, best_idx] = max(report_means);
    best_method = method_names{best_idx};
    best_vec = Data.(best_method).(metric).vec;
    best_sd = report_sds(best_idx);

    % ---------- decide compared methods ----------
    if compare_only_strong_baselines
        compared_methods = strong_baselines;
        compared_methods = setdiff(compared_methods, best_method, 'stable');
    else
        compared_methods = setdiff(method_names, best_method, 'stable');
    end

    raw_p = nan(numel(compared_methods), 1);
    adj_p = nan(numel(compared_methods), 1);
    tstat = nan(numel(compared_methods), 1);
    df    = nan(numel(compared_methods), 1);

    compared_mean = nan(numel(compared_methods), 1);
    compared_sd   = nan(numel(compared_methods), 1);
    mean_diff     = nan(numel(compared_methods), 1);
    n_pair        = nan(numel(compared_methods), 1);
    sig_each      = false(numel(compared_methods), 1);

    % ---------- paired t-tests ----------
    for j = 1:numel(compared_methods)

        cname = compared_methods{j};
        comp_vec = Data.(cname).(metric).vec;

        [x2, y2] = align_pair_vectors(best_vec, comp_vec);

        valid = ~isnan(x2) & ~isnan(y2);
        x2 = x2(valid);
        y2 = y2(valid);

        n_pair(j) = numel(x2);
        compared_mean(j) = Data.(cname).(metric).report_mean;
        compared_sd(j)   = Data.(cname).(metric).report_sd;
        mean_diff(j)     = mean(x2 - y2, 'omitnan');

        if numel(x2) < 2
            warning('Not enough paired samples for %s vs %s on %s.', ...
                best_method, cname, metric);
            continue;
        end

        diff_vec = x2 - y2;

        if all(abs(diff_vec - diff_vec(1)) < 1e-15)
            warning('Constant paired differences for %s vs %s on %s. t-test may be undefined.', ...
                best_method, cname, metric);
        end

        [~, p, ~, stats] = ttest(x2, y2, ...
            'Alpha', alpha, ...
            'Tail', 'both');

        raw_p(j) = p;
        tstat(j) = stats.tstat;
        df(j)    = stats.df;
    end

    % ---------- multiple comparison correction ----------
    switch lower(correction_method)
        case 'holm'
            adj_p = holm_adjust(raw_p);
        case 'bonferroni'
            adj_p = min(raw_p * sum(~isnan(raw_p)), 1);
        otherwise
            error('Unknown correction method: %s', correction_method);
    end

    % 显著性的方向必须是最优方法均值更高
    sig_each = (adj_p < alpha) & (mean_diff > 0);

    valid_sig = ~isnan(adj_p);
    if isempty(sig_each) || ~any(valid_sig)
        add_star = false;
    else
        if require_significant_against_all
            add_star = all(sig_each(valid_sig));
        else
            add_star = any(sig_each(valid_sig));
        end
    end

    % ---------- print summary ----------
    fprintf('\nMetric: %s\n', metric);
    fprintf('  Best method: %s | report_mean = %.6f | report_SD = %.6f | n = %d\n', ...
        best_method, best_mean, best_sd, ns(best_idx));

    if add_star
        fprintf('  Suggested table mark: ADD $^{*}$ to %s for %s\n', ...
            best_method, metric);
    else
        fprintf('  Suggested table mark: no $^{*}$ for this row\n');
    end

    for j = 1:numel(compared_methods)
        fprintf('    vs %-8s: compared_mean = %.6f, mean_diff_vec = %.6f, p_raw = %.4g, p_adj = %.4g, sig = %d, n_pair = %d\n', ...
            compared_methods{j}, compared_mean(j), mean_diff(j), raw_p(j), adj_p(j), sig_each(j), n_pair(j));
    end

    % ---------- save detailed results ----------
    T = table();
    T.Dataset       = repmat({dataset_name}, numel(compared_methods), 1);
    T.Test          = repmat({test_name}, numel(compared_methods), 1);
    T.Metric        = repmat({metric}, numel(compared_methods), 1);
    T.BestMethod    = repmat({best_method}, numel(compared_methods), 1);
    T.BestMean      = repmat(best_mean, numel(compared_methods), 1);
    T.BestSD        = repmat(best_sd, numel(compared_methods), 1);
    T.ComparedWith  = compared_methods(:);
    T.ComparedMean  = compared_mean;
    T.ComparedSD    = compared_sd;
    T.MeanDiff      = mean_diff;
    T.N_pairs       = n_pair;
    T.t_stat        = tstat;
    T.df            = df;
    T.p_raw         = raw_p;
    T.p_adjusted    = adj_p;
    T.Significant   = sig_each;
    T.AddStarToBest = repmat(add_star, numel(compared_methods), 1);

    all_results = [all_results; T]; %#ok<AGROW>

    Srow = table();
    Srow.Dataset     = {dataset_name};
    Srow.Test        = {test_name};
    Srow.Metric      = {metric};
    Srow.BestMethod  = {best_method};
    Srow.BestMean    = best_mean;
    Srow.BestSD      = best_sd;
    Srow.AddStar     = add_star;
    Srow.TableString = {make_latex_value(best_mean, best_sd, add_star)};

    star_summary = [star_summary; Srow]; %#ok<AGROW>
end

%% ===================== Save results =====================

out_csv_detail = fullfile(data_dir, sprintf('%s_%s_best_method_paired_ttest_%s_detail.csv', ...
    dataset_name, test_name, correction_method));

out_csv_summary = fullfile(data_dir, sprintf('%s_%s_best_method_star_summary_%s.csv', ...
    dataset_name, test_name, correction_method));

writetable(all_results, out_csv_detail);
writetable(star_summary, out_csv_summary);

fprintf('\n\n==================== Manuscript table suggestion ====================\n');
disp(star_summary);

fprintf('\n[SAVED DETAIL]  %s\n', out_csv_detail);
fprintf('[SAVED SUMMARY] %s\n', out_csv_summary);

fprintf('\nDone.\n');

%% ========================================================================
%% Local functions
%% ========================================================================

function [vec, report_mean, report_sd, source_msg] = extract_denovo_metric(S, metric, method_name)
% Extract Denovo metric.
%
% Main manuscript values:
%   AUC       -> Result_AUC_value
%   AUPR      -> R_m_A_AUPR_value
%   Precision -> Result_Precision_value
%
% Per-ePos/test-drug vectors for paired t-test:
%   AUC       -> RowAucValue / rowAucValue / n_RowAucValue
%   AUPR      -> n_RowAuPRValue / n_rowAuPRValue / RowAuPRValue / rowAuPRValue
%   Precision -> RowPrecisionValue / rowPrecisionValue / n_RowPrecisionValue

    switch lower(metric)

        case 'auc'
            report_mean = get_required_scalar(S, {'Result_AUC_value', 'Denovo_AUC_mean'}, ...
                method_name, 'AUC report mean');
            vec = get_first_vector(S, {'RowAucValue', 'rowAucValue', 'n_RowAucValue', 'n_rowAucValue'}, ...
                method_name, 'AUC vector');
            report_sd = get_scalar_or_std(S, {'Denovo_AUC_SD', 'R_m_A_AUC_SD'}, vec);
            source_msg = 'main=Result_AUC_value; vector=RowAucValue';

        case 'aupr'
            % 用户指定：Denovo AUPR 主值使用 R_m_A_AUPR_value
            report_mean = get_required_scalar(S, {'R_m_A_AUPR_value'}, ...
                method_name, 'AUPR report mean');
            vec = get_first_vector(S, {'n_RowAuPRValue', 'n_rowAuPRValue', ...
                                       'RowAuPRValue', 'rowAuPRValue'}, ...
                method_name, 'AUPR vector');
            report_sd = get_scalar_or_std(S, {'R_m_A_AUPR_SD', 'Denovo_AUPR_SD', 'Denovo_raw_AUPR_SD'}, vec);
            source_msg = 'main=R_m_A_AUPR_value; vector=n_RowAuPRValue';

        case 'precision'
            report_mean = get_required_scalar(S, {'Result_Precision_value', 'Denovo_Precision_mean'}, ...
                method_name, 'Precision report mean');
            vec = get_first_vector(S, {'RowPrecisionValue', 'rowPrecisionValue', ...
                                       'n_RowPrecisionValue', 'n_rowPrecisionValue'}, ...
                method_name, 'Precision vector');
            report_sd = get_scalar_or_std(S, {'Denovo_Precision_SD', 'R_m_A_Precision_SD'}, vec);
            source_msg = 'main=Result_Precision_value; vector=RowPrecisionValue';

        otherwise
            error('Unknown metric: %s', metric);
    end

    vec = vec(:);
end

function value = get_required_scalar(S, candidate_names, method_name, description)

    value = [];

    for i = 1:numel(candidate_names)
        fname = candidate_names{i};
        if isfield(S, fname)
            tmp = S.(fname);
            if isnumeric(tmp) && isscalar(tmp)
                value = tmp;
                return;
            end
        end
    end

    error('[%s] Missing required scalar for %s. Candidate fields: %s', ...
        method_name, description, strjoin(candidate_names, ', '));
end

function vec = get_first_vector(S, candidate_names, method_name, description)

    vec = [];

    for i = 1:numel(candidate_names)
        fname = candidate_names{i};
        if isfield(S, fname)
            tmp = S.(fname);
            if isnumeric(tmp) && numel(tmp) >= 2
                vec = tmp(:);
                return;
            end
        end
    end

    error('[%s] Missing required vector for %s. Candidate fields: %s', ...
        method_name, description, strjoin(candidate_names, ', '));
end

function sd_value = get_scalar_or_std(S, candidate_names, vec)

    sd_value = NaN;

    for i = 1:numel(candidate_names)
        fname = candidate_names{i};
        if isfield(S, fname)
            tmp = S.(fname);
            if isnumeric(tmp) && isscalar(tmp)
                sd_value = tmp;
                return;
            end
        end
    end

    sd_value = std(vec, 'omitnan');
end

function [x2, y2] = align_pair_vectors(x, y)
% Align paired vectors.
% If lengths differ, truncate to the smaller length.

    x = x(:);
    y = y(:);

    nx = numel(x);
    ny = numel(y);

    if nx ~= ny
        n = min(nx, ny);
        warning('Paired vectors have different lengths: %d vs %d. Truncated to %d.', ...
            nx, ny, n);
        x2 = x(1:n);
        y2 = y(1:n);
    else
        x2 = x;
        y2 = y;
    end
end

function p_adj = holm_adjust(p)
% Holm-Bonferroni adjusted p-values.
%
% Input:
%   p      raw p-values
% Output:
%   p_adj  adjusted p-values in original order

    p = p(:);
    p_adj = nan(size(p));

    valid = ~isnan(p);
    p_valid = p(valid);

    m = numel(p_valid);

    if m == 0
        return;
    end

    [p_sorted, order] = sort(p_valid, 'ascend');

    adj_sorted = nan(size(p_sorted));

    for i = 1:m
        adj_sorted(i) = (m - i + 1) * p_sorted(i);
    end

    % enforce monotonicity
    for i = 2:m
        adj_sorted(i) = max(adj_sorted(i), adj_sorted(i - 1));
    end

    adj_sorted = min(adj_sorted, 1);

    temp = nan(m, 1);
    temp(order) = adj_sorted;

    idx_valid = find(valid);
    p_adj(idx_valid) = temp;
end

function s = make_latex_value(mean_value, sd_value, add_star)
% Generate LaTeX string for manuscript table.

    if isnan(sd_value)
        if add_star
            s = sprintf('$\\mathbf{%.3f}^{*}$ (--)', mean_value);
        else
            s = sprintf('$\\mathbf{%.3f}$ (--)', mean_value);
        end
    else
        if add_star
            s = sprintf('$\\mathbf{%.3f}^{*}$ (%.3f)', mean_value, sd_value);
        else
            s = sprintf('$\\mathbf{%.3f}$ (%.3f)', mean_value, sd_value);
        end
    end
end
