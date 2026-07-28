%% best_method_paired_ttest_holm_for_tables.m
% =========================================================================
% Purpose:
%   Identify the best-performing method for each metric and test whether
%   the best method is significantly better than strong baselines using
%   paired two-sided t-tests + Holm/Bonferroni correction.
%
% Important:
%   This optimized version strictly extracts the 10CV repeated-run vectors:
%       AUC       -> A_AUC_values
%       AUPR      -> m_A_AUPR_values
%       Precision -> Precision_values
%
%   These vectors match the mean ± SD values reported in Tables 4/7/8.
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

test_name    = '10CV';

correction_method = 'holm';

alpha = 0.05;

compare_only_strong_baselines = false;

strong_baselines = {'MSBMF', 'HGIMC', 'BNNR'};

require_significant_against_all = true;

%% ===================== Method files =====================

files = struct();

%Fdataset
dataset_name = 'Fdataset';      % 可改为 'Cdataset' 或 'CTDdataset2023'
files.multiGMF = fullfile(data_dir, 'Fdataset_multiGMF_10CV_fold_results.mat');
files.DDASKF   = fullfile(data_dir, 'Fdataset_STresult_DDASKF_10CV_fold_results.mat');
files.MSBMF    = fullfile(data_dir, 'Fdataset_STresult_MSBMF_10CV_fold_results.mat');
files.HGIMC    = fullfile(data_dir, 'Fdataset_STresult_HGIMC_10CV_fold_results.mat');
files.BNNR     = fullfile(data_dir, 'Fdataset_STresult_BNNR_10CV_fold_results.mat');
files.MBiRW    = fullfile(data_dir, 'Fdataset_STresult_MBiRW_10CV_fold_results.mat');


% %Cdataset
% dataset_name = 'Cdataset';
% 
% files.multiGMF = fullfile(data_dir, 'Cdataset_multiGMF_10CV_fold_results.mat');
% files.DDASKF   = fullfile(data_dir, 'Cdataset_STresult_DDASKF_10CV_fold_results.mat');
% files.MSBMF    = fullfile(data_dir, 'Cdataset_STresult_MSBMF_10CV_fold_results.mat');
% files.HGIMC    = fullfile(data_dir, 'Cdataset_STresult_HGIMC_10CV_fold_results.mat');
% files.BNNR     = fullfile(data_dir, 'Cdataset_STresult_BNNR_10CV_fold_results.mat');
% files.MBiRW    = fullfile(data_dir, 'Cdataset_STresult_MBiRW_10CV_fold_results.mat');

%%CTDdataset2023
% dataset_name = 'CTDdataset2023';

% files.multiGMF = fullfile(data_dir, 'CTDdataset2023_multiGMF_10CV_fold_results.mat');
% files.DDASKF   = fullfile(data_dir, 'CTDdataset2023_STresult_DDASKF_10CV_fold_results.mat');
% files.MSBMF    = fullfile(data_dir, 'CTDdataset2023_STresult_MSBMF_10CV_fold_results.mat');
% files.HGIMC    = fullfile(data_dir, 'CTDdataset2023_STresult_HGIMC_10CV_fold_results.mat');
% files.BNNR     = fullfile(data_dir, 'CTDdataset2023_STresult_BNNR_10CV_fold_results.mat');
% files.MBiRW    = fullfile(data_dir, 'CTDdataset2023_STresult_MBiRW_10CV_fold_results.mat');


metric_names = {'AUC', 'AUPR', 'Precision'};

%% ===================== Load data =====================

method_names = fieldnames(files);
Data = struct();

fprintf('\n============================================================\n');
fprintf('Dataset: %s | Test: %s | Correction: %s\n', ...
    dataset_name, test_name, correction_method);
fprintf('Vector source: A_AUC_values / m_A_AUPR_values / Precision_values\n');
fprintf('============================================================\n');

for i = 1:numel(method_names)

    mname = method_names{i};
    fpath = files.(mname);

    if ~isfile(fpath)
        error('File not found for method %s:\n%s', mname, fpath);
    end

    fprintf('\n[LOAD] %s\n%s\n', mname, fpath);

    S = load(fpath);

    for k = 1:numel(metric_names)

        metric = metric_names{k};
        vec = extract_metric_vector_strict(S, metric, mname);
        vec = vec(:);

        if isempty(vec)
            error('Cannot extract %s from %s.', metric, fpath);
        end

        Data.(mname).(metric) = vec;

        fprintf('  %-10s: n = %d, mean = %.6f, SD = %.6f\n', ...
            metric, numel(vec), mean(vec, 'omitnan'), std(vec, 'omitnan'));
    end
end

%% ===================== Best-method paired t-test =====================

all_results = table();
star_summary = table();

fprintf('\n\n==================== Best method significance suggestions ====================\n');

for k = 1:numel(metric_names)

    metric = metric_names{k};

    % ---------- compute mean and SD for each method ----------
    means = nan(numel(method_names), 1);
    sds   = nan(numel(method_names), 1);
    ns    = nan(numel(method_names), 1);

    for i = 1:numel(method_names)
        mname = method_names{i};
        x = Data.(mname).(metric);
        means(i) = mean(x, 'omitnan');
        sds(i)   = std(x, 'omitnan');
        ns(i)    = sum(~isnan(x));
    end

    % ---------- identify best method ----------
    [best_mean, best_idx] = max(means);
    best_method = method_names{best_idx};
    best_vec = Data.(best_method).(metric);

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
        comp_vec = Data.(cname).(metric);

        [x2, y2] = align_pair_vectors(best_vec, comp_vec);

        valid = ~isnan(x2) & ~isnan(y2);
        x2 = x2(valid);
        y2 = y2(valid);

        n_pair(j) = numel(x2);
        compared_mean(j) = mean(y2, 'omitnan');
        compared_sd(j)   = std(y2, 'omitnan');
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
            adj_p = min(raw_p * numel(raw_p), 1);
        otherwise
            error('Unknown correction method: %s', correction_method);
    end

    %
    sig_each = (adj_p < alpha) & (mean_diff > 0);

    if isempty(sig_each)
        add_star = false;
    else
        if require_significant_against_all
            add_star = all(sig_each);
        else
            add_star = any(sig_each);
        end
    end

    % ---------- print summary ----------
    fprintf('\nMetric: %s\n', metric);
    fprintf('  Best method: %s | mean = %.6f | SD = %.6f | n = %d\n', ...
        best_method, best_mean, sds(best_idx), ns(best_idx));

    if add_star
        fprintf('  Suggested table mark: ADD $^{*}$ to %s for %s\n', ...
            best_method, metric);
    else
        fprintf('  Suggested table mark: no $^{*}$ for this row\n');
    end

    for j = 1:numel(compared_methods)
        fprintf('    vs %-8s: compared_mean = %.6f, mean_diff = %.6f, p_raw = %.4g, p_adj = %.4g, sig = %d\n', ...
            compared_methods{j}, compared_mean(j), mean_diff(j), raw_p(j), adj_p(j), sig_each(j));
    end

    % ---------- save detailed results ----------
    T = table();
    T.Dataset       = repmat({dataset_name}, numel(compared_methods), 1);
    T.Test          = repmat({test_name}, numel(compared_methods), 1);
    T.Metric        = repmat({metric}, numel(compared_methods), 1);
    T.BestMethod    = repmat(best_method, numel(compared_methods), 1);
    T.BestMean      = repmat(best_mean, numel(compared_methods), 1);
    T.BestSD        = repmat(sds(best_idx), numel(compared_methods), 1);
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
    Srow.BestSD      = sds(best_idx);
    Srow.AddStar     = add_star;
    Srow.TableString = {make_latex_value(best_mean, sds(best_idx), add_star)};

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

function vec = extract_metric_vector_strict(S, metric, method_name)
% Strict extraction of repeated-run 10CV metric vectors.
%
% Correct variables:
%   AUC       -> A_AUC_values
%   AUPR      -> m_A_AUPR_values
%   Precision -> Precision_values
%
% These match the reported 10CV mean ± SD in the manuscript.

    switch lower(metric)

        case 'auc'
            if isfield(S, 'A_AUC_values')
                vec = S.A_AUC_values;
            else
                error('[%s] Missing A_AUC_values for AUC.', method_name);
            end

        case 'aupr'
            if isfield(S, 'm_A_AUPR_values')
                vec = S.m_A_AUPR_values;
            else
                error('[%s] Missing m_A_AUPR_values for AUPR.', method_name);
            end

        case 'precision'
            if isfield(S, 'Precision_values')
                vec = S.Precision_values;
            else
                error('[%s] Missing Precision_values for Precision.', method_name);
            end

        otherwise
            error('Unknown metric: %s', metric);
    end

    vec = vec(:);

    % Basic sanity check
    if numel(vec) ~= 10
        warning('[%s] %s vector length is %d, expected 10 repeated 10CV values.', ...
            method_name, metric, numel(vec));
    end
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

    if add_star
        s = sprintf('$\\mathbf{%.3f}^{*}$ (%.3f)', mean_value, sd_value);
    else
        s = sprintf('$\\mathbf{%.3f}$ (%.3f)', mean_value, sd_value);
    end
end
