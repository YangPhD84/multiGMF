%% inspect_fold_result_mat_files.m
clear; clc;

data_dir = 'C:\Users\biny9\Desktop\MultiGMF_Rresponses\t-test';

files = {
    'Fdataset_best_multiGMF_10CV_fold_results.mat'
    'Fdataset_STresult_HGIMC_10CV_fold_results.mat'
    'Fdataset_STresult_DDASKF_10CV_fold_results.mat'
    'Fdataset_STresult_BNNR_10CV_fold_results.mat'
    'Fdataset_STresult_MBiRW_10CV_fold_results.mat'
    'Fdataset_STresult_MSBMF_10CV_fold_results.mat'
};

for f = 1:numel(files)

    mat_path = fullfile(data_dir, files{f});

    fprintf('\n\n============================================================\n');
    fprintf('[FILE] %s\n', mat_path);
    fprintf('============================================================\n');

    if ~isfile(mat_path)
        warning('File not found: %s', mat_path);
        continue;
    end

    S = load(mat_path);
    print_struct_info(S, '');
end


function print_struct_info(S, prefix)

    names = fieldnames(S);

    for i = 1:numel(names)

        name = names{i};
        value = S.(name);

        full_name = name;
        if ~isempty(prefix)
            full_name = [prefix '.' name];
        end

        if isnumeric(value)

            sz = size(value);
            value_vec = value(:);
            value_vec = value_vec(~isnan(value_vec));

            fprintf('%-50s | numeric | size = %-15s', ...
                full_name, mat2str(sz));

            if ~isempty(value_vec)
                fprintf(' | min = %.6f | mean = %.6f | max = %.6f | SD = %.6f', ...
                    min(value_vec), mean(value_vec), max(value_vec), std(value_vec));
            end

            fprintf('\n');

        elseif isstruct(value)

            fprintf('%-50s | struct  | size = %s\n', full_name, mat2str(size(value)));

            if numel(value) == 1
                print_struct_info(value, full_name);
            else
                fprintf('%-50s | struct array, skipped deep print\n', full_name);
            end

        elseif istable(value)

            fprintf('%-50s | table   | size = %s\n', full_name, mat2str(size(value)));
            disp(value.Properties.VariableNames');

        elseif iscell(value)

            fprintf('%-50s | cell    | size = %s\n', full_name, mat2str(size(value)));

        else

            fprintf('%-50s | %s | size = %s\n', ...
                full_name, class(value), mat2str(size(value)));
        end
    end
end