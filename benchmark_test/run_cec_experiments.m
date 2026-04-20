function results = run_cec_experiments(functions, Dim, max_fes, runs)
    if nargin < 1 || isempty(functions)
        functions = 1:30;
    end
    if nargin < 2 || isempty(Dim)
        Dim = 30;
    end
    if nargin < 3 || isempty(max_fes)
        max_fes = Dim * 10000;
    end
    if nargin < 4 || isempty(runs)
        runs = 30;
    end

    original_dir = pwd;
    cec_base_dir = fullfile(original_dir, 'cec2017_data', 'cec2017');
    if ~exist(cec_base_dir, 'dir')
        error('Cannot find CEC data folder: %s', cec_base_dir);
    end
    cd(cec_base_dir);
    addpath(cec_base_dir);

    Lb = -100 * ones(1, Dim);
    Ub = 100 * ones(1, Dim);
    results = struct();

    for f = functions
        fprintf('Starting F%d...\n', f);
        run_fits = zeros(runs, 1);

        parfor r = 1:runs
            [best_fit, ~, ~] = main_cec(f, Dim, max_fes, Lb, Ub);
            run_fits(r) = best_fit;
        end

        results(f).mean = mean(run_fits);
        results(f).std = std(run_fits);
        results(f).raw = run_fits;
        fprintf('F%d Done. Mean: %.4e\n', f, results(f).mean);
    end

    cd(original_dir);
    save(['results_cec2017_D', num2str(Dim), '.mat'], 'results');
end