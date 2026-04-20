function results = run_comparison_para(varargin)
    fprintf('========== 对比实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'comparison_algorithms')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
    addParameter(p, 'verbose', false);
    addParameter(p, 'n_workers', 3);
    parse(p, varargin{:});
    n_runs = p.Results.n_runs;
    map_name = p.Results.map_name;
    n_workers = p.Results.n_workers;

    maps_dir = fullfile(project_dir, 'maps');
    map_file = fullfile(maps_dir, [map_name, '.mat']);

    if ~exist(map_file, 'file')
        error('地图文件不存在: %s', map_file);
    end

    fprintf('--- 加载固定地图: %s ---\n', map_name);
    map_data = load(map_file);
    User = map_data.User;
    priorities = map_data.priorities;
    N_User = map_data.N_User;

    N_User_actual = size(User, 1);
    if N_User_actual == 200
        N_UAV = 8;
        config_suffix = 'Small';
    elseif N_User_actual == 500
        N_UAV = 15;
        config_suffix = 'Medium';
    else
        N_UAV = 25;
        config_suffix = 'Large';
    end

    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);

    config_file = fullfile(maps_dir, ['Map_', config_suffix, '_Config.mat']);
    if exist(config_file, 'file')
        config_data = load(config_file);
        RRH = config_data.RRH;
        RRH_type = config_data.RRH_type;
        N_RRH = config_data.N_RRH;
        N_eRRH = config_data.N_eRRH;
        UAV_type = config_data.UAV_type;
        params.D = config_data.D;
        params.C = config_data.C;
        params.DT = config_data.DT;
    else
        Lb = [0, 0];
        Ub = [1000, 1000];
        N_RRH = 10;
        N_eRRH = 4;
        RRH = Lb + (Ub - Lb) .* rand(N_RRH, 2);
        RRH_type = zeros(N_RRH, 1);
        RRH_type(1:N_eRRH) = 1;
        UAV_type = zeros(N_UAV, 1);
        UAV_type(1:floor(N_UAV * 0.3)) = 1;
        params.D = [];
        params.C = [];
        params.DT = [];
    end

    params.RRH = RRH;

    Lb = [0, 0];
    Ub = [1000, 1000];

    params.max_latency = 1.0;
    params.B_total = 10e6;
    params.F_total = 10e9;
    params.noise = 1e-13;
    params.P_tx = 0.5;
    params.PtxU = 0.5;
    params.B_total_relay = 20e6;
    params.f_BBU = 50e9;
    params.kappa = 1e-28;
    params.k_move = 10;
    params.cover_radius = 150;
    params.RRH_radius = 150;
    params.E_max = 50000;
    params.D_UU = 10;
    params.D_RU = 10;

    if ~isfield(params, 'D') || isempty(params.D) || length(params.D) ~= N_User
        params.D = ones(N_User, 1) * 2e6;
    end
    if ~isfield(params, 'C') || isempty(params.C) || length(params.C) ~= N_User
        params.C = ones(N_User, 1) * 0.5e9;
    end

    params.FES_max = 300;
    params.K = 40;

    params.G_weights = [0.4, 0.3, 0.3];
    params.subpop_params = struct();
    params.subpop_params.mu0 = [500, 500];
    params.subpop_params.sigma0 = [150, 150; 120, 120; 80, 80];
    params.subpop_params.sigma_min = [5, 8, 3];
    params.subpop_params.w_inertia = [0.7, 0.6, 0.8];
    params.subpop_params.c = [0.15, 0.10, 0.08];
    params.subpop_params.q = [0.6, 0.5, 0.4];
    params.subpop_params.beta = [0.8, 0.7, 0.6];

    hyperparam_file = fullfile(project_dir, 'experiments', 'best_algo_hyperparams.mat');
    if exist(hyperparam_file, 'file')
        opt_params = load(hyperparam_file);
        params.K = opt_params.best_K;
        params.subpop_params.q = [opt_params.best_q, opt_params.best_q, opt_params.best_q];
        fprintf('🎯 [动态注入] 已加载最优超参数 K=%d, q=%.1f\n', params.K, opt_params.best_q);
    end

    params.enable_early_stop = false;
    params.enable_smart_stop = false;

    algorithms = {
        'cSA_GOA', 'cSA-GOA (Proposed)';
        'PSO', 'PSO (Particle Swarm Optimization)';
        'GA', 'GA (Genetic Algorithm)';
        'GOA', 'GOA (Grasshopper Optimization)';
        'cSA', 'cSA (Compact Sine Algorithm)';
        'GWO', 'GWO (Grey Wolf Optimizer)'
    };

    reference_point = [1.0, 100000];

    if isempty(gcp('nocreate'))
        fprintf('启动并行池 (%d workers)...\n', n_workers);
        parpool('local', n_workers);
    end

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for alg_idx = 1:size(algorithms, 1)
        alg_name = algorithms{alg_idx, 1};
        alg_desc = algorithms{alg_idx, 2};
        fprintf('\n--- 运行算法: %s (并行) ---\n', alg_desc);

        best_fits = zeros(n_runs, 1);
        energies = zeros(n_runs, 1);
        cov_high = zeros(n_runs, 1);
        cov_total = zeros(n_runs, 1);
        iter_counts = zeros(n_runs, 1);
        convergence_curves = cell(n_runs, 1);
        hv_values = zeros(n_runs, 1);
        igd_values = zeros(n_runs, 1);
        spread_values = zeros(n_runs, 1);
        pareto_sizes = zeros(n_runs, 1);
        pareto_fronts_cell = cell(n_runs, 1);
        runtimes = zeros(n_runs, 1);
        success_rates = zeros(n_runs, 1);

        parfor run = 1:n_runs
            stream = RandStream('mt19937ar', 'Seed', run * 200 + alg_idx * 1000);
            RandStream.setGlobalStream(stream);

            t_start = tic;

            best_fit = 0;
            bestUAV = zeros(N_UAV, 2);
            cg_curve = zeros(1, 300);
            pareto_archive = [];

            switch alg_name
                case 'cSA_GOA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, ~, ~, ~, ~, ~, pareto_archive] = ...
                        cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'PSO'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        PSO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GOA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GOA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'cSA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GWO'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GWO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
            end

            runtimes(run) = toc(t_start);

            [~, ~, ~, s_rate] = calcMEC_Objectives(bestUAV, User, priorities, params);
            success_rates(run) = s_rate;

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                arch_util = [pareto_archive.Utility];
                arch_lat = [pareto_archive.Latency];
                arch_energy = [pareto_archive.Energy];

                max_U = sum(priorities);
                max_L = N_User * params.max_latency;
                max_E = params.E_max * N_UAV;

                norm_util = (max_U - arch_util) / (max_U + 1e-6);
                norm_lat = arch_lat / (max_L + 1e-6);
                norm_eng = arch_energy / (max_E + 1e-6);

                W_U = 0.70; W_L = 0.15; W_E = 0.15;
                distances_to_ideal = sqrt(W_U * norm_util.^2 + W_L * norm_lat.^2 + W_E * norm_eng.^2);
                [~, idx_knee] = min(distances_to_ideal);

                best_fits(run) = best_fit;
                energies(run) = pareto_archive(idx_knee).Energy;

                if isfield(pareto_archive(idx_knee), 'UAV_pos')
                    bestUAV = pareto_archive(idx_knee).UAV_pos;
                end
            else
                best_fits(run) = best_fit;
                center_point = repmat([500, 500], N_UAV, 1);
                fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
                energies(run) = sum(params.k_move * fly_dist);
            end

            cov_high(run) = calcCoverageWithRRH(bestUAV, User(priorities>=3,:), params.cover_radius, RRH, params.RRH_radius) * 100;
            cov_total(run) = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius) * 100;

            iter_counts(run) = length(cg_curve);
            convergence_curves{run} = cg_curve;

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                pareto_front = zeros(length(pareto_archive), 3);
                for p_idx = 1:length(pareto_archive)
                    pareto_front(p_idx, 1) = pareto_archive(p_idx).Utility;
                    pareto_front(p_idx, 2) = pareto_archive(p_idx).Latency;
                    pareto_front(p_idx, 3) = pareto_archive(p_idx).Energy;
                end
                pareto_sizes(run) = length(pareto_archive);
                pareto_fronts_cell{run} = pareto_front;

                try
                    max_U = sum(priorities); min_U = 0;
                    max_L = N_User * params.max_latency; min_L = 0;
                    max_E = params.E_max * N_UAV; min_E = 0;

                    pf_norm = zeros(size(pareto_front));
                    pf_norm(:, 1) = (max_U - pareto_front(:, 1)) / (max_U - min_U + 1e-6);
                    pf_norm(:, 2) = (pareto_front(:, 2) - min_L) / (max_L - min_L + 1e-6);
                    pf_norm(:, 3) = (pareto_front(:, 3) - min_E) / (max_E - min_E + 1e-6);

                    ref_point_norm = [1.1, 1.1, 1.1];
                    metrics = calculate_all_metrics(pf_norm, [], ref_point_norm);
                    hv_values(run) = metrics.hv;
                    spread_values(run) = metrics.spread;
                catch
                    hv_values(run) = 0;
                    spread_values(run) = NaN;
                end
            else
                pareto_sizes(run) = 0;
                hv_values(run) = 0;
                spread_values(run) = NaN;
                pareto_fronts_cell{run} = [];
            end
        end

        pareto_fronts{alg_idx} = pareto_fronts_cell;

        results.(alg_name) = struct();
        results.(alg_name).description = alg_desc;
        results.(alg_name).best_fits = best_fits;
        results.(alg_name).energies = energies;
        results.(alg_name).cov_high = cov_high;
        results.(alg_name).cov_total = cov_total;
        results.(alg_name).iter_counts = iter_counts;
        results.(alg_name).convergence_curves = convergence_curves;
        results.(alg_name).mean_fitness = mean(best_fits);
        results.(alg_name).std_fitness = std(best_fits);
        results.(alg_name).mean_energy = mean(energies);
        results.(alg_name).mean_cov_high = mean(cov_high);
        results.(alg_name).mean_cov_total = mean(cov_total);
        results.(alg_name).hv_values = hv_values;
        results.(alg_name).mean_hv = mean(hv_values);
        results.(alg_name).std_hv = std(hv_values);
        results.(alg_name).igd_values = igd_values;
        results.(alg_name).mean_igd = mean(igd_values);
        results.(alg_name).spread_values = spread_values;
        results.(alg_name).mean_spread = mean(spread_values);
        results.(alg_name).mean_pareto_size = mean(pareto_sizes);
        results.(alg_name).pareto_fronts = pareto_fronts_cell;
        results.(alg_name).runtimes = runtimes;
        results.(alg_name).mean_runtime = mean(runtimes);
        results.(alg_name).success_rates = success_rates;
        results.(alg_name).mean_success_rate = mean(success_rates);

        fprintf('  >> %s 平均结果:\n', alg_desc);
        fprintf('     Mean Fitness: %.2f +/- %.2f\n', mean(best_fits), std(best_fits));
        fprintf('     Mean Energy: %.2f J\n', mean(energies));
        fprintf('     Mean High-Priority Coverage: %.2f%%\n', mean(cov_high));
        fprintf('     Mean HV: %.4f +/- %.4f\n', mean(hv_values), std(hv_values));
        fprintf('     Mean Runtime: %.2f s\n', mean(runtimes));
        fprintf('     Mean Success Rate: %.2f%%\n', mean(success_rates) * 100);
    end

    fprintf('\n========== 开始执行全局 Min-Max 归一化与指标修正 ==========\n');
    all_raw_points = [];
    alg_names = fieldnames(results);
    valid_algs = {};
    for i = 1:length(alg_names)
        alg = alg_names{i};
        if ~strcmp(alg, 'map_name') && ~strcmp(alg, 'N_User') && ~strcmp(alg, 'N_UAV')
            valid_algs{end+1} = alg;
            for run = 1:length(results.(alg).pareto_fronts)
                pf = results.(alg).pareto_fronts{run};
                if ~isempty(pf)
                    all_raw_points = [all_raw_points; pf];
                end
            end
        end
    end

    max_U = sum(priorities); min_U = 0;
    max_L = N_User * params.max_latency; min_L = 0;
    max_E = params.E_max * N_UAV; min_E = 0;

    fprintf('效用理论边界: [%.2f, %.2f]\n', min_U, max_U);
    fprintf('时延物理边界: [%.2f, %.2f] s\n', min_L, max_L);
    fprintf('能耗物理边界: [%.2f, %.2f] J\n', min_E, max_E);

    all_points_norm = [];
    for i = 1:length(valid_algs)
        alg = valid_algs{i};
        for run = 1:length(results.(alg).pareto_fronts)
            pf = results.(alg).pareto_fronts{run};
            if ~isempty(pf)
                pf_norm = zeros(size(pf));
                pf_norm(:, 1) = (max_U - pf(:, 1)) / (max_U - min_U + 1e-6);
                pf_norm(:, 2) = (pf(:, 2) - min_L) / (max_L - min_L + 1e-6);
                pf_norm(:, 3) = (pf(:, 3) - min_E) / (max_E - min_E + 1e-6);
                all_points_norm = [all_points_norm; pf_norm];
                results.(alg).pareto_fronts_norm{run} = pf_norm;
            else
                results.(alg).pareto_fronts_norm{run} = [];
            end
        end
    end

    all_points_norm = unique(all_points_norm, 'rows');
    true_front_norm = extractNonDominated(all_points_norm);
    ref_point_norm = [1.1, 1.1, 1.1];

    fprintf('\n========== 对比实验最终汇总表格 (Min-Max 归一化) ==========\n');
    fprintf('%-12s | %-10s | %-10s | %-13s | %-13s | %-8s\n', ...
        'Algorithm', 'Fitness', 'Energy(J)', 'HV(Norm)↑', 'IGD(Norm)↓', 'Spread↑');
    fprintf('%s\n', repmat('-', 1, 95));

    for i = 1:length(valid_algs)
        alg = valid_algs{i};
        r = results.(alg);
        hvs = zeros(n_runs, 1); igds = zeros(n_runs, 1); spreads = zeros(n_runs, 1);
        for run = 1:length(r.pareto_fronts_norm)
            pf_norm = r.pareto_fronts_norm{run};
            if ~isempty(pf_norm)
                try
                    metrics = calculate_all_metrics(pf_norm, true_front_norm, ref_point_norm);
                    hvs(run) = metrics.hv; igds(run) = metrics.igd; spreads(run) = metrics.spread;
                catch
                    hvs(run) = NaN; igds(run) = NaN; spreads(run) = NaN;
                end
            end
        end
        results.(alg).mean_hv_norm = nanmean(hvs); results.(alg).std_hv_norm = nanstd(hvs);
        results.(alg).mean_igd_norm = nanmean(igds); results.(alg).std_igd_norm = nanstd(igds);
        results.(alg).mean_spread_norm = nanmean(spreads); results.(alg).std_spread_norm = nanstd(spreads);

        fprintf('%-12s | %-10.2f | %-10.2f | %-5.4f±%-5.4f | %-5.4f±%-5.4f | %-8.4f\n', ...
            alg, r.mean_fitness, r.mean_energy, nanmean(hvs), nanstd(hvs), nanmean(igds), nanstd(igds), nanmean(spreads));
    end
    fprintf('%s\n注: HV↑越大越好 | IGD↓越小越好 | Spread↑分布越均匀\n', repmat('-', 1, 95));

    results_file = fullfile(project_dir, 'experiments', ['comparison_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results', 'map_data', 'params');
    fprintf('\n完美结果已保存至: %s\n', results_file);
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos, 1)
        user = User_pos(i, :);
        dist_UAV = min(sqrt(sum((UAV_pos - repmat(user, size(UAV_pos, 1), 1)).^2, 2)));
        dist_RRH = min(sqrt(sum((RRH - repmat(user, size(RRH, 1), 1)).^2, 2)));
        if dist_UAV <= UAV_radius || dist_RRH <= RRH_radius
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos, 1);
end

function pf = extractNonDominated(points)
    n = size(points, 1);
    is_dominated = false(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j
                if all(points(j, :) <= points(i, :)) && any(points(j, :) < points(i, :))
                    is_dominated(i) = true;
                    break;
                end
            end
        end
    end
    pf = points(~is_dominated, :);
end