function results = run_ablation_para(varargin)
    fprintf('========== 消融实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
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

    variants = {
        'proposed', 'Proposed - cSA-GOA + Multi-Subpop + Pareto Leader';
        'no_pareto', 'w/o Pareto Leader - 使用标量最优引导';
        'no_subpop', 'w/o Multi-Subpop - 单种群搜索';
        'no_goa', 'w/o GOA - 去掉排斥机制，纯随机游走'
    };

    reference_point = [1.0, 100000];

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for v_idx = 1:size(variants, 1)
        variant_name = variants{v_idx, 1};
        variant_desc = variants{v_idx, 2};
        fprintf('\n--- 运行变体: %s ---\n', variant_desc);

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
        temp_pareto_fronts = cell(n_runs, 1);

        if isempty(gcp('nocreate'))
            fprintf('启动并行池 (%d workers)...\n', n_workers);
            parpool('local', n_workers);
        end

        parfor run = 1:n_runs
            stream = RandStream('mt19937ar', 'Seed', run * 100 + v_idx * 1000);
            RandStream.setGlobalStream(stream);

            best_fit = 0;
            bestUAV = zeros(N_UAV, 2);
            cg_curve = zeros(1, 300);
            pareto_archive = [];

            switch variant_name
                case 'proposed'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'proposed');
                case 'no_pareto'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'no_pareto');
                case 'no_subpop'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'no_subpop');
                case 'no_goa'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'no_goa');
            end

            best_fits(run) = best_fit;

            center_point = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
            energies(run) = sum(params.k_move * fly_dist);

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
                temp_pareto_fronts{run} = pareto_front;

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
                temp_pareto_fronts{run} = [];
            end
        end

        pareto_fronts{v_idx} = temp_pareto_fronts;

        results.(variant_name) = struct();
        results.(variant_name).description = variant_desc;
        results.(variant_name).best_fits = best_fits;
        results.(variant_name).energies = energies;
        results.(variant_name).cov_high = cov_high;
        results.(variant_name).cov_total = cov_total;
        results.(variant_name).iter_counts = iter_counts;
        results.(variant_name).convergence_curves = convergence_curves;
        results.(variant_name).mean_fitness = mean(best_fits);
        results.(variant_name).std_fitness = std(best_fits);
        results.(variant_name).mean_energy = mean(energies);
        results.(variant_name).mean_cov_high = mean(cov_high);
        results.(variant_name).mean_cov_total = mean(cov_total);
        results.(variant_name).hv_values = hv_values;
        results.(variant_name).mean_hv = mean(hv_values);
        results.(variant_name).std_hv = std(hv_values);
        results.(variant_name).igd_values = igd_values;
        results.(variant_name).mean_igd = mean(igd_values);
        results.(variant_name).spread_values = spread_values;
        results.(variant_name).mean_spread = mean(spread_values);
        results.(variant_name).pareto_fronts = temp_pareto_fronts;
        results.(variant_name).mean_pareto_size = mean(pareto_sizes);

        fprintf('  >> %s 平均结果:\n', variant_desc);
        fprintf('     平均适应度: %.2f +/- %.2f\n', mean(best_fits), std(best_fits));
        fprintf('     平均能耗: %.2f J\n', mean(energies));
        fprintf('     平均高优先级覆盖率: %.2f%%\n', mean(cov_high));
        fprintf('     平均全局覆盖率: %.2f%%\n', mean(cov_total));
        fprintf('     平均HV: %.4f +/- %.4f\n', mean(hv_values), std(hv_values));
        fprintf('     平均Spread: %.4f\n', mean(spread_values));
    end

    fprintf('\n========== 开始执行消融实验 Min-Max 归一化 ==========\n');
    all_raw_points = [];
    var_keys = fieldnames(results);
    valid_vars = {};
    for i = 1:length(var_keys)
        v = var_keys{i};
        if ~strcmp(v, 'map_name') && ~strcmp(v, 'N_User') && ~strcmp(v, 'N_UAV')
            valid_vars{end+1} = v;
            for run = 1:length(results.(v).pareto_fronts)
                pf = results.(v).pareto_fronts{run};
                if ~isempty(pf)
                    all_raw_points = [all_raw_points; pf];
                end
            end
        end
    end

    max_U = sum(priorities); min_U = 0;
    max_L = N_User * params.max_latency; min_L = 0;
    max_E = params.E_max * N_UAV; min_E = 0;

    all_points_norm = [];
    for i = 1:length(valid_vars)
        v = valid_vars{i};
        for run = 1:length(results.(v).pareto_fronts)
            pf = results.(v).pareto_fronts{run};
            if ~isempty(pf)
                pf_norm = zeros(size(pf));
                pf_norm(:, 1) = (max_U - pf(:, 1)) / (max_U - min_U + 1e-6);
                pf_norm(:, 2) = (pf(:, 2) - min_L) / (max_L - min_L + 1e-6);
                pf_norm(:, 3) = (pf(:, 3) - min_E) / (max_E - min_E + 1e-6);
                all_points_norm = [all_points_norm; pf_norm];
                results.(v).pareto_fronts_norm{run} = pf_norm;
            else
                results.(v).pareto_fronts_norm{run} = [];
            end
        end
    end

    all_points_norm = unique(all_points_norm, 'rows');
    true_front_norm = extractNonDominated(all_points_norm);
    ref_point_norm = [1.1, 1.1, 1.1];

    fprintf('\n========== 消融实验最终汇总表格 (Min-Max 归一化) ==========\n');
    fprintf('%-12s | %-10s | %-10s | %-13s | %-13s | %-8s\n', ...
        'Variant', 'Fitness', 'Energy(J)', 'HV(Norm)↑', 'IGD(Norm)↓', 'Spread↑');
    fprintf('%s\n', repmat('-', 1, 95));

    for i = 1:length(valid_vars)
        v = valid_vars{i};
        r = results.(v);
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
        results.(v).mean_hv_norm = nanmean(hvs); results.(v).mean_igd_norm = nanmean(igds); results.(v).mean_spread_norm = nanmean(spreads);

        fprintf('%-12s | %-10.2f | %-10.2f | %-5.4f±%-5.4f | %-5.4f±%-5.4f | %-8.4f\n', ...
            v, r.mean_fitness, r.mean_energy, nanmean(hvs), nanstd(hvs), nanmean(igds), nanstd(igds), nanmean(spreads));
    end
    fprintf('%s\n', repmat('-', 1, 95));

    results_file = fullfile(project_dir, 'experiments', ['ablation_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results', 'map_data', 'params');
    fprintf('\n完美消融结果已保存至: %s\n', results_file);
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos,1)
        dists_uav = sqrt(sum((UAV_pos - User_pos(i,:)).^2, 2));
        covered_by_uav = any(dists_uav <= UAV_radius);

        if size(RRH,1) > 0
            dists_rrh = sqrt(sum((RRH - User_pos(i,:)).^2, 2));
            covered_by_rrh = any(dists_rrh <= RRH_radius);
        else
            covered_by_rrh = false;
        end

        if covered_by_uav || covered_by_rrh
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos,1);
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