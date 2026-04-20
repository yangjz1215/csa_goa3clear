function results = run_ablation(varargin)
    fprintf('========== 消融实验开始 ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
    parse(p, varargin{:});
    n_runs = p.Results.n_runs;
    map_name = p.Results.map_name;

    rng('shuffle');

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

    Lb = [0, 0];
    Ub = [1000, 1000];

    params.D_UU = 10;
    params.D_RU = 10;
    params.cover_radius = 150;
    params.RRH_radius = 150;
    params.E_max = 50000;
    params.k_move = 15;
    params.bandwidth = 1e6;

    params.Pho = 100;
    params.ki = 1e-27;
    params.PtxU = 10;
    params.PtxEU = 5;
    params.Ptx = 1;
    params.PtxR = 10;
    params.alpha0 = 1.42e-4;
    params.sigma2 = 3.98e-12;
    params.f_eUAV = 2e9;
    params.f_eRRH = 2e9;
    params.f_BBU = 4e9;

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

    params.enable_early_stop = false;
    params.enable_smart_stop = false;

    variants = {
        'proposed', 'Proposed - cSA-GOA + Multi-Subpop + Pareto Leader';
        'no_pareto', 'w/o Pareto Leader - 使用标量最优引导';
        'w/o_subpop', 'w/o Multi-Subpop - 单种群搜索';
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
        pareto_fronts{v_idx} = cell(n_runs, 1);

        for run = 1:n_runs
            fprintf('  Run %d/%d... ', run, n_runs);
            rng(run * 100 + v_idx * 1000);

            switch variant_name
                case 'proposed'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'proposed');
                case 'no_pareto'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'w/o_pareto_leader');
                case 'w/o_subpop'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'w/o_subpop');
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
                pareto_front = zeros(length(pareto_archive), 2);
                for p_idx = 1:length(pareto_archive)
                    pareto_front(p_idx, 1) = pareto_archive(p_idx).Coverage / N_User;
                    pareto_front(p_idx, 2) = pareto_archive(p_idx).Energy;
                end
                pareto_sizes(run) = length(pareto_archive);
                pareto_fronts{v_idx}{run} = pareto_front;

                try
                    metrics = calculate_all_metrics(pareto_front, [], reference_point);
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
                pareto_fronts{v_idx}{run} = [];
            end

            fprintf('Fitness=%.2f, CovHigh=%.2f%%, HV=%.4f\n', ...
                best_fit, cov_high(run), hv_values(run));
        end

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
        results.(variant_name).mean_pareto_size = mean(pareto_sizes);

        fprintf('  >> %s 平均结果:\n', variant_desc);
        fprintf('     平均适应度: %.2f +/- %.2f\n', mean(best_fits), std(best_fits));
        fprintf('     平均能耗: %.2f J\n', mean(energies));
        fprintf('     平均高优先级覆盖率: %.2f%%\n', mean(cov_high));
        fprintf('     平均全局覆盖率: %.2f%%\n', mean(cov_total));
        fprintf('     平均HV: %.4f +/- %.4f\n', mean(hv_values), std(hv_values));
        fprintf('     平均Spread: %.4f\n', mean(spread_values));
    end

    fprintf('\n========== 计算IGD (使用合并Pareto前沿作为参考) ==========\n');
    all_pareto_points = [];
    for v_idx = 1:size(variants, 1)
        for run = 1:n_runs
            pf = pareto_fronts{v_idx}{run};
            if ~isempty(pf)
                all_pareto_points = [all_pareto_points; pf];
            end
        end
    end
    
    if ~isempty(all_pareto_points)
        true_front = extractNonDominated(all_pareto_points);
        fprintf('  合并Pareto解数量: %d, 非支配解数量: %d\n', size(all_pareto_points, 1), size(true_front, 1));
        
        for v_idx = 1:size(variants, 1)
            variant_name = variants{v_idx, 1};
            igd_values = zeros(n_runs, 1);
            for run = 1:n_runs
                pf = pareto_fronts{v_idx}{run};
                if ~isempty(pf)
                    igd_values(run) = igd(pf, true_front);
                else
                    igd_values(run) = NaN;
                end
            end
            results.(variant_name).igd_values = igd_values;
            results.(variant_name).mean_igd = mean(igd_values);
        end
    end

    fprintf('\n========== 消融实验完成 (地图: %s) ==========\n', map_name);
    fprintf('\n========== 结果汇总表格 ==========\n');
    fprintf('%-35s | %-10s | %-10s | %-10s | %-10s | %-8s\n', ...
        '变体', '适应度', '能耗(J)', '高优%', '全局%', 'HV');
    fprintf('%s\n', repmat('-', 1, 95));
    for v_idx = 1:size(variants, 1)
        variant_name = variants{v_idx, 1};
        r = results.(variant_name);
        fprintf('%-35s | %-10.2f | %-10.2f | %-10.2f | %-10.2f | %-8.4f\n', ...
            variant_name, r.mean_fitness, r.mean_energy, ...
            r.mean_cov_high, r.mean_cov_total, r.mean_hv);
    end
    fprintf('================================\n');

    fprintf('\n========== 多目标优化指标 ==========\n');
    fprintf('%-35s | %-12s | %-12s | %-12s\n', ...
        '变体', 'HV(mean±std)', 'IGD(mean)', 'Spread(mean)');
    fprintf('%s\n', repmat('-', 1, 75));
    for v_idx = 1:size(variants, 1)
        variant_name = variants{v_idx, 1};
        r = results.(variant_name);
        fprintf('%-35s | %-12.4f | %-12.4f | %-12.4f\n', ...
            variant_name, r.mean_hv, r.mean_igd, r.mean_spread);
    end
    fprintf('================================\n');

    results_file = fullfile(project_dir, 'experiments', ['ablation_results_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results');
    fprintf('结果已保存: %s\n', results_file);
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