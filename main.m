% main.m - 单目标优化主程序
% 优化目标：最大化适应度（优先级覆盖之和）

clear; clc; close all;
rng('shuffle');  % 根据当前时间生成真随机数

%% 1. 场景参数配置
Lb = [0, 0];
Ub = [1000, 1000];

N_User = 500;
User = Lb + (Ub - Lb) .* rand(N_User, 2);

ratios = [0.25, 0.25, 0.25, 0.25];
num_levels = length(ratios);
counts = floor(ratios * N_User);
remainder = N_User - sum(counts);
extra_indices = randperm(num_levels, remainder);
for i = 1:remainder
    counts(extra_indices(i)) = counts(extra_indices(i)) + 1;
end
task_levels = [];
for i = 1:num_levels
    task_levels = [task_levels; i * ones(counts(i), 1)];
end
priorities = task_levels(randperm(N_User));
U_high_idx = find(priorities >= 3);

%% 1.5 创新点：长尾用户机会成本评估 (孤岛衰减机制)
% 计算每个用户的有效优先级 (Effective Priority)
center_point = [500, 500];
cover_radius = 150;  % 覆盖半径（与params保持一致）
effective_priorities = zeros(N_User, 1);

for i = 1:N_User
    dists_to_others = sqrt(sum((User - User(i,:)).^2, 2));
    neighbor_count = sum(dists_to_others <= cover_radius) - 1;

    dist_to_center = norm(User(i,:) - center_point);

    if neighbor_count == 0
        discount_factor = exp(-dist_to_center / 500);
    elseif neighbor_count < 3
        discount_factor = 0.8 + 0.2 * exp(-dist_to_center / 800);
    else
        discount_factor = 1.0;
    end

    effective_priorities(i) = priorities(i) * discount_factor;
end

priorities = effective_priorities;
U_high_idx = find(priorities >= 3);

N_RRH = 10;
N_eRRH = 4;
RRH = GenerateRRH(N_RRH, Ub, Lb);
RRH_type = zeros(N_RRH, 1);
RRH_type(1:N_eRRH) = 1;

N_UAV = 15;
N_eUAV = 5;
UAV_type = zeros(N_UAV, 1);
UAV_type(1:N_eUAV) = 1;

params = struct();
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

D_max = 500;
D_min = 100;
params.D = ((D_max-D_min)*rand(N_User,1)+D_min)*8192;

C_max = 0.05;
C_min = 0.01;
params.C = ((C_max-C_min)*rand(N_User,1)+C_min)*10^9;

DT_max = 1.2;
DT_min = 0.8;
params.DT = ((DT_max-DT_min)*rand(N_User,1)+DT_min);

params.B_total = 20e6;
params.F_total = 10e9;
params.max_latency = 1.0;
params.noise = 1e-13;
params.kappa = 1e-27;
params.P_tx = 1;
params.PtxU = 0.5;
params.f_BBU = 50e9;
params.k_move = 15;
params.enable_bilevel = true;
params.RRH = RRH;
params.B_total_relay = 5e6;

%% 2. 算法参数配置
optimization_method = 'fixed';

params.FES_max = 300;
params.K = 40;

hyperparam_file = 'best_algo_hyperparams.mat';
if exist(hyperparam_file, 'file')
    opt_params = load(hyperparam_file);
    params.K = opt_params.best_K;
    params.subpop_params.q = [opt_params.best_q, opt_params.best_q, opt_params.best_q];
    fprintf('🎯 [动态注入] 已加载最优超参数 K=%d, q=%.1f\n', params.K, opt_params.best_q);
elseif exist(fullfile('experiments', hyperparam_file), 'file')
    opt_params = load(fullfile('experiments', hyperparam_file));
    params.K = opt_params.best_K;
    params.subpop_params.q = [opt_params.best_q, opt_params.best_q, opt_params.best_q];
    fprintf('🎯 [动态注入] 已加载最优超参数 K=%d, q=%.1f\n', params.K, opt_params.best_q);
end

params.enable_early_stop = true;
params.enable_smart_stop = true;

switch optimization_method
    case 'fixed'
        params.G_weights = [0.35, 0.45, 0.2];
        params.subpop_params = struct(...
            'q', [0.5, 0.4, 0.3], ...
            'c', [0.3, 0.25, 0.2], ...
            'beta', [0.2, 0.3, 0.1], ...
            'sigma0', [25, 30, 20]);
        fprintf('使用固定参数方法\n');
end

%% 3. 运行算法
fprintf('\n========== 开始运行 cSA-GOA 算法 ==========\n');
[best_fit, bestUAV, cg_curve, energy_consumption, E_remaining_history, final_E_remaining, curr_curve, actual_iter, ~, pareto_archive] = ...
    cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, ...
    Ub, Lb, params, priorities);

%% 4. 结果可视化
figure('Name','收敛曲线','Position',[100,100,800,500]);
plot(1:actual_iter, curr_curve(1:actual_iter), 'b-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('当前适应度');
title('算法收敛曲线'); 
grid on;
xticks(0:10:actual_iter);
xlim([1, actual_iter]);

figure('Name','能耗变化','Position',[200,200,800,500]);
plot(1:actual_iter, energy_consumption(1:actual_iter), 'r-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('总能耗（J）');
title('无人机总能耗变化'); grid on;

% 计算能耗用于Pareto图
center_point = repmat([500, 500], N_UAV, 1);
fly_distances = sqrt(sum((bestUAV - center_point).^2, 2));
total_energy_consumed = params.k_move * sum(fly_distances);
final_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

% Pareto前沿分析 (3D)
if ~isempty(pareto_archive)
    arch_util = [pareto_archive.Utility];
    arch_lat = [pareto_archive.Latency];
    arch_energy = [pareto_archive.Energy];

    figure('Name','3D Pareto 前沿分析','Position',[400,200,800,600]);
    scatter3(arch_lat, arch_energy, arch_util, 80, arch_util, 'filled', 'MarkerEdgeColor', 'k');
    colormap(jet);
    colorbar;

    grid on; view(45, 30);
    title('Bilevel MEC 联合调度 3D Pareto 前沿');
    xlabel('总计算与传输延迟 (s) (越小越好)');
    ylabel('系统总能耗 (J) (越小越好)');
    zlabel('系统总效用 (越大越好)');

    norm_util = (max(arch_util) - arch_util) / (max(arch_util) - min(arch_util) + 1e-6);
    norm_lat = (arch_lat - min(arch_lat)) / (max(arch_lat) - min(arch_lat) + 1e-6);
    norm_eng = (arch_energy - min(arch_energy)) / (max(arch_energy) - min(arch_energy) + 1e-6);
    distances_to_ideal = sqrt(norm_util.^2 + norm_lat.^2 + norm_eng.^2);
    [~, idx_knee] = min(distances_to_ideal);

    hold on;
    scatter3(arch_lat(idx_knee), arch_energy(idx_knee), arch_util(idx_knee), ...
        250, 'pentagram', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    legend('Pareto 非支配解', 'Knee Point (最佳权衡解)', 'Location', 'best');

    fprintf('\n========== 3D Pareto 前沿分析 ==========\n');
    fprintf('Knee Point 解 - Utility: %.2f, Latency: %.4f s, Energy: %.2f J\n', ...
        arch_util(idx_knee), arch_lat(idx_knee), arch_energy(idx_knee));
    fprintf('=======================================\n');
end

figure('Name','最终部署','Position',[300,300,900,700]);
plot(User(priorities<3,1), User(priorities<3,2), 'bo','MarkerSize',6,'DisplayName','低优先级用户');
hold on;
plot(User(priorities>=3,1), User(priorities>=3,2), 'ro','MarkerSize',8,'DisplayName','高优先级用户');
plot(RRH(:,1), RRH(:,2), 'g^','MarkerSize',10,'DisplayName','RRH');
plot(bestUAV(UAV_type==0,1), bestUAV(UAV_type==0,2), 'k*','MarkerSize',12,'DisplayName','普通UAV');
plot(bestUAV(UAV_type==1,1), bestUAV(UAV_type==1,2), 'm*','MarkerSize',12,'DisplayName','增强UAV');
xlabel('X坐标（米）'); ylabel('Y坐标（米）');
legend('Location','best'); axis equal; grid on;

fprintf('\n========== 最终优化结果汇总 ==========\n');
fprintf('优化方法: %s\n', optimization_method);
fprintf('1. 历史最优综合得分：%.4f\n', best_fit);
fprintf('2. 最终总能耗：%.2f J（剩余：%.2f J）\n', energy_consumption(actual_iter), sum(final_E_remaining));
fprintf('3. Pareto归档非支配解数量：%d 个\n', length(pareto_archive));
fprintf('-------------------------------------\n');

if ~isempty(pareto_archive)
    arch_util = [pareto_archive.Utility];
    arch_lat = [pareto_archive.Latency];
    arch_energy = [pareto_archive.Energy];

    norm_util = (max(arch_util) - arch_util) / (max(arch_util) - min(arch_util) + 1e-6);
    norm_lat = (arch_lat - min(arch_lat)) / (max(arch_lat) - min(arch_lat) + 1e-6);
    norm_eng = (arch_energy - min(arch_energy)) / (max(arch_energy) - min(arch_energy) + 1e-6);

    W_U = 0.70;
    W_L = 0.15;
    W_E = 0.15;
    distances = sqrt(W_U * norm_util.^2 + W_L * norm_lat.^2 + W_E * norm_eng.^2);
    [~, knee_idx] = min(distances);

    [~, max_u_idx] = max(arch_util);
    [~, min_e_idx] = min(arch_energy);

    fprintf('🏆 最优效用方案 (激进救援)：\n');
    fprintf('   效用: %.1f | 时延: %.2f s | 能耗: %.2f J\n', arch_util(max_u_idx), arch_lat(max_u_idx), arch_energy(max_u_idx));

    fprintf('🔋 最省电方案 (保守续航)：\n');
    fprintf('   效用: %.1f | 时延: %.2f s | 能耗: %.2f J\n', arch_util(min_e_idx), arch_lat(min_e_idx), arch_energy(min_e_idx));

    fprintf('🌟 最佳折中方案 (Knee Point / 推荐执行)：\n');
    fprintf('   效用: %.1f | 时延: %.2f s | 能耗: %.2f J\n', arch_util(knee_idx), arch_lat(knee_idx), arch_energy(knee_idx));
end
fprintf('=====================================\n');

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
