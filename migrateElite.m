% migrateElite.m - 精英迁移函数（按设计方案修改）
% 当子种群连续20轮适应度无提升时，注入其他子种群的精英解
function new_mem = migrateElite(mem_matrix, target_g, Ub, Lb, User, priorities, ...
    E_remaining, E_max, k_move, params, N_UAV, cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params_full)
% 按设计方案实现三种重组策略
% 输入参数：
%   params: subpop_params（子种群参数）
%   params_full: 完整的params结构体（包含所有参数）

    new_mem = mem_matrix{target_g};
    K = size(new_mem, 1);
    
    % 按设计方案的重组策略
    if target_g == 1
        % G1注入G2的全局探索经验
        % X_migrate^{1,t+1} = X_elite^{2,t} + 0.1·randn(1,2)·(Ub - Lb)
        src_g = 2;
        elite_count = max(1, round(0.1 * K));  % 前10%精英解
        [~, elite_indices] = getEliteSolutions(mem_matrix{src_g}, User, priorities, ...
            E_remaining, E_max, k_move, src_g, params, N_UAV, cover_radius, RRH, capturability_g(src_g), N_RRH, RRH_type, UAV_type, params_full);
        elites = mem_matrix{src_g}(elite_indices(1:elite_count), :, :);
        
    elseif target_g == 2
        % G2注入G1/G3的精英解
        % X_migrate^{2,t+1} = 0.5·X_elite^{1,t} + 0.5·X_elite^{3,t} + 0.1·randn(1,2)·(Ub - Lb)
        elite_count = max(1, round(0.1 * K));
        
        % 获取G1精英解
        [~, elite_indices_1] = getEliteSolutions(mem_matrix{1}, User, priorities, ...
            E_remaining, E_max, k_move, 1, params, N_UAV, cover_radius, RRH, capturability_g(1), N_RRH, RRH_type, UAV_type, params_full);
        elites_1 = mem_matrix{1}(elite_indices_1(1:elite_count), :, :);
        
        % 获取G3精英解
        [~, elite_indices_3] = getEliteSolutions(mem_matrix{3}, User, priorities, ...
            E_remaining, E_max, k_move, 3, params, N_UAV, cover_radius, RRH, capturability_g(3), N_RRH, RRH_type, UAV_type, params_full);
        elites_3 = mem_matrix{3}(elite_indices_3(1:elite_count), :, :);
        
        % 融合精英解
        elites = zeros(elite_count, N_UAV, 2);
        for i = 1:elite_count
            elites(i, :, :) = 0.5 * squeeze(elites_1(i, :, :)) + 0.5 * squeeze(elites_3(i, :, :));
        end
        
    else  % target_g == 3
        % G3注入G2的低能耗位置
        % X_migrate^{3,t+1} = X_elite,energy^{2,t} + 0.1·randn(1,2)·(Ub - Lb)
        src_g = 2;
        elite_count = max(1, round(0.1 * K));
        [~, elite_indices] = getEliteSolutionsByEnergy(mem_matrix{src_g}, E_remaining, E_max, N_UAV);
        elites = mem_matrix{src_g}(elite_indices(1:elite_count), :, :);
    end
    
    % 替换目标子种群的较差解（后10%）
    [~, worst_indices] = getEliteSolutions(new_mem, User, priorities, ...
        E_remaining, E_max, k_move, target_g, params, N_UAV, cover_radius, RRH, capturability_g(target_g), N_RRH, RRH_type, UAV_type, params_full);
    replace_indices = worst_indices(end-elite_count+1:end);
    
    for i = 1:elite_count
        elite_pos = squeeze(elites(i, :, :));
        % 添加随机扰动：0.1·randn(1,2)·(Ub - Lb)
        perturbed_elite = elite_pos + 0.1 * (Ub - Lb) .* randn(size(elite_pos));
        perturbed_elite = max(Lb, min(Ub, perturbed_elite));
        new_mem(replace_indices(i), :, :) = perturbed_elite;
    end
end

% 辅助函数：获取精英解索引（按适应度）
function [fits, sorted_indices] = getEliteSolutions(mem_matrix, User, priorities, ...
    E_remaining, E_max, k_move, g, params, N_UAV, cover_radius, RRH, cap, N_RRH, RRH_type, UAV_type, params_full)
    fits = zeros(1, size(mem_matrix, 1));
    for i = 1:size(mem_matrix, 1)
        candidate = squeeze(mem_matrix(i, :, :));
        [fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
            E_remaining, E_max, k_move, g, params, N_UAV, cover_radius, RRH, cap, N_RRH, RRH_type, UAV_type, params_full);
    end
    [~, sorted_indices] = sort(fits, 'descend');
end

% 辅助函数：按能耗获取精英解索引
function [energies, sorted_indices] = getEliteSolutionsByEnergy(mem_matrix, E_remaining, E_max, N_UAV)
    % 简化：基于剩余能量评估
    energies = E_remaining;  % 使用当前剩余能量
    [~, sorted_indices] = sort(energies, 'descend');  % 剩余能量高的优先
end