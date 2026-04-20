function [best_fit, bestUAV, cg_curve, energy_consumption, E_remaining_history, final_E_remaining, curr_curve, actual_iter, weighted_best_curve, pareto_archive] = ...
    cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)

if ~isfield(params, 'enable_early_stop')
    params.enable_early_stop = true;
end
if ~isfield(params, 'enable_smart_stop')
    params.enable_smart_stop = true;
end

params.RRH = RRH;

%% 初始化
% 传递覆盖半径参数给初始化函数
subpops = initSubpopulations(N_UAV, User, RRH, priorities, params.subpop_params, Ub, Lb, params.cover_radius, params.D_RU);

% 记忆矩阵初始化
% 确保K是整数
if isfield(params, 'K')
    params.K = round(params.K);
    params.K = max(20, min(60, params.K));  % 限制范围
end
mem_matrix = cell(3,1);
for g = 1:3
    mem_matrix{g} = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, ...
        params.D_UU, params.D_RU);
end

% 记录数组
cg_curve = zeros(1, params.FES_max);  % 历史最优适应度
curr_curve = zeros(1, params.FES_max);  % 当前迭代适应度
energy_consumption = zeros(1, params.FES_max);
E_remaining_history = zeros(params.FES_max, N_UAV);
E_remaining = params.E_max * ones(N_UAV, 1);
weighted_best_curve = zeros(1, params.FES_max);

% Pareto归档器初始化
pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

% 停滞计数器
stagnation_counter = zeros(3, 1);
prev_fits = zeros(3, 1);

% 提前停止条件相关变量
% 用于检测连续收敛（适应度变化在阈值内）
convergence_window = 10;  % 连续10代
convergence_threshold = 0.2;  % 适应度变化阈值
fit_history = zeros(convergence_window, 1);  % 记录最近N代的适应度

% 用于检测连续下降
decline_window = 5;  % 连续5代下降
decline_counter = 0;  % 下降计数器

% 计算初始捕获能力
capturability_g = zeros(3, 1);
for g = 1:3
    capturability_g(g) = calcCapturability(subpops{g}, 1, params.FES_max, g);
end
    if ~isfield(params, 'RRH_radius'); params.RRH_radius = 120; end

% 初始全局最优
[initial_fit, initial_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
    User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
    N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
best_fit = initial_fit;
bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
    E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
[~, weighted_initial, ~, ~] = calcFitness(bestUAV, User, priorities, ...
    E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
    N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
cg_curve(1) = best_fit;  % 历史最优适应度
curr_curve(1) = initial_fit;  % 当前适应度
weighted_best_curve(1) = weighted_initial;
energy_consumption(1) = initial_energy;
E_remaining_history(1,:) = E_remaining';

    % 计算初始覆盖率
    initial_cov_high = calcCoverageWithRRH(bestUAV, User(priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
    initial_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);
    
    % 初始化子种群适应度（用于停滞检测）
    for g = 1:3
        subpop_fits = zeros(1, size(mem_matrix{g},1));
        for i = 1:size(mem_matrix{g},1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end
        prev_fits(g) = max(subpop_fits);  % 记录初始子种群最优适应度
    end

    [init_util, init_lat, init_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);

    fprintf('初始化完成：综合适应度=%.4f | 真实效用(优先级和)=%.1f | 时延=%.2fs | 能耗=%.1f J\n', ...
    initial_fit, init_util, init_lat, init_nrg);

% MO-Aware Smart Stop 变量
mo_stagnation_counter = 0;
archive_size_history = zeros(1, params.FES_max);

%% 主迭代
for iter = 2:params.FES_max
    t = 1 - iter/params.FES_max;  % 俯冲强度系数
     
    % 计算当前迭代的捕获能力
    for g = 1:3
        capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
    end
    
    % 阶段2：探索阶段
    for g = 1:3
        % 候选坐标初步采样
        candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, ...
            params.D_UU, params.D_RU);
        candidates = zeros(params.K, N_UAV, 2);
        
        % 确保 mem_matrix 大小与 params.K 匹配（如果K被自适应调整）
        current_mem_size = size(mem_matrix{g}, 1);
        if current_mem_size < params.K
            % 如果K增加了，需要扩展mem_matrix
            % 通过复制现有解来扩展
            additional_needed = params.K - current_mem_size;
            additional_candidates = sampleCandidates(subpops{g}, additional_needed, N_UAV, Ub, Lb, RRH, ...
                params.D_UU, params.D_RU);
            mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
        elseif current_mem_size > params.K
            % 如果K减少了，截取前K个
            mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
        end
        
        % GOA俯冲行为调整
        for i = 1:params.K
            cand_i = squeeze(candidates_init(i, :, :));
            if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                cand_i = cand_i';
            end
            X_mean_g = mean(cand_i, 1);  % 1×2 行向量
            
            for uav_idx = 1:N_UAV
                X_init = cand_i(uav_idx, :);  % 直接从cand_i获取
                pos = X_init;
                
                % 从记忆矩阵中获取参考候选坐标（确保索引在范围内）
                mem_idx = min(i, size(mem_matrix{g}, 1));
                mem_candidate = squeeze(mem_matrix{g}(mem_idx, :, :));
                % 确保mem_candidate是N_UAV×2的矩阵
                if size(mem_candidate, 1) == 2 && size(mem_candidate, 2) == N_UAV
                    mem_candidate = mem_candidate';
                end
                % mem_candidate应该是N_UAV×2，这里取第uav_idx个UAV的参考位置
                mem_ref_pos = mem_candidate(uav_idx, :);
                
                % 选择U型/V型俯冲（使用外部函数）
                if rand >= params.subpop_params.q(g)
                    pos = goaUShape(pos, subpops{g}, mem_ref_pos, t, X_init, g);
                else
                    pos = goaVShape(pos, subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                end

                pos = max(Lb, min(Ub, pos));

                % 约束检查（被动拒绝方式）
                valid_pos = true;
                % 检查与RRH的间距
                for rrh_idx = 1:N_RRH
                    if norm(pos - RRH(rrh_idx,:)) < params.D_RU
                        valid_pos = false;
                        break;
                    end
                end
                % 检查与其他UAV的间距
                if valid_pos
                    for other_uav = 1:N_UAV
                        if other_uav ~= uav_idx
                            other_pos = cand_i(other_uav, :);
                            if norm(pos - other_pos) < params.D_UU
                                valid_pos = false;
                                break;
                            end
                        end
                    end
                end
                % 如果不满足约束，保留原位置
                if ~valid_pos
                    pos = X_init;
                end
                candidates(i, uav_idx, :) = pos(:)';
            end
        end
        
        % 阶段3：开发阶段
        
        % === 手术1：去中心化领导 - 各子种群找自己的Leader ===
        local_bests = zeros(3, N_UAV, 2);
        for g_idx = 1:3
            best_f_local = -inf;
            for i_idx = 1:size(mem_matrix{g_idx}, 1)
                cand_temp = squeeze(mem_matrix{g_idx}(i_idx,:,:));
                if size(cand_temp,1)==1; cand_temp = reshape(cand_temp, N_UAV, 2); end
                [f_val, ~] = calcFitness(cand_temp, User, priorities, E_remaining, params.E_max, params.k_move, g_idx, params.subpop_params, N_UAV, params.cover_radius, RRH, capturability_g(g_idx), N_RRH, RRH_type, UAV_type, params);
                if f_val > best_f_local
                    best_f_local = f_val;
                    local_bests(g_idx, :, :) = cand_temp;
                end
            end
        end
        
        % 确保 mem_matrix 大小与 params.K 匹配
        current_mem_size = size(mem_matrix{g}, 1);
        if current_mem_size < params.K
            additional_needed = params.K - current_mem_size;
            additional_candidates = sampleCandidates(subpops{g}, additional_needed, N_UAV, Ub, Lb, RRH, ...
                params.D_UU, params.D_RU);
            mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
        elseif current_mem_size > params.K
            mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
        end

        if length(pareto_archive) >= 3
            arch_utils = [pareto_archive.Utility];
            arch_lats = [pareto_archive.Latency];
            arch_nrgs = [pareto_archive.Energy];

            [~, max_u_idx] = max(arch_utils);
            leader_G1 = reshape(pareto_archive(max_u_idx).UAV_pos, N_UAV, 2);

            [~, min_l_idx] = min(arch_lats);
            leader_G2 = reshape(pareto_archive(min_l_idx).UAV_pos, N_UAV, 2);

            [~, min_e_idx] = min(arch_nrgs);
            leader_G3 = reshape(pareto_archive(min_e_idx).UAV_pos, N_UAV, 2);
        else
            leader_G1 = bestUAV; leader_G2 = bestUAV; leader_G3 = bestUAV;
        end

        for i = 1:params.K
            cand_i = squeeze(candidates(i, :, :));
            if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                cand_i = cand_i';
            end

            for uav_idx = 1:N_UAV
                pos = cand_i(uav_idx, :);

                if g == 1
                    subpop_best_uav = leader_G1(uav_idx, :);
                elseif g == 2
                    subpop_best_uav = leader_G2(uav_idx, :);
                else
                    subpop_best_uav = leader_G3(uav_idx, :);
                end

                pos = goaTurn(pos, subpop_best_uav, capturability_g(g), t);
                pos = max(Lb, min(Ub, pos));
                behavior_types{i} = 'turn';

                % 约束检查（被动拒绝方式）
                valid_pos = true;
                % 检查与RRH的间距
                for rrh_idx = 1:N_RRH
                    if norm(pos - RRH(rrh_idx,:)) < params.D_RU
                        valid_pos = false;
                        break;
                    end
                end
                % 检查与其他UAV的间距
                if valid_pos
                    for other_uav = 1:N_UAV
                        if other_uav ~= uav_idx
                            other_pos = cand_i(other_uav, :);
                            if norm(pos - other_pos) < params.D_UU
                                valid_pos = false;
                                break;
                            end
                        end
                    end
                end
                % 如果不满足约束，保留原位置
                if ~valid_pos
                    pos = cand_i(uav_idx, :);
                end
                candidates(i, uav_idx, :) = pos(:)';
            end
            
            % 开发阶段后，再次检查整个候选解的UAV间距约束
            candidate_pos = squeeze(candidates(i, :, :));
            if size(candidate_pos, 1) == 2 && size(candidate_pos, 2) == N_UAV
                candidate_pos = candidate_pos';
            end
            valid_candidate = true;
            for uav_a = 1:N_UAV
                for uav_b = uav_a+1:N_UAV
                    if norm(candidate_pos(uav_a,:) - candidate_pos(uav_b,:)) < params.D_UU
                        valid_candidate = false;
                        break;
                    end
                end
                if ~valid_candidate
                    break;
                end
            end
            % 如果不满足全局UAV间距约束，使用原位置
            if ~valid_candidate
                cand_init_i = squeeze(candidates_init(i, :, :));
                if size(cand_init_i, 1) == 2 && size(cand_init_i, 2) == N_UAV
                    cand_init_i = cand_init_i';
                end
                candidate_pos = cand_init_i;
            end
            candidates(i, :, :) = reshape(candidate_pos, 1, N_UAV, 2);
        end
        
        % 更新记忆矩阵
        mem_matrix{g} = updateMemory(mem_matrix{g}, candidates, User, priorities, ...
            E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
    end
    
    % 阶段4：PV矩阵更新
    local_mus = zeros(3, N_UAV, 2);
    for g = 1:3
        for uav_idx = 1:N_UAV
            local_mus(g, uav_idx, :) = mean(squeeze(mem_matrix{g}(:, uav_idx, :)), 1);
        end
    end
    % 更新各子种群PV矩阵
    for g = 1:3
        subpops{g} = updateSubpopPV(subpops{g}, mem_matrix{g}, squeeze(local_mus(g,:,:)), ...
            params.subpop_params, g, iter, params.FES_max, N_UAV);
    end
    
    % 阶段5：适应度评估与全局最优更新
    % 计算当前迭代的平均适应度（基于当前记忆矩阵的平均值，更真实反映当前状态）
    curr_fit_mean = 0;
    total_candidates = 0;
    for g = 1:3
        for i = 1:size(mem_matrix{g},1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [fit_val, ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
            curr_fit_mean = curr_fit_mean + fit_val;
            total_candidates = total_candidates + 1;
        end
    end
    curr_fit_mean = curr_fit_mean / total_candidates;
    
    % 同时计算最优适应度（用于更新全局最优）
    [curr_fit_best, curr_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
        User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
    
    % 记录当前适应度（使用最优值而非均值，避免早熟收敛误判）
    curr_curve(iter) = curr_fit_best;
    
    % 更新历史最优（基于最优适应度）
    % 适应度现在是实际优先级和，改进阈值需要相应调整
    improvement_threshold = 1e-6;
    if curr_fit_best > best_fit + improvement_threshold
        best_fit = curr_fit_best;
        bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
            E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
        [~, weighted_best, ~, ~] = calcFitness(bestUAV, User, priorities, ...
            E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
    else
        weighted_best = weighted_best_curve(iter-1);
    end
    % 确保 best_fit 是标量
    if ~isscalar(best_fit)
        best_fit = best_fit(1);
    end
    cg_curve(iter) = best_fit;  % 历史最优适应度
    
    % 确保 weighted_best 是标量
    if ~isscalar(weighted_best)
        weighted_best = weighted_best(1);
    end
    weighted_best_curve(iter) = weighted_best;
    
    % 确保 curr_energy 是标量
    if ~isscalar(curr_energy)
        curr_energy = curr_energy(1);
    end
    energy_consumption(iter) = curr_energy;
    
    % 多样性机制：如果长时间无改进，增加探索（自适应增强探索能力）
    % 调整：更频繁的触发，更积极的探索增强
    % 注意：使用curr_fit_best而不是cg_curve，因为cg_curve(iter)还未赋值
    if iter > 10
        % 检查历史最优适应度是否改进
        if curr_fit_best > best_fit - improvement_threshold
            last_improve_iter = iter;  % 当前迭代有改进
        else
            % 查找上一次改进的迭代
            if iter > 1
                last_improve_iter = find(cg_curve(1:iter-1) > best_fit - improvement_threshold, 1, 'last');
                if isempty(last_improve_iter)
                    last_improve_iter = 1;
                end
            else
                last_improve_iter = 1;
            end
        end
        stagnation_generations = iter - last_improve_iter;
        
        % 自适应触发阈值：前期更宽松（15代），后期更严格（8代），更敏感地检测停滞
        threshold = max(8, 15 - floor(iter / 40));
        
        if stagnation_generations > threshold
            % ==========================================================
            % 论文创新点：基于未覆盖感知的精准Levy跃迁 (Smart Targeted Jump)
            % ==========================================================
            high_pri_users = User(priorities >= 3, :);
            uncovered_targets = [];

            if ~isempty(high_pri_users)
                for ui = 1:size(high_pri_users, 1)
                    dists_to_uavs = sqrt(sum((bestUAV - high_pri_users(ui, :)).^2, 2));
                    if min(dists_to_uavs) > params.cover_radius
                        uncovered_targets = [uncovered_targets; high_pri_users(ui, :)];
                    end
                end
            end

            if isempty(uncovered_targets)
                if ~isempty(high_pri_users)
                    uncovered_targets = high_pri_users(randperm(size(high_pri_users, 1), min(3, size(high_pri_users, 1))), :);
                else
                    uncovered_targets = Lb + (Ub - Lb) .* rand(3, 2);
                end
            end

            for g = 1:3
                boost_factor = 1.0 + 0.5 * min(1.0, stagnation_generations / 20);
                subpops{g}.sigma = min(subpops{g}.sigma * boost_factor, 10 * mean(params.subpop_params.sigma0(:)));

                if rand < 0.5
                    perturb_count = max(1, round(size(mem_matrix{g}, 1) * 0.15));
                    perturb_idx = randperm(size(mem_matrix{g}, 1), perturb_count);

                    for idx = perturb_idx
                        temp_pos = reshape(mem_matrix{g}(idx, :, :), N_UAV, 2);

                        uavs_to_jump = randperm(N_UAV, randi([1, 2]));
                        for u = uavs_to_jump
                            u_fly_dist = norm(temp_pos(u,:) - [500, 500]);
                            real_E_remaining_u = params.E_max - params.k_move * u_fly_dist;
                            if real_E_remaining_u > 0.6 * params.E_max
                                target_c = uncovered_targets(randi(size(uncovered_targets, 1)), :);
                                temp_pos(u, :) = target_c + 15 * randn(1, 2);
                            end
                        end

                        temp_pos = max(Lb, min(Ub, temp_pos));
                        mem_matrix{g}(idx, :, :) = reshape(temp_pos, 1, N_UAV, 2);
                    end
                end
            end
        end
    end
    
    % 计算各子种群适应度
    curr_subpop_fits = zeros(3, 1);
    for g = 1:3
        subpop_fits = zeros(1, size(mem_matrix{g},1));
        for i = 1:size(mem_matrix{g},1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [subpop_fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end
        curr_subpop_fits(g) = max(subpop_fits);
    end
    
    % 注意：自适应参数调整将在覆盖率计算之后进行（见下方）
    
    % 更新停滞计数器和触发精英迁移
    % 适应度现在是实际优先级和，改进阈值需要相应调整
    % priorities_sum 已在前面计算，直接重用
    improvement_threshold_subpop = 1e-6;

    for g = 1:3
        if curr_subpop_fits(g) > prev_fits(g) + improvement_threshold_subpop
            % 有改进，重置停滞计数器
            stagnation_counter(g) = 0;
            prev_fits(g) = curr_subpop_fits(g);
        else
            % 无改进，增加停滞计数器
            stagnation_counter(g) = stagnation_counter(g) + 1;
        end
        
        % 如果子种群连续20代无改进，触发精英迁移
        if iter > 20 && stagnation_counter(g) >= 20
            fprintf('[精英迁移] 迭代 %d: 子种群 G%d 连续 %d 代无改进，触发精英迁移\n', ...
                iter, g, stagnation_counter(g));
            mem_matrix{g} = migrateElite(mem_matrix, g, Ub, Lb, User, priorities, ...
                E_remaining, params.E_max, params.k_move, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
            % 重置停滞计数器（迁移后重新开始计数）
            stagnation_counter(g) = 0;
            % 重新计算该子种群的适应度（因为记忆矩阵已更新）
            subpop_fits_new = zeros(1, size(mem_matrix{g},1));
            for i = 1:size(mem_matrix{g},1)
                candidate = squeeze(mem_matrix{g}(i,:,:));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                [subpop_fits_new(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
                    E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                    N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
            end
            curr_subpop_fits(g) = max(subpop_fits_new);
            prev_fits(g) = curr_subpop_fits(g);  % 更新prev_fits
        end
    end
    
    % 更新能耗（不再随迭代扣减，保持满电状态）
    % 适应度函数内部会根据无人机部署位置计算"如果飞到那里需要耗多少电"
    E_remaining = params.E_max * ones(N_UAV, 1);
    E_remaining_history(iter,:) = E_remaining';
    
    % 记录历史位置
    for g = 1:3
        subpops{g}.prev_mu = subpops{g}.mu;
    end
    
    % 记录覆盖率（使用固定半径，不再缩放）
    % 注意：适应度计算和覆盖率显示都使用固定半径（cover_radius=150米）
    % 覆盖率应该随着sigma扩散逐渐提高：初期聚拢（低覆盖率）→ 中期扩散（覆盖率提高）→ 后期收敛（高覆盖率）
    % 适应度 = 实际优先级加权成功数（不再归一化），基于固定150米半径
    % 覆盖率 = 覆盖用户数 / 总用户数，基于固定150米半径
    % 两者应该一致：适应度考虑优先级权重，覆盖率是简单比例
    curr_cov_high = calcCoverageWithRRH(bestUAV, User(priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
    curr_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);
    
    % 自适应参数调整（如果启用）- 移到覆盖率计算之后
    if isfield(params, 'enable_adaptive') && params.enable_adaptive
        % 构建算法状态
        algorithm_state = struct();
        algorithm_state.convergence_rate = abs(curr_fit_mean - curr_curve(max(1, iter-1))) / (curr_curve(max(1, iter-1)) + 1e-6);
        algorithm_state.diversity = std(curr_subpop_fits) / (mean(curr_subpop_fits) + 1e-6);
        algorithm_state.stagnation_count = max(stagnation_counter);
        algorithm_state.fitness_improvement = curr_fit_best - best_fit;
        algorithm_state.subpop_fits = curr_subpop_fits;
        algorithm_state.coverage_ratio = curr_cov_total;
        
        % 计算搜索范围（简化版）
        all_positions = [];
        for g = 1:3
            for i = 1:min(5, size(mem_matrix{g},1))  % 只采样部分解以加快速度
                candidate = squeeze(mem_matrix{g}(i,:,:));
                if size(candidate, 1) == 1
                    candidate = reshape(candidate, N_UAV, 2);
                end
                all_positions = [all_positions; candidate];
            end
        end
        if size(all_positions, 1) > 1
            pos_range = max(all_positions) - min(all_positions);
            algorithm_state.search_range = norm(pos_range);
        else
            algorithm_state.search_range = 0;
        end
        
        % 自适应调整参数
        if isfield(params, 'enable_adaptive') && params.enable_adaptive
            safe_sigmas = cell(3, 1);
            for g_sync = 1:3
                safe_sigmas{g_sync} = subpops{g_sync}.sigma;
            end
            safe_c = params.subpop_params.c;

            params = adaptiveParams(params, algorithm_state, iter, params.FES_max);

            for g_sync = 1:3
                subpops{g_sync}.sigma = min(subpops{g_sync}.sigma, safe_sigmas{g_sync} * 1.05);
            end

            if 3 >= 3
                subpops{3}.sigma = min(subpops{3}.sigma, safe_sigmas{3});
                params.subpop_params.c(3) = min(params.subpop_params.c(3), safe_c(3));
            end

            params.subpop_params.c = min(params.subpop_params.c, safe_c * 1.1);
        end
    end
    
    % 更新适应度历史（用于收敛检测）
    fit_history = circshift(fit_history, -1);
    fit_history(end) = curr_fit_mean;

    % 显示进度（每5次迭代或前10次）
    % 计算当前迭代的剩余能量（基于bestUAV到起飞点的飞行距离）
    fly_distances = sqrt(sum((bestUAV - repmat([500, 500], N_UAV, 1)).^2, 2));
    curr_avg_energy = mean(fly_distances) * params.k_move;
    curr_avg_remaining = params.E_max - curr_avg_energy;
    if mod(iter, 5) == 0 || iter <= 10
        [best_util, best_lat, best_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);

        fprintf('迭代 %d/%d, 综合得分: %.4f | 真实效用: %.1f | 时延: %.2f s | 能耗: %.2f J | 档案解数量: %d\n', ...
            iter, params.FES_max, best_fit, best_util, best_lat, best_nrg, length(pareto_archive));
    end
    
    % 更新Pareto归档（检查记忆矩阵中的所有候选解）
    old_archive_size = length(pareto_archive);
    pareto_updated_this_iter = false;
    for g = 1:3
        for i = 1:size(mem_matrix{g}, 1)
            candidate = squeeze(mem_matrix{g}(i, :, :));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(candidate, User, priorities, params);
            [pareto_archive, is_updated] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
            if is_updated
                pareto_updated_this_iter = true;
            end
        end
    end

    % ===== MO-Aware Smart Stop (基于 Pareto 前沿扩展的智能停止) =====
    if params.enable_smart_stop && params.enable_early_stop && iter > 100
        if pareto_updated_this_iter
            mo_stagnation_counter = 0;
        else
            mo_stagnation_counter = mo_stagnation_counter + 1;
        end

        max_mo_stagnation = 25;

        if mo_stagnation_counter >= max_mo_stagnation
            fprintf('\n[MO-Aware Smart Stop] 迭代 %d/%d: 连续 %d 代未能发现更好的 Pareto 折中解\n', ...
                iter, params.FES_max, max_mo_stagnation);
            break;
        end
    end
end

% 记录实际迭代次数
actual_iter = iter;

% 计算最终剩余能量（基于历史最优位置的飞行能耗）
center_point = repmat([500, 500], N_UAV, 1);
fly_distances = sqrt(sum((bestUAV - center_point).^2, 2));
fly_energy_costs = params.k_move * fly_distances;
final_E_remaining = max(0, params.E_max - fly_energy_costs');

% 计算最终覆盖率（使用完整半径，不再缩放）
final_cov_high = calcCoverageWithRRH(bestUAV, User(priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
final_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

[final_util, final_lat, final_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
[pareto_archive, ~] = updateParetoArchive3D(pareto_archive, bestUAV, final_util, final_lat, final_nrg);

% 输出最终结果
fprintf('\n========== 算法完成 ==========\n');
if iter < params.FES_max
    fprintf('提前停止于迭代 %d/%d\n', iter, params.FES_max);
else
    fprintf('达到最大迭代次数 %d\n', params.FES_max);
end

% 计算最终总能耗（基于历史最优位置到起飞点的飞行距离）
center_point = [500, 500];
total_fly_distance = sum(sqrt(sum((bestUAV - repmat(center_point, N_UAV, 1)).^2, 2)));
total_energy_consumed = params.k_move * total_fly_distance;
avg_energy_consumed = total_energy_consumed / N_UAV;

fprintf('最终结果：\n');
fprintf('  1. 历史最优综合得分：%.4f\n', best_fit);
fprintf('  2. 最终真实系统效用：%.1f（即原理的优先级和）\n', final_util);
fprintf('  3. 最终系统总时延：%.4f s\n', final_lat);
fprintf('  4. 最终系统总能耗：%.2f J\n', final_nrg);
fprintf('  5. Pareto 存档非支配解数量：%d 个\n', length(pareto_archive));
fprintf('============================\n');

end

% 辅助函数：计算覆盖率（考虑UAV和RRH）
function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos,1)
        % 检查UAV覆盖
        dists_uav = sqrt(sum((UAV_pos - User_pos(i,:)).^2, 2));
        covered_by_uav = any(dists_uav <= UAV_radius);
        
        % 检查RRH覆盖
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

% ========== GOA行为函数（局部函数）==========
function new_pos = goaUShape(pos, subpop, mem_ref_pos, t, X_init, g)
    r2 = 2*pi*rand;
    a_coeffs = [0.6, 0.7, 0.5];
    r4 = rand;
    A_g = (2*r4 - 1) * a_coeffs(g);

    X_r = mem_ref_pos(:)';  % 确保是行向量
    X_init = X_init(:)';    % 确保是行向量

    sigma_mean = mean(subpop.sigma(:));
    new_pos = X_init + 3*cos(r2)*t*sigma_mean + A_g*(X_r - X_init);
end

function new_pos = goaVShape(pos, subpop, mem_ref_pos, t, X_init, X_mean, g)
    r3 = rand;
    x = 2*pi*r3;

    if x >= 0 && x < pi
        V_x = -x/pi + 1;
    else
        V_x = x/pi - 1;
    end

    b_coeffs = [0.5, 0.6, 0.4];
    r5 = rand;
    B_g = (2*r5 - 1) * b_coeffs(g);

    sigma_mean = mean(subpop.sigma(:));
    X_init = X_init(:)';    % 确保是行向量
    X_mean = X_mean(:)';    % 确保是行向量
    new_pos = X_init + 3*V_x*t*sigma_mean + B_g*(X_mean - X_init);
end

function new_pos = goaTurn(pos, global_best_uav, cap, t)
    delta = cap * norm(pos - global_best_uav);
    direction = sign(global_best_uav - pos);
    theta = randn * 0.2;
    rot_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    direction = direction * rot_matrix;
    new_pos = pos(:)' + t*delta.*direction;
end