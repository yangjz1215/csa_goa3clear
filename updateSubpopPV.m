% updateSubpopPV.m - 子种群概率向量更新（按设计方案修改）
% 实现均值协同更新和标准差行为适配更新（加入退火与子种群差异化）
function subpop = updateSubpopPV(subpop, mem_matrix, global_best, params, g, iter, FES_max, N_UAV)
% 输入参数：
%   无 behavior_type（统一使用 'turn' 行为）

    % 1. 均值更新：分"局部更新 + 协同更新"（加入退火，限制步长，避免早期指数级提升）
    % 1.1 子种群内局部更新（基础集中趋势）
    local_mu = zeros(N_UAV, 2);
    for uav_idx = 1:N_UAV
        local_mu(uav_idx, :) = mean(squeeze(mem_matrix(:, uav_idx, :)), 1);
    end

    % 退火：降低早期的β，迭代后期逐步恢复（更平滑的退火曲线）
    % 调整：使用更平缓的曲线，让前期探索能力保持更久，避免快速上升
    progress = iter / FES_max;
    % 使用更平缓的分段曲线：前期极慢（0.2倍），中期缓慢（0.2→0.6），后期加速（0.6→1.0）
    if progress < 0.4
        anneal = 0.2 + 0.2 * (progress / 0.4);  % 0.2 → 0.4（前期极慢）
    elseif progress < 0.75
        anneal = 0.4 + 0.3 * ((progress - 0.4) / 0.35);  % 0.4 → 0.7（中期缓慢）
    else
        anneal = 0.7 + 0.3 * ((progress - 0.75) / 0.25);  % 0.7 → 1.0（后期加速）
    end
    beta_eff = params.beta(g) * anneal;

    % 目标增量
    delta = global_best - local_mu;

    % 限制每次更新的最大步长（防止早期大跳跃，后期逐步放宽）
    % 调整：前期更小的步长，避免快速上升
    max_step_base = 80;  % 基础最大步长
    % 使用更平缓的曲线：前期极小（0.3倍），后期才增大
    max_step = max_step_base * (0.3 + 0.7 * (progress^2.0));  % 前期极小，后期逐步增大
    step_norm = sqrt(sum(delta.^2, 2));
    scale = min(1, max_step ./ (step_norm + 1e-9));
    delta_capped = delta .* scale;

    subpop.prev_mu = subpop.mu;
    subpop.mu = local_mu + beta_eff * delta_capped;

    % 2. 标准差更新：分"基础衰减 + 行为适配"（子种群差异化）
    % 标准差更新策略：先增大后减小，实现从聚拢→扩散→收敛的转换
    % 初期：小sigma（初始聚拢，低覆盖率）→ 中期：大sigma（mu扩散，覆盖率提高）→ 后期：小sigma（收敛，高覆盖率）
    % 注意：sigma先增大让mu扩散提高覆盖率，然后减小实现收敛
    progress = iter / FES_max;
    
    % sigma先增大后减小：使用先增后减的曲线
    % 前期增大（0→0.5）：从初始小值逐渐增大，延长扩散阶段
    % 中期保持（0.5→0.8）：保持较大值，继续扩散
    % 后期减小（0.8→1.0）：逐渐减小，实现收敛
    % 注意：初始sigma是0.1*sigma0，growth_factor最高控制在6~8以内，避免过快扩散
    if progress < 0.5
        % 前期：从1.0逐渐增大到6.0（增长阶段拉长到50%迭代）
        growth_factor = 1.0 + 5.0 * (progress / 0.5);  % 1.0 → 6.0
    elseif progress < 0.8
        % 中期：保持较大值，轻微增大以维持探索
        growth_factor = 6.0 + 0.5 * ((progress - 0.5) / 0.3);  % 6.0 → 6.5
    else
        % 后期：逐渐减小（回落到约1.0，方便收敛）
        decay_progress = (progress - 0.8) / 0.2;  % 0 → 1
        growth_factor = 6.5 - 5.5 * (decay_progress^1.2);  % 6.5 → 1.0
    end
    
    % 根据子种群调整
    if g == 2
        % G2需要更强的探索能力，允许更大的growth_factor（上限约8）
        if progress < 0.5
            growth_factor = 1.0 + 7.0 * (progress / 0.5);  % 1.0 → 8.0
        elseif progress < 0.8
            growth_factor = 8.0 + 0.5 * ((progress - 0.5) / 0.3);  % 8.0 → 8.5
        else
            decay_progress = (progress - 0.8) / 0.2;
            growth_factor = 8.5 - 7.5 * (decay_progress^1.2);  % 8.5 → 1.0
        end
    end
    
    sigma_base = subpop.sigma * growth_factor;

    % 子种群差异化的行为适配系数
    % 注意：不再使用exploration_boost，让sigma自然衰减
    % G1：偏收敛；G2：偏外扩；G3：中性微收敛
    if g == 2
        adapt = 0.95;    % G2-转向也尽量保持外扩能力
    elseif g == 1
        adapt = 0.85;    % G1-转向收敛
    else
        adapt = 0.90;    % G3-适度收敛
    end

    subpop.sigma = sigma_base * adapt;

    % 约束sigma范围，避免过小/过大
    % 调整：提高sigma下限，防止过度收敛；根据进度动态调整下限
    min_sigma = 0.15 + 0.05 * (1 - progress);  % 前期0.2，后期0.15，保持一定探索能力
    max_sigma = 5 * mean(params.sigma0(:));  % 上限防炸裂
    subpop.sigma = max(subpop.sigma, min_sigma);
    subpop.sigma = min(subpop.sigma, max_sigma);
end