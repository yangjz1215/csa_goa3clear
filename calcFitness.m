function [fitness, F1, F2, F3] = calcFitness(uav_pos, User, priorities, E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params)
    % =====================================================================
    % 统一桥接评价器 - 将所有单目标/多目标算法导向向量化 MEC 物理引擎
    % =====================================================================

    if ~isfield(params, 'enable_bilevel') || ~params.enable_bilevel
        error('calcFitness 桥接器要求 params.enable_bilevel = true');
    end

    [F1, F2, F3] = calcMEC_Objectives(uav_pos, User, priorities, params);

    N_User = size(User, 1);

    max_utility = sum(priorities);
    max_latency = N_User * params.max_latency;
    max_energy = E_max * N_UAV * 2;

    if isfield(params, 'test_weights')
        w1 = params.test_weights(g, 1);
        w2 = params.test_weights(g, 2);
        w3 = params.test_weights(g, 3);
    else
        w1_def = [0.6, 0.4, 0.3];
        w2_def = [0.1, 0.4, 0.2];
        w3_def = [0.3, 0.2, 0.5];
        w1 = w1_def(g);
        w2 = w2_def(g);
        w3 = w3_def(g);
    end

    norm_util = F1 / max_utility;
    norm_lat = F2 / max_latency;
    norm_energy = F3 / max_energy;

    fitness = w1 * norm_util - w2 * norm_lat - w3 * norm_energy;
end
