function [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, variant)
    if nargin < 13
        variant = 'proposed';
    end

    params.variant = variant;
    params.RRH = RRH;

    switch variant
        case 'proposed'
            params.enable_pareto_leader = true;
            params.enable_multi_subpop = true;
            params.enable_goa_repulsion = true;
        case 'no_pareto'
            params.enable_pareto_leader = false;
            params.enable_multi_subpop = true;
            params.enable_goa_repulsion = true;
        case 'no_subpop'
            params.enable_pareto_leader = true;
            params.enable_multi_subpop = false;
            params.enable_goa_repulsion = true;
        case 'no_goa'
            params.enable_pareto_leader = true;
            params.enable_multi_subpop = true;
            params.enable_goa_repulsion = false;
    end

    subpops = initSubpopulations(N_UAV, User, RRH, priorities, params.subpop_params, Ub, Lb, params.cover_radius, params.D_RU);

    if ~params.enable_multi_subpop
        subpops = {subpops{1}};
        params.G_weights = [1.0];
        params.subpop_params.sigma0 = params.subpop_params.sigma0(1,:);
        params.subpop_params.sigma_min = params.subpop_params.sigma_min(1,:);
        params.subpop_params.w_inertia = params.subpop_params.w_inertia(1);
        params.subpop_params.c = params.subpop_params.c(1);
        params.subpop_params.q = params.subpop_params.q(1);
        params.subpop_params.beta = params.subpop_params.beta(1);
        n_subpops = 1;
    else
        n_subpops = 3;
    end

    if isfield(params, 'K')
        params.K = round(params.K);
        params.K = max(10, min(60, params.K));
    end
    mem_matrix = cell(n_subpops,1);
    for g = 1:n_subpops
        mem_matrix{g} = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
    end

    cg_curve = zeros(1, params.FES_max);
    energy_consumption = zeros(1, params.FES_max);
    E_remaining = params.E_max * ones(N_UAV, 1);
    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});

    capturability_g = zeros(n_subpops, 1);
    for g = 1:n_subpops
        capturability_g(g) = calcCapturability(subpops{g}, 1, params.FES_max, g);
    end
    if ~isfield(params, 'RRH_radius'); params.RRH_radius = 120; end

    if ~isfield(params, 'enable_bilevel'); params.enable_bilevel = true; end
    if ~isfield(params, 'B_total'); params.B_total = 10e6; end
    if ~isfield(params, 'F_total'); params.F_total = 10e9; end
    if ~isfield(params, 'max_latency'); params.max_latency = 1.0; end
    if ~isfield(params, 'kappa'); params.kappa = 1e-27; end
    if ~isfield(params, 'P_tx'); params.P_tx = 1; end
    if ~isfield(params, 'noise'); params.noise = 1e-13; end

    [initial_fit, initial_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
        User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
    best_fit = initial_fit;
    bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
        E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

    cg_curve(1) = best_fit;
    energy_consumption(1) = initial_energy;

    for iter = 2:params.FES_max
        t = 1 - iter/params.FES_max;

        for g = 1:n_subpops
            capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
        end

        for g = 1:n_subpops
            candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
            candidates = zeros(params.K, N_UAV, 2);

            current_mem_size = size(mem_matrix{g}, 1);
            if current_mem_size < params.K
                additional_candidates = sampleCandidates(subpops{g}, params.K - current_mem_size, N_UAV, Ub, Lb, RRH, params.D_UU, params.D_RU);
                mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
            elseif current_mem_size > params.K
                mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
            end

            for i = 1:params.K
                cand_i = squeeze(candidates_init(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV; cand_i = cand_i'; end
                X_mean_g = mean(cand_i, 1);

                for uav_idx = 1:N_UAV
                    X_init = cand_i(uav_idx, :);
                    pos = X_init;
                    mem_idx = min(i, size(mem_matrix{g}, 1));
                    mem_candidate = squeeze(mem_matrix{g}(mem_idx, :, :));
                    if size(mem_candidate, 1) == 2 && size(mem_candidate, 2) == N_UAV; mem_candidate = mem_candidate'; end
                    mem_ref_pos = mem_candidate(uav_idx, :);

                    if params.enable_goa_repulsion
                        if rand >= params.subpop_params.q(g)
                            pos = goaUShape(pos, subpops{g}, mem_ref_pos, t, X_init, g);
                        else
                            pos = goaVShape(pos, subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                        end
                    else
                        pos = X_init + randn(1,2) * mean(subpops{g}.sigma(:));
                    end
                    candidates(i, uav_idx, :) = pos(:)';
                end
            end

            for i = 1:params.K
                cand_i = squeeze(candidates(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV; cand_i = cand_i'; end

                if params.enable_pareto_leader && length(pareto_archive) >= 3
                    arch_utils = [pareto_archive.Utility];
                    arch_lats = [pareto_archive.Latency];
                    arch_nrgs = [pareto_archive.Energy];

                    [~, max_u_idx] = max(arch_utils);
                    leader_G1 = reshape(pareto_archive(max_u_idx).UAV_pos, N_UAV, 2);

                    [~, min_l_idx] = min(arch_lats);
                    leader_G2 = reshape(pareto_archive(min_l_idx).UAV_pos, N_UAV, 2);

                    [~, min_e_idx] = min(arch_nrgs);
                    leader_G3 = reshape(pareto_archive(min_e_idx).UAV_pos, N_UAV, 2);

                    if g == 1; global_leader = leader_G1;
                    elseif g == 2; global_leader = leader_G2;
                    else; global_leader = leader_G3; end
                else
                    global_leader = bestUAV;
                end

                for uav_idx = 1:N_UAV
                    pos = cand_i(uav_idx, :);
                    subpop_best_uav = global_leader(uav_idx, :);
                    pos = goaTurn(pos, subpop_best_uav, capturability_g(g), t);

                    if ~params.enable_goa_repulsion
                        pos = cand_i(uav_idx, :);
                    end
                    candidates(i, uav_idx, :) = pos(:)';
                end
            end

            if ~params.enable_multi_subpop
                g_eval = 1;
            else
                g_eval = g;
            end

            mem_matrix{g} = updateMemory(mem_matrix{g}, candidates, User, priorities, ...
                E_remaining, params.E_max, params.k_move, g_eval, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end

        local_mus = zeros(n_subpops, N_UAV, 2);
        for g = 1:n_subpops
            for uav_idx = 1:N_UAV
                local_mus(g, uav_idx, :) = mean(squeeze(mem_matrix{g}(:, uav_idx, :)), 1);
            end
            subpops{g} = updateSubpopPV(subpops{g}, mem_matrix{g}, squeeze(local_mus(g,:,:)), params.subpop_params, g, iter, params.FES_max, N_UAV);
        end

        [curr_fit_best, curr_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
            User, priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

        if curr_fit_best > best_fit
            best_fit = curr_fit_best;
            bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, priorities, ...
                E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
        end

        cg_curve(iter) = best_fit;
        energy_consumption(iter) = curr_energy;

        for g = 1:n_subpops
            for i = 1:size(mem_matrix{g}, 1)
                candidate = squeeze(mem_matrix{g}(i, :, :));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                [cand_util, cand_lat, cand_nrg] = calcMEC_Objectives(candidate, User, priorities, params);
                [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, candidate, cand_util, cand_lat, cand_nrg);
            end
        end
    end
    [final_util, final_lat, final_nrg] = calcMEC_Objectives(bestUAV, User, priorities, params);
    [pareto_archive, ~] = updateParetoArchive3D(pareto_archive, bestUAV, final_util, final_lat, final_nrg);
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos,1)
        dists_uav = sqrt(sum((UAV_pos - User_pos(i,:)).^2, 2));
        if any(dists_uav <= UAV_radius) || (size(RRH,1) > 0 && any(sqrt(sum((RRH - User_pos(i,:)).^2, 2)) <= RRH_radius))
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos,1);
end

function new_pos = goaUShape(~, subpop, mem_ref_pos, t, X_init, g)
    A_g = (2*rand - 1) * [0.6, 0.7, 0.5];
    new_pos = X_init(:)' + 3*cos(2*pi*rand)*t*mean(subpop.sigma(:)) + A_g(g)*(mem_ref_pos(:)' - X_init(:)');
end

function new_pos = goaVShape(~, subpop, mem_ref_pos, t, X_init, X_mean, g)
    x = 2*pi*rand;
    V_x = (x < pi) * (-x/pi + 1) + (x >= pi) * (x/pi - 1);
    B_g = (2*rand - 1) * [0.5, 0.6, 0.4];
    new_pos = X_init(:)' + 3*V_x*t*mean(subpop.sigma(:)) + B_g(g)*(X_mean(:)' - X_init(:)');
end

function new_pos = goaTurn(pos, global_best_uav, cap, t)
    direction = sign(global_best_uav - pos) * [cos(randn*0.2), -sin(randn*0.2); sin(randn*0.2), cos(randn*0.2)];
    new_pos = pos(:)' + t * (cap * norm(pos - global_best_uav)) .* direction;
end