function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = GWO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    pop_size = 40;
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
        pop_size = 40;
    end

    if ~isfield(params, 'enable_bilevel'); params.enable_bilevel = true; end
    if ~isfield(params, 'B_total'); params.B_total = 20e6; end
    if ~isfield(params, 'F_total'); params.F_total = 10e9; end
    if ~isfield(params, 'max_latency'); params.max_latency = 1.0; end
    if ~isfield(params, 'kappa'); params.kappa = 1e-27; end
    if ~isfield(params, 'P_tx'); params.P_tx = 1; end
    if ~isfield(params, 'noise'); params.noise = 1e-13; end
    if ~isfield(params, 'RRH_radius'); params.RRH_radius = 120; end
    params.RRH = RRH;

    n_vars = N_UAV * 2;

    population = zeros(pop_size, n_vars);

    center_point = [500, 500];
    jitter = 10;

    for i = 1:pop_size
        for j = 1:N_UAV
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            population(i, (j-1)*2+1) = max(Lb(1), min(Ub(1), init_x));
            population(i, (j-1)*2+2) = max(Lb(2), min(Ub(2), init_y));
        end
    end

    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});
    fitness_pop = zeros(pop_size, 1);

    for i = 1:pop_size
        uav_pos = reshape(population(i, :), N_UAV, 2);
        if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
            uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
            population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
        end
        center_point = repmat([500, 500], N_UAV, 1);
        fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
        fly_energy = params.k_move * fly_dist;
        E_curr = max(0, params.E_max - fly_energy);

        [fitness, cand_util, cand_lat, cand_nrg] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
            params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
        fitness_pop(i) = fitness;

        pareto_archive = updateParetoArchive3D(pareto_archive, uav_pos, cand_util, cand_lat, cand_nrg);
    end

    [best_fit, best_idx] = max(fitness_pop);
    bestUAV = reshape(population(best_idx, :), N_UAV, 2);
    alpha = population(best_idx, :);
    alpha_score = best_fit;

    [~, beta_idx] = max(fitness_pop);
    beta = population(beta_idx, :);
    beta_score = fitness_pop(beta_idx);

    [~, delta_idx] = max(fitness_pop);
    delta = population(delta_idx, :);
    delta_score = fitness_pop(delta_idx);

    cg_curve = zeros(1, max_iter);
    cg_curve(1) = best_fit;

    for iter = 1:max_iter
        a = 2 - iter * (2 / max_iter);

        for i = 1:pop_size
            for d = 1:n_vars
                r1_alpha = rand();
                r2_alpha = rand();
                A_alpha = 2 * a * r1_alpha - a;
                C_alpha = 2 * r2_alpha;
                D_alpha = abs(C_alpha * alpha(d) - population(i, d));
                X1_alpha = alpha(d) - A_alpha * D_alpha;

                r1_beta = rand();
                r2_beta = rand();
                A_beta = 2 * a * r1_beta - a;
                C_beta = 2 * r2_beta;
                D_beta = abs(C_beta * beta(d) - population(i, d));
                X1_beta = beta(d) - A_beta * D_beta;

                r1_delta = rand();
                r2_delta = rand();
                A_delta = 2 * a * r1_delta - a;
                C_delta = 2 * r2_delta;
                D_delta = abs(C_delta * delta(d) - population(i, d));
                X1_delta = delta(d) - A_delta * D_delta;

                population(i, d) = (X1_alpha + X1_beta + X1_delta) / 3;
            end

            dim_idx = 1;
            for j = 1:N_UAV
                population(i, dim_idx) = max(Lb(1), min(Ub(1), population(i, dim_idx)));
                population(i, dim_idx+1) = max(Lb(2), min(Ub(2), population(i, dim_idx+1)));
                dim_idx = dim_idx + 2;
            end

            uav_pos = reshape(population(i, :), N_UAV, 2);
            if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
                uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
                population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
            end

            center_point = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
            fly_energy = params.k_move * fly_dist;
            E_curr = max(0, params.E_max - fly_energy);

            [fitness, cand_util, cand_lat, cand_nrg] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

            if fitness > alpha_score
                alpha_score = fitness;
                alpha = population(i, :);
            end

            if fitness > beta_score && fitness < alpha_score
                beta_score = fitness;
                beta = population(i, :);
            end

            if fitness > delta_score && fitness < beta_score && fitness < alpha_score
                delta_score = fitness;
                delta = population(i, :);
            end

            fitness_pop(i) = fitness;

            if fitness > best_fit
                best_fit = fitness;
                bestUAV = uav_pos;
            end

            pareto_archive = updateParetoArchive3D(pareto_archive, uav_pos, cand_util, cand_lat, cand_nrg);
        end

        cg_curve(iter) = best_fit;

        if mod(iter, 50) == 0
            fprintf('GWO iter %d/%d, Best fitness: %.4f\n', iter, max_iter, best_fit);
        end
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
    best_energy = sum(params.k_move * fly_dist);
end

function feasible = checkConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    feasible = true;

    for i = 1:N_UAV
        for j = i+1:N_UAV
            dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
            if dist < D_UU
                feasible = false;
                return;
            end
        end
    end

    for i = 1:N_UAV
        for j = 1:size(RRH, 1)
            dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
            if dist < D_RU
                feasible = false;
                return;
            end
        end
    end
end

function uav_pos = enforceConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    max_iter = 100;

    for iter = 1:max_iter
        changed = false;

        for i = 1:N_UAV
            for j = i+1:N_UAV
                dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
                if dist < D_UU && dist > 0
                    direction = (uav_pos(i, :) - uav_pos(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_UU - dist) * 0.5;
                    uav_pos(j, :) = uav_pos(j, :) - direction * (D_UU - dist) * 0.5;
                    changed = true;
                end
            end
        end

        for i = 1:N_UAV
            for j = 1:size(RRH, 1)
                dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
                if dist < D_RU && dist > 0
                    direction = (uav_pos(i, :) - RRH(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_RU - dist);
                    changed = true;
                end
            end
        end

        if ~changed
            break;
        end
    end
end