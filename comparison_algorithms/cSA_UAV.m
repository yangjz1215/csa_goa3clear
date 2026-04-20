function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
    end

    params.RRH = RRH;

    lamda = 10;
    Np = 300;

    mu = zeros(N_UAV, 2);
    sicma = lamda * ones(N_UAV, 2);

    center_point = [500, 500];
    jitter = 10;

    population = zeros(N_UAV, 2);
    for j = 1:N_UAV
        valid = false;
        attempts = 0;
        while ~valid && attempts < 100
            attempts = attempts + 1;
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            pos = [max(Lb(1), min(Ub(1), init_x)), max(Lb(2), min(Ub(2), init_y))];

            dist_uu_ok = true;
            for k = 1:j-1
                if norm(pos - population(k,:)) < params.D_UU
                    dist_uu_ok = false; break;
                end
            end
            dist_ru_ok = true;
            for k = 1:size(RRH,1)
                if norm(pos - RRH(k,:)) < params.D_RU
                    dist_ru_ok = false; break;
                end
            end

            if dist_uu_ok && dist_ru_ok
                population(j,:) = pos;
                valid = true;
            end
        end
        if ~valid
            population(j,:) = pos;
        end
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((population - center_point).^2, 2));
    fly_energy = params.k_move * fly_dist;
    E_curr = max(0, params.E_max - fly_energy);

    [best_fit, best_util, best_lat, best_nrg] = calcFitness(population, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
        params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
    bestUAV = population;

    cg_curve = zeros(1, max_iter);
    cg_curve(1) = best_fit;

    pareto_archive = struct('UAV_pos', {}, 'Utility', {}, 'Latency', {}, 'Energy', {});
    pareto_archive = updateParetoArchive3D(pareto_archive, bestUAV, best_util, best_lat, best_nrg);

    iter_count = 0;

    while iter_count < max_iter
        a = 2 - iter_count * (2 / max_iter);

        new_population = zeros(N_UAV, 2);
        for j = 1:N_UAV
            new_population(j, 1) = Lb(1) + rand() * (Ub(1) - Lb(1));
            new_population(j, 2) = Lb(2) + rand() * (Ub(2) - Lb(2));
        end

        for i = 1:N_UAV
            for j = 1:2
                r1 = 2 * a * rand() - a;
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                new_population(i, j) = new_population(i, j) + ...
                    (r1 * sin(r2) * (r3 * bestUAV(i, j) - new_population(i, j)));
            end

            flagub = new_population(i, :) > Ub;
            new_population(i, flagub) = 2 * Ub(flagub) - new_population(i, flagub);
            flaglb = new_population(i, :) < Lb;
            new_population(i, flaglb) = 2 * Lb(flaglb) - new_population(i, flaglb);
        end

        for i = 1:N_UAV
            tmp_UAV = bestUAV;
            tmp_UAV(i, :) = new_population(i, :);
            uavindex = 1:N_UAV;
            uavindex(i) = [];

            if all(sqrt(sum((tmp_UAV(i, :) - tmp_UAV(uavindex, :)).^2, 2)) >= params.D_UU) && ...
               all(sqrt(sum((tmp_UAV(i, :) - RRH(1:end, :)).^2, 2)) >= params.D_RU)

                iter_count = iter_count + 1;

                fly_dist = sqrt(sum((tmp_UAV - center_point).^2, 2));
                fly_energy = params.k_move * fly_dist;
                E_curr = max(0, params.E_max - fly_energy);

                [tmp_fitness, tmp_util, tmp_lat, tmp_nrg] = calcFitness(tmp_UAV, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                    params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

                winner = 2 * (bestUAV(i, :) - Lb) ./ (Ub - Lb) - 1;
                loser = 2 * (new_population(i, :) - Lb) ./ (Ub - Lb) - 1;

                if tmp_fitness > best_fit
                    winner = 2 * (new_population(i, :) - Lb) ./ (Ub - Lb) - 1;
                    loser = 2 * (bestUAV(i, :) - Lb) ./ (Ub - Lb) - 1;
                    bestUAV(i, :) = new_population(i, :);
                    best_fit = tmp_fitness;
                end

                for k = 1:2
                    mut = mu(i, k);
                    mu(i, k) = mut + (1 / Np) * (winner(k) - loser(k));
                    tt = sicma(i, k)^2 + mut^2 - mu(i, k)^2 + (1 / Np) * (winner(k)^2 - loser(k)^2);
                    if tt > 0
                        sicma(i, k) = sqrt(tt);
                    else
                        sicma(i, k) = 10;
                    end
                end

                if iter_count > max_iter
                    break;
                end

                cg_curve(iter_count) = best_fit;

                [~, best_util, best_lat, best_nrg] = calcFitness(bestUAV, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                    params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
                pareto_archive = updateParetoArchive3D(pareto_archive, bestUAV, best_util, best_lat, best_nrg);
            end
        end

        if mod(iter_count, 50) == 0
            fprintf('cSA iter %d/%d, Best fitness: %.4f\n', iter_count, max_iter, best_fit);
        end
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
    best_energy = sum(params.k_move * fly_dist);
end
