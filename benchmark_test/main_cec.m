function [best_fit, best_x, cg_curve] = main_cec(func_num, Dim, max_fes, Lb, Ub)
    N_UAV = Dim / 2;
    lamda = 10;
    Np = 300;

    mu = zeros(N_UAV, 2);
    sicma = lamda * ones(N_UAV, 2);

    population = Lb(1) + (Ub(1) - Lb(1)) * rand(N_UAV, 2);
    best_fit = benchmark_wrapper(population, func_num);
    bestUAV = population;
    best_x = population(:)';

    cg_curve = zeros(1, max_fes);
    cg_curve(1) = best_fit;

    fes_count = 1;
    stagnation_counter = 0;
    prev_best = best_fit;

    while fes_count < max_fes
        a = 2 - fes_count * (2 / max_fes);
        new_population = Lb(1) + (Ub(1) - Lb(1)) * rand(N_UAV, 2);

        for i = 1:N_UAV
            for j = 1:2
                r1 = 2 * a * rand() - a;
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                new_population(i, j) = new_population(i, j) + ...
                    (r1 * sin(r2) * (r3 * bestUAV(i, j) - new_population(i, j)));
            end

            curr_Ub = Ub(1:2);
            curr_Lb = Lb(1:2);

            flagub = new_population(i, :) > curr_Ub;
            new_population(i, flagub) = 2 * curr_Ub(flagub) - new_population(i, flagub);
            new_population(i, :) = max(curr_Lb, min(curr_Ub, new_population(i, :)));

            flaglb = new_population(i, :) < curr_Lb;
            new_population(i, flaglb) = 2 * curr_Lb(flaglb) - new_population(i, flaglb);
            new_population(i, :) = max(curr_Lb, min(curr_Ub, new_population(i, :)));
        end

        for i = 1:N_UAV
            tmp_UAV = bestUAV;
            tmp_UAV(i, :) = new_population(i, :);

            tmp_fit = benchmark_wrapper(tmp_UAV, func_num);
            fes_count = fes_count + 1;

            if tmp_fit < best_fit
                winner = 2 * (new_population(i, :) - curr_Lb) ./ (curr_Ub - curr_Lb) - 1;
                loser = 2 * (bestUAV(i, :) - curr_Lb) ./ (curr_Ub - curr_Lb) - 1;

                bestUAV(i, :) = new_population(i, :);
                best_fit = tmp_fit;
                best_x = tmp_UAV(:)';

                for k = 1:2
                    mut = mu(i, k);
                    mu(i, k) = mut + (1 / Np) * (winner(k) - loser(k));
                    tt = sicma(i, k)^2 + mut^2 - mu(i, k)^2 + (1 / Np) * (winner(k)^2 - loser(k)^2);
                    sicma(i, k) = sqrt(max(tt, 1e-10));
                end
            end

            if fes_count > max_fes, break; end
            cg_curve(fes_count) = best_fit;
        end

        if abs(best_fit - prev_best) < 1e-12
            stagnation_counter = stagnation_counter + 1;
        else
            stagnation_counter = 0;
        end
        prev_best = best_fit;
        if stagnation_counter > 200 && fes_count > max_fes * 0.8, break; end
    end
    cg_curve(fes_count+1:end) = [];
end