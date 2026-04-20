function coverage_metrics = calculate_coverage_metrics(UAV_pos, User, priorities, params, RRH)
    if nargin < 6
        RRH = [];
    end

    N_User = size(User, 1);
    cover_radius = params.cover_radius;
    RRH_radius = params.RRH_radius;

    covered_by_uav = false(N_User, 1);
    covered_by_rrh = false(N_User, 1);
    covered_count = 0;
    high_priority_covered = 0;
    total_priority_sum = 0;
    covered_priority_sum = 0;

    for i = 1:N_User
        dists_uav = sqrt(sum((UAV_pos - User(i, :)).^2, 2));
        covered_by_uav(i) = any(dists_uav <= cover_radius);

        if ~isempty(RRH)
            dists_rrh = sqrt(sum((RRH - User(i, :)).^2, 2));
            covered_by_rrh(i) = any(dists_rrh <= RRH_radius);
        end

        if covered_by_uav(i) || covered_by_rrh(i)
            covered_count = covered_count + 1;
            covered_priority_sum = covered_priority_sum + priorities(i);

            if priorities(i) >= 3
                high_priority_covered = high_priority_covered + 1;
            end
        end

        total_priority_sum = total_priority_sum + priorities(i);
    end

    coverage_metrics.total_coverage_ratio = covered_count / N_User;
    coverage_metrics.high_priority_coverage_ratio = high_priority_covered / sum(priorities >= 3);
    coverage_metrics.priority_coverage_ratio = covered_priority_sum / total_priority_sum;
    coverage_metrics.covered_user_count = covered_count;
    coverage_metrics.high_priority_covered_count = high_priority_covered;
end

function energy_metrics = calculate_energy_metrics(UAV_pos, k_move, center_point)
    if nargin < 3
        center_point = [500, 500];
    end

    fly_distances = sqrt(sum((UAV_pos - center_point).^2, 2));
    total_energy = k_move * sum(fly_distances);
    avg_energy_per_uav = total_energy / size(UAV_pos, 1);

    energy_metrics.total_energy = total_energy;
    energy_metrics.avg_energy_per_uav = avg_energy_per_uav;
    energy_metrics.total_fly_distance = sum(fly_distances);
    energy_metrics.avg_fly_distance = mean(fly_distances);
    energy_metrics.max_fly_distance = max(fly_distances);
end

function convergence_metrics = calculate_convergence_metrics(fitness_curve, threshold)
    if nargin < 2
        threshold = 0.95;
    end

    [best_fitness, best_idx] = max(fitness_curve);
    target_fitness = best_fitness * threshold;

    convergence_idx = find(fitness_curve >= target_fitness, 1, 'first');

    if isempty(convergence_idx)
        convergence_metrics.convergence_FES = length(fitness_curve);
        convergence_metrics.converged = false;
    else
        convergence_metrics.convergence_FES = convergence_idx;
        convergence_metrics.converged = true;
    end

    convergence_metrics.final_fitness = best_fitness;
    convergence_metrics.improvement_rate = (best_fitness - fitness_curve(1)) / fitness_curve(1);
end
