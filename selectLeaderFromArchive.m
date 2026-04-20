function leader_pos = selectLeaderFromArchive(pareto_archive, N_UAV)
    if isempty(pareto_archive)
        leader_pos = rand(N_UAV, 2) * 1000;
        return;
    elseif length(pareto_archive) < 3
        idx = randi(length(pareto_archive));
        leader_pos = pareto_archive(idx).UAV_pos;
        return;
    end

    N = length(pareto_archive);
    objs = zeros(N, 3);
    for i = 1:N
        objs(i, 1) = -pareto_archive(i).Utility;
        objs(i, 2) = pareto_archive(i).Latency;
        objs(i, 3) = pareto_archive(i).Energy;
    end

    crowding_distance = zeros(N, 1);
    for m = 1:3
        [sorted_obj, sort_idx] = sort(objs(:, m));
        crowding_distance(sort_idx(1)) = inf;
        crowding_distance(sort_idx(end)) = inf;

        f_min = sorted_obj(1);
        f_max = sorted_obj(end);

        if f_max == f_min
            continue;
        end

        for i = 2:(N-1)
            distance = (sorted_obj(i+1) - sorted_obj(i-1)) / (f_max - f_min);
            crowding_distance(sort_idx(i)) = crowding_distance(sort_idx(i)) + distance;
        end
    end

    tournament_size = min(3, N);
    candidates = randperm(N, tournament_size);
    best_candidate_idx = candidates(1);

    for i = 2:tournament_size
        if crowding_distance(candidates(i)) > crowding_distance(best_candidate_idx)
            best_candidate_idx = candidates(i);
        end
    end

    leader_pos = pareto_archive(best_candidate_idx).UAV_pos;
end