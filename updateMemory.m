function mem = updateMemory(old_mem, new_candidates, User, priorities, E_remaining, E_max, ...
    k_move, g, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params)
    K_old = size(old_mem, 1);
    K_new = size(new_candidates, 1);
    all_candidates = zeros(K_old + K_new, N_UAV, 2);
    all_candidates(1:K_old, :, :) = old_mem;
    all_candidates(K_old+1:end, :, :) = new_candidates;

    if nargin < 13 || isempty(capturability)
        capturability = 0.5;
    end

    if isfield(params, 'variant') && strcmp(params.variant, 'no_subpop')
        g_eval = 1;
    else
        g_eval = g;
    end

    fits = zeros(1, K_old + K_new);
    for i = 1:(K_old + K_new)
        candidate = squeeze(all_candidates(i, :, :));
        [fits(i), ~, ~, ~] = calcFitness(candidate, User, priorities, ...
            E_remaining, E_max, k_move, g_eval, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params);
    end

    [~, idx] = sort(fits, 'descend');
    mem = all_candidates(idx(1:K_old), :, :);
end