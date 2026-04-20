function [utility, latency, energy, success_rate] = calcMEC_Objectives(UAV_pos, User, priorities, params)
    N_User = size(User, 1);
    N_UAV = size(UAV_pos, 1);

    utility = 0;
    latency = 0;
    comp_energy = 0;
    success_count = 0;

    B_available = ones(N_UAV, 1) * params.B_total;
    F_available = ones(N_UAV, 1) * params.F_total;
    B_relay_available = ones(N_UAV, 1) * params.B_total_relay;

    [sorted_prio, sorted_idx] = sort(priorities, 'descend');

    for i = 1:N_User
        u_idx = sorted_idx(i);
        prio = sorted_prio(i);

        task_D_user = params.D(u_idx);
        task_C_user = params.C(u_idx);

        dists = sqrt(sum((UAV_pos - User(u_idx, :)).^2, 2));
        [min_dist, best_uav] = min(dists);
        min_dist = max(min_dist, 10.0);

        if min_dist <= params.cover_radius
            snr = params.P_tx / (min_dist^2 * params.noise);
            spectral_efficiency = log2(1 + snr);

            max_t_trans = params.max_latency * 0.5;
            max_t_comp = params.max_latency * 0.5;

            req_B = task_D_user / (max_t_trans * spectral_efficiency);
            req_F = task_C_user / max_t_comp;

            if B_available(best_uav) >= req_B && F_available(best_uav) >= req_F
                B_available(best_uav) = B_available(best_uav) - req_B;
                F_available(best_uav) = F_available(best_uav) - req_F;

                utility = utility + prio;
                success_count = success_count + 1;

                t_trans = task_D_user / (req_B * spectral_efficiency);
                t_comp = task_C_user / req_F;
                latency = latency + (t_trans + t_comp);

                comp_energy = comp_energy + params.kappa * req_F^2 * task_C_user;
            else
                if isfield(params, 'RRH') && ~isempty(params.RRH) && isfield(params, 'PtxU') && isfield(params, 'B_total_relay') && isfield(params, 'f_BBU')
                    dists_rrh = sqrt(sum((params.RRH - UAV_pos(best_uav, :)).^2, 2));
                    [min_dist_rrh, best_rrh] = min(dists_rrh);
                    min_dist_rrh = max(min_dist_rrh, 10.0);

                    if min_dist_rrh <= 300
                        req_B_access = task_D_user / ((params.max_latency * 0.4) * spectral_efficiency);

                        if B_available(best_uav) >= req_B_access
                            snr_relay = params.PtxU / (min_dist_rrh^2 * params.noise);
                            spectral_eff_relay = log2(1 + snr_relay);

                            req_B_relay = task_D_user / ((params.max_latency * 0.4) * spectral_eff_relay);

                            if B_relay_available(best_uav) >= req_B_relay
                                B_available(best_uav) = B_available(best_uav) - req_B_access;
                                B_relay_available(best_uav) = B_relay_available(best_uav) - req_B_relay;

                                t_trans1 = task_D_user / (req_B_access * spectral_efficiency);
                                t_trans2 = task_D_user / (req_B_relay * spectral_eff_relay);
                                t_comp_bbu = task_C_user / params.f_BBU;
                                total_relay_latency = t_trans1 + t_trans2 + t_comp_bbu;

                                if total_relay_latency <= params.max_latency
                                    utility = utility + prio;
                                    latency = latency + total_relay_latency;
                                    comp_energy = comp_energy + (params.PtxU * t_trans2);
                                    success_count = success_count + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    fly_dist = sqrt(sum((UAV_pos - repmat([500, 500], N_UAV, 1)).^2, 2));
    flight_energy = sum(fly_dist * params.k_move);

    energy = flight_energy + comp_energy;
    success_rate = success_count / N_User;
end