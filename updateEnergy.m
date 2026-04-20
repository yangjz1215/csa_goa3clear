% updateEnergy.m - 能耗更新函数
% 参考 Code_main0\UserOffloading.m 的能耗模型
% 能耗包括：移动能耗 + 悬停能耗 + 计算能耗
function E_remaining = updateEnergy(E_remaining, curr_pos, prev_pos, k_move, N_UAV, E_max, varargin)
    Pho = 100;
    ki = 1e-27;
    hover_time = 1;
    C_uav = 0;
    
    if nargin >= 7
        params = varargin{1};
        if isfield(params, 'Pho')
            Pho = params.Pho;
        end
        if isfield(params, 'ki')
            ki = params.ki;
        end
        if isfield(params, 'hover_time')
            hover_time = params.hover_time;
        end
        if isfield(params, 'C_uav')
            C_uav = params.C_uav;
        end
    end
    
    for uav_idx = 1:N_UAV
        move_dist = norm(curr_pos(uav_idx,:) - prev_pos(uav_idx,:));
        move_energy = k_move * move_dist;
        hover_energy = Pho * hover_time;
        compute_energy = ki * C_uav^2 * C_uav;
        max_consume = 0.05 * E_remaining(uav_idx);
        total_energy = move_energy + hover_energy + compute_energy;
        energy_cost = min(total_energy, max_consume);
        E_remaining(uav_idx) = max(0, E_remaining(uav_idx) - energy_cost);
    end
end