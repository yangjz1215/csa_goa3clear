function subpops = initSubpopulations(N_UAV, User, RRH, priorities, subpop_params, Ub, Lb, cover_radius, D_RU)
    subpops = cell(3,1);
    high_users = User(priorities>=3,:);
    all_users = User;

    scene_size = Ub - Lb;
    D_UU = 10;

    start_point = (Ub + Lb) / 2;

    g1_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    subpops{1}.mu = g1_pos;
    subpops{1}.sigma = repmat(subpop_params.sigma0(1) * 0.3, N_UAV, 2);
    subpops{1}.prev_mu = subpops{1}.mu;
    subpops{1}.id = 1;

    g2_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    subpops{2}.mu = g2_pos;
    subpops{2}.sigma = repmat(subpop_params.sigma0(2) * 0.3, N_UAV, 2);
    subpops{2}.prev_mu = subpops{2}.mu;
    subpops{2}.id = 2;

    if size(RRH,1) > 0
        ref_rrh_idx = randi(size(RRH, 1));
        ref_rrh = RRH(ref_rrh_idx, :);
        g3_pos = generateCenteredUAV(N_UAV, ref_rrh, RRH, D_UU, D_RU, Ub, Lb, 2);
    else
        g3_pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, 1);
    end
    subpops{3}.mu = g3_pos;
    subpops{3}.sigma = repmat(subpop_params.sigma0(3) * 0.3, N_UAV, 2);
    subpops{3}.prev_mu = subpops{3}.mu;
    subpops{3}.id = 3;
end

function pos = generateCenteredUAV(N_UAV, start_point, RRH, D_UU, D_RU, Ub, Lb, jitter)
    pos = repmat(start_point, N_UAV, 1);
    max_attempts = 100;
    for i = 1:N_UAV
        valid = false;
        attempts = 0;
        while ~valid && attempts < max_attempts
            candidate = start_point + jitter * randn(1, 2);
            candidate = max(Lb, min(Ub, candidate));
            if i == 1
                valid = true;
            else
                dist_to_others = sqrt(sum((pos(1:i-1,:) - repmat(candidate, i-1, 1)).^2, 2));
                dist_to_rrh = sqrt(sum((RRH - candidate).^2, 2));
                valid = all(dist_to_others >= D_UU) && all(dist_to_rrh >= D_RU);
            end
            attempts = attempts + 1;
        end
        if valid
            pos(i,:) = candidate;
        else
            pos(i,:) = start_point;
        end
    end
end

% 辅助函数：生成分散的位置
function pos = generateDispersedPositions(N, Lb, Ub, min_dist)
    pos = zeros(N, 2);
    max_attempts = 1000;
    
    for i = 1:N
        attempts = 0;
        valid = false;
        while ~valid && attempts < max_attempts
            candidate = Lb + (Ub - Lb) .* rand(1, 2);
            if i == 1
                valid = true;
            else
                dists = sqrt(sum((pos(1:i-1,:) - candidate).^2, 2));
                valid = all(dists >= min_dist);
            end
            attempts = attempts + 1;
        end
        if ~valid
            % 如果找不到，使用网格点
            grid_size = ceil(sqrt(N));
            x_idx = mod(i-1, grid_size);
            y_idx = floor((i-1) / grid_size);
            candidate = Lb + (Ub - Lb) .* [x_idx/(grid_size-1), y_idx/(grid_size-1)];
        end
        pos(i,:) = candidate;
    end
end

% 辅助函数：在场景角落初始化
function pos = initializeAtCorners(N, Lb, Ub, min_dist)
    pos = zeros(N, 2);
    corners = [Lb; [Ub(1), Lb(2)]; [Lb(1), Ub(2)]; Ub];  % 四个角落
    
    for i = 1:N
        corner_idx = randi(4);
        corner = corners(corner_idx, :);
        offset = 200 * (rand(1,2) - 0.5);
        candidate = corner + offset;
        candidate = max(Lb, min(Ub, candidate));
        pos(i,:) = candidate;
    end
    pos = enforceDispersion(pos, min_dist, Ub, Lb);
end

% 辅助函数：强制分散处理
function pos = enforceDispersion(pos, min_dist, Ub, Lb)
    N = size(pos, 1);
    max_iter = 50;
    
    for iter = 1:max_iter
        moved = false;
        for i = 1:N
            for j = i+1:N
                dist = norm(pos(i,:) - pos(j,:));
                if dist < min_dist
                    dir = (pos(i,:) - pos(j,:)) / (dist + 1e-6);
                    delta = (min_dist - dist) / 2;
                    pos(i,:) = pos(i,:) + dir * delta;
                    pos(j,:) = pos(j,:) - dir * delta;
                    moved = true;
                end
            end
            pos(i,:) = max(Lb, min(Ub, pos(i,:)));
        end
        if ~moved
            break;
        end
    end
end