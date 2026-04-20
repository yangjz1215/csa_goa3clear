function candidates = sampleCandidates(subpop, K, N_UAV, Ub, Lb, RRH, D_UU, D_RU)
    candidates = zeros(K, N_UAV, 2);
    for i = 1:K
        uav_pos = zeros(N_UAV, 2);
        for uav_idx = 1:N_UAV
            if subpop.id == 3  % 能耗优化子种群
                pos = subpop.mu(uav_idx,:) + 0.5*subpop.sigma(uav_idx,:) .* randn(1,2);
            else  % 覆盖优化子种群
                pos = subpop.mu(uav_idx,:) + subpop.sigma(uav_idx,:) .* randn(1,2);
            end
            uav_pos(uav_idx, :) = max(Lb, min(Ub, pos));  % 边界约束
        end
        % 整合的UAV距离约束处理（原enforceUAVDistance函数）
        uav_pos = enforceUAVDistance(uav_pos, D_UU);
        candidates(i, :, :) = reshape(uav_pos, 1, N_UAV, 2);
    end
end

% 整合的辅助函数：强制UAV间最小距离（随机方向双向推开，增强版）
function pos = enforceUAVDistance(pos, D_UU)
    N = size(pos,1);
    max_iter = 20;  % 增加最大迭代次数，确保严格满足约束
    for iter = 1:max_iter
        moved = false;
        for i = 1:N
            for j = i+1:N
                dist = norm(pos(i,:) - pos(j,:));
                if dist < D_UU && dist > 0
                    % 随机方向向量（打破固定垂直方向）
                    dir = randn(1,2);  % 随机方向
                    if norm(dir) < 1e-6
                        dir = [1, 0];  % 避免零向量
                    end
                    dir = dir / norm(dir);  % 单位化
                    % 双向推开（各移动一半距离，确保至少达到最小距离）
                    delta = (D_UU * 1.1 - dist) / 2;  % 拉开到1.1倍最小距离，确保严格满足
                    pos(i,:) = pos(i,:) - dir * delta;
                    pos(j,:) = pos(j,:) + dir * delta;
                    moved = true;
                end
            end
        end
        if ~moved
            break;
        end
    end
    % 最终检查：确保所有UAV对都满足最小距离约束
    for i = 1:N
        for j = i+1:N
            dist = norm(pos(i,:) - pos(j,:));
            if dist < D_UU && dist > 0
                % 如果仍然太近，强制拉开
                dir = (pos(i,:) - pos(j,:)) / (dist + 1e-6);
                if norm(dir) < 1e-6
                    dir = [1, 0];
                end
                dir = dir / norm(dir);
                mid_point = (pos(i,:) + pos(j,:)) / 2;
                pos(i,:) = mid_point + dir * (D_UU * 1.1 / 2);
                pos(j,:) = mid_point - dir * (D_UU * 1.1 / 2);
            end
        end
    end
end